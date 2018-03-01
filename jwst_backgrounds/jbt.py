"""
This is a module to predict the background levels for JWST observations, 
for use in JWST proposal planning. 

It accesses a precompiled background cache prepared by STScI, to do the following:
- Plot the background versus calendar day.
- Compute the number of days per year that a target is observable at low background,
  for a given wavelength and a selectable threshold.

Software is provided as-is, with no warranty. Use the latest versions of APT and ETC to confirm 
the observability of any JWST targets. 
 
"""

import os
import struct
import urllib
import healpy
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d

import astropy.units as u

from jwst_backgrounds.version import __version__

# JWST Imaging / WFSS pixel units
WFC3IR_PIXEL = u.def_unit('WFC3/IR.pixel', (0.12825*u.arcsec)**2)
NIRISS_PIXEL = u.def_unit('NIRISS.pixel', (0.065*u.arcsec)**2)
NIRCAM_SW_PIXEL = u.def_unit('NIRCam.SW.pixel', (0.031*u.arcsec)**2)
NIRCAM_LW_PIXEL = u.def_unit('NIRCam.LW.pixel', (0.063*u.arcsec)**2)
MIRI_PIXEL = u.def_unit('MIRI.pixel', (0.11*u.arcsec)**2)
FLAMBDA_CGS = u.erg/u.second/u.cm**2/u.AA
FNU_CGS = u.erg/u.second/u.cm**2/u.Hz
PHOTLAM = u.photon/u.second/u.cm**2/u.AA
PHOTNU = u.photon/u.second/u.cm**2/u.Hz
# From Pandeia, JWST and HST primary area in cm
JWST_AREA = u.def_unit('JWST.Primary', 254009.0*u.cm**2)
HST_AREA = u.def_unit('HST.Primary', 38990.0*u.cm**2)

class JWSTBackground():
    '''
    Main background class. It is initialized with all background data for a specific
    position (RA, DEC). The wavelength at which the bathtub curve is calculated 
    can be updated as needed.
    
    Parameters
    ----------
    ra: float
        Right ascension in decimal degrees
    dec: float 
        Declination in decimal degrees
    wavelength: float 
        Wavelength at which the bathtub curve is calculated, in micron
    thresh: float
        the background threshold, relative to the minimum.  Default=1.1, which corresponds to <5% above the minimum background noise.
        Note that the actual noise difference will be even smaller, as there are often other significant sources of noise than just the
        background (source shot noise, detector noise, etc.).
        
    Attributes
    ----------
    bkg_data: dict
        Contains all the data for the background for the input (RA,DEC) position 
    bathtub: 
        Contains the (RA,DEC) background information as a function of calendar day, interpolated at wavelength
    '''
    def __init__(self, ra=53.122751, dec=-27.805089, wavelength=2., thresh=1.1):
        # global attributes
        self.cache_url = 'https://archive.stsci.edu/missions/jwst/simulations/straylight/sl_cache/' # Path to the online location of the background cache
        self.local_path = os.path.join(os.path.dirname(__file__),'refdata')
        self.wave_file = 'std_spectrum_wavelengths.txt' # The wavelength grid of the background cache
        self.thermal_file = 'thermal_curve_jwst_jrigby_cchen_1.1a.csv' # The constant (not time variable) thermal self-emission curve
        self.nside = 128  # Healpy parameter, from generate_backgroundmodel_cache.c .  
        self.wave_array,self.thermal_bg = self.read_static_data()
        self.sl_nwave = self.wave_array.size  # Size of wavelength array
        
        # input parameters
        self.ra = ra
        self.dec = dec
        self.wavelength = wavelength
        self.thresh = thresh

        # Load variable content
        self.cache_file = self.myfile_from_healpix(ra, dec)
        self.bkg_data = self.read_bkg_data(self.cache_file)
                
        # Interpolate bathtub curve and package it    
        self.make_bathtub(wavelength)


    def myfile_from_healpix(self, ra, dec):
        # old versions of healpy don't have lonlat
        healpix = healpy.pixelfunc.ang2pix(self.nside, ra, dec, nest=False, lonlat=True)  
        # we have to pad with 0s up to 6 characters to the left to match the file name convention
        healpix_str_pad = str(healpix).zfill(6)
        file = healpix_str_pad[0:4] + "/sl_pix_" + healpix_str_pad + ".bin" 
        return file

    def read_static_data(self):
        # Standard wavelength array. 
        abs_wave_file = os.path.join(self.local_path, self.wave_file)  
        wave_array = np.loadtxt(abs_wave_file)    
                
        thermal = np.genfromtxt(os.path.join(self.local_path, self.thermal_file), delimiter=',')
        # interpolate to same wavelength_array as others.
        thermal_int = self.interpolate_spec(thermal[:, 0], thermal[:,1],  wave_array, fill=0.0)  
        
        return wave_array,thermal_int

    def read_bkg_data(self, cache_file, verbose=False):
        """
        Method for reading one JWST background file, and parsing it.    
        
        Schema of each binary file in the cache:
        ----------------------------------------
        JRR verified the schema against the source code, generate_stray_light_with_threads.c. 
        The cache uses a Healpix RING tesselation, with NSIDE=128.  Every point on the sky 
        (tesselated tile) corresponds to one binary file, whose name includes its healpix 
        pixel number, in a directory corresponding to the first 4 digits of the healpix number.  

        - double RA
        - double DEC
        - double pos[3]
        - double nonzodi_bg[SL_NWAVE]
        - int[366] date_map  : This maps dates to indices.  NOTE: There are 366 days, not 365!
        - for each day in FOR:
            - double zodi_bg[SL_NWAVE]
            - double stray_light_bg[SL_NWAVE]
        
        parameters
        ----------
        cache_file: string
        
        attributes
        ----------
        
        """

        # Read the background file via http 
        try:
            # Python 3
            # Read the background cache version
            version_file = urllib.request.urlopen(self.cache_url + 'VERSION')
            sbet_file = urllib.request.urlopen(self.cache_url + cache_file)
        except:
            # Python 2
            # Read the background cache version
            version_file = urllib.urlopen(self.cache_url + 'VERSION')
            sbet_file = urllib.urlopen(self.cache_url + cache_file)

        self.cache_version = version_file.readlines()[0].decode('utf-8')[:-1]
        sbet_data = sbet_file.read()
        
        # Unpack the constant first part
        if verbose: 
            print("File has", len(sbet_data), "bytes, which is", len(sbet_data)/8., "doubles")
            
        size_calendar = struct.calcsize("366i") # bytes, not doubles
        partA = struct.unpack(str(5 + self.sl_nwave)+'d', sbet_data[0: (5 + self.sl_nwave)*8])
        ra = partA[0]
        dec = partA[1]
        pos = partA[2:5]
        nonzodi_bg = np.array(partA[5:5+self.sl_nwave])

        # Unpack the calendar dates - the dates go from 0 to 365 days.
        date_map = np.array(struct.unpack('366i', sbet_data[(5 + self.sl_nwave)*8  : (5 + self.sl_nwave)*8 + size_calendar]))
        if verbose: 
            print("Out of", len(date_map), "days, these many are legal:", np.sum(date_map >=0))

        calendar = np.where(date_map >=0)[0]

        Ndays = len(calendar) 
        if verbose: 
            print(len(date_map), Ndays)

        # Unpack part B, the time-variable part
        zodi_bg        = np.zeros((Ndays,self.sl_nwave))
        stray_light_bg = np.zeros((Ndays,self.sl_nwave))
        perday = self.sl_nwave*2
        partB= struct.unpack(str((len(calendar))*self.sl_nwave*2)+'d', sbet_data[perday*Ndays*-8 : ])

        # The index dd in zodi_bg[dd, : ] corresponds to the calendar day lookup[dd]
        for dd in range(0, int(Ndays)):
            br1 = dd*perday
            br2 = br1 + self.sl_nwave
            br3 = br2 + self.sl_nwave
            zodi_bg[dd, ] = partB[br1 : br2]
            stray_light_bg[dd, ] = partB[br2 : br3]

        # Expand static background components to the same shape as zodi_bg
        total_bg = np.tile(nonzodi_bg + self.thermal_bg,(Ndays,1)) + stray_light_bg + zodi_bg

        # pack everything up as a dict
        return {'calendar':calendar, 'ra':ra, 'dec':dec, 'pos':pos, 'wave_array':self.wave_array, 'nonzodi_bg':nonzodi_bg, 
                'thermal_bg':self.thermal_bg, 'zodi_bg':zodi_bg, 'stray_light_bg':stray_light_bg, 'total_bg':total_bg} 

    def make_bathtub(self, wavelength):
        """
        This method interpolates a bathtub curve at a given wavelength.
        It also uses the threshold fraction ("thresh") above the minimum background, to calculate number of good days.
        
        parameters
        ----------
        wavelength: float
        
        """
        
        self.wavelength = wavelength
        wave_array = self.bkg_data['wave_array']

        # Use linear interpolation to provide the background at any wavelength
        total_thiswave = (interp1d(wave_array, self.bkg_data['total_bg'], bounds_error=True))(wavelength)
        stray_thiswave = (interp1d(wave_array, self.bkg_data['stray_light_bg'], bounds_error=True))(wavelength)
        zodi_thiswave = (interp1d(wave_array, self.bkg_data['zodi_bg'], bounds_error=True))(wavelength)
        thermal_thiswave = (interp1d(wave_array, self.bkg_data['thermal_bg'], bounds_error=True))(wavelength)
        nonzodi_thiswave = (interp1d(wave_array, self.bkg_data['nonzodi_bg'], bounds_error=True))(wavelength)
            
        themin = np.min(total_thiswave)
        good_days =  int(np.sum(total_thiswave < themin * self.thresh)*1.0)
        
        self.bathtub = {'wavelength':wavelength,'themin':themin,'good_days':good_days,
                        'total_thiswave':total_thiswave,'stray_thiswave':stray_thiswave,'zodi_thiswave':zodi_thiswave,
                        'thermal_thiswave':thermal_thiswave,'nonzodi_thiswave':nonzodi_thiswave}

    def interpolate_spec(self, wave, specin, new_wave, fill=np.nan):
         # With these settings, writes NaN to extrapolated regions
        f = interp1d(wave, specin, bounds_error=False, fill_value=fill)
        new_spec = f(new_wave)
        return new_spec
    
    def get_spectrum(self, date='2019-05-01', area_unit=u.steradian, flux_unit=u.MJy, wavelength_unit=u.micron):
        from astropy.time import Time
        
        wave_array = self.bkg_data['wave_array']
        
        t = Time(date)
        doy = int(t.yday.split(':')[1])
        if doy-1 in self.bkg_data['calendar']:
            spw = wave_array*(1*wavelength_unit).unit*u.micron.to(wavelength_unit)
            
            try:
                # Fnu
                spf_unit = (1*u.MJy).to(flux_unit)/(1*area_unit).unit*(1*area_unit.to(u.steradian))
            except:
                # Flam and others need wavelength for conversion
                spf_unit = (1*u.MJy).to(u.erg/u.second/u.cm**2/u.Hz).to(flux_unit, equivalencies=u.spectral_density(spw))/(1*area_unit).unit*area_unit.to(u.steradian)
                
            spf = self.bkg_data['total_bg'][doy-1, :]*spf_unit
            
            return spw, spf
        else:
            print('Date "{0}" (yday={1}) not available'.format(date, doy))
            
            
    def plot_background(self, fontsize=16, xrange=(0.6,30), yrange=(1e-4,1e4), thisday=None):

        wave_array = self.bkg_data['wave_array']
        calendar = self.bkg_data['calendar']

        if thisday in calendar:
            thisday_index = np.where(thisday == calendar)[0][0]
        else:
            print("The input calendar day {}".format(thisday)+" is not available")
            return
                
        plt.plot(wave_array, self.bkg_data['nonzodi_bg'], label="ISM")
        plt.plot(wave_array, self.bkg_data['zodi_bg'][thisday_index, :], label="Zodi")
        plt.plot(wave_array, self.bkg_data['stray_light_bg'][thisday_index, :], label="Stray light")
        plt.plot(wave_array, self.bkg_data['thermal_bg'], label = "Thermal")
        plt.plot(wave_array, self.bkg_data['total_bg'][thisday_index, :], label = "Total", color='black', lw=3)
        plt.xlim(xrange)
        plt.ylim(yrange)
        
        plt.xlabel("wavelength (micron)", fontsize=fontsize)
        plt.ylabel("Equivalent in-field radiance (MJy/sr)", fontsize=fontsize)
        plt.title("Background for calendar day "+str(thisday))
        plt.legend()
        plt.yscale('log')
        plt.show()

    def plot_bathtub(self,showthresh=True, showplot=False, showsubbkgs=False, showannotate=True, title=False, label=False, showdate=False):
        
        from astropy.time import Time
        
        bathtub = self.bathtub # local link
        
        if not label:
            label="Total " + str(bathtub['wavelength']) + " micron"
        
        calendar = self.bkg_data['calendar']
        
        if showdate:
            t = Time(2019.+calendar/365.242, format='decimalyear')
            tdata = t.datetime
        else:
            tdata = calendar
            
        plt.scatter(tdata, bathtub['total_thiswave'], s=20, label=label)
            
        if showannotate:
            annotation = str(bathtub['good_days']) + " good days out of " + str(calendar.size) + \
                         " days observable, for threshold " + str(self.thresh)
            plt.title(annotation)
            plt.ylabel("bkg at " + str(bathtub['wavelength']) + " um (MJy/sr)", fontsize=12)
        else: 
            plt.ylabel("bkg (MJy/SR)", fontsize=fontsize)

        if showsubbkgs:
            plt.scatter(tdata, bathtub['zodi_thiswave'], s=20, label="Zodiacal")
            plt.scatter(tdata, bathtub['stray_thiswave'], s=20, label="Stray light")
            plt.scatter(tdata, bathtub['nonzodi_thiswave']*np.ones_like(calendar), s=20, label="ISM+CIB")
            plt.scatter(tdata, bathtub['thermal_thiswave']*np.ones_like(calendar), s=20, label="Thermal")
            plt.legend(fontsize=10, frameon=False, labelspacing=0)
            plt.grid()
            plt.locator_params(axis='x', nbins=10)
            plt.locator_params(axis='y', nbins=10)

        if showthresh: 
            percentiles = (bathtub['themin'], bathtub['themin']*self.thresh)
            plt.hlines(percentiles, tdata[0], tdata[-1], color='black')
                    
        if title: 
            plt.title(title)
        
        if showdate:
            plt.xlabel("Date", fontsize=12)
            plt.gcf().autofmt_xdate()  # orient date labels at a slant
            plt.xlim(tdata[0],tdata[-1])
        else:
            plt.xlabel("Day of the year", fontsize=12)
            plt.xlim(0,366)
            
        plt.show()
        
    def write_bathtub(self,bathtub_file='background_versus_day.txt'):
        f = open(bathtub_file,'w')
        header_text = ["# Output of JWST_backgrounds version " + str(__version__) + "\n",
                       "# background cache version " + str(self.cache_version) + '\n',
                       "\n"
                       "# for RA="+str(self.ra) + ", DEC=" + str(self.dec) + " at wavelength=" + str(self.wavelength) + " micron \n",
                       "# Columns: \n",
                       "# - Calendar day (Jan1=0) \n",
                       "# - Total background (MJy/sr)\n"] 
        for line in header_text:
            f.write(line)               
        
        for i,calendar_day in enumerate(self.bkg_data['calendar']):
            f.write('{0}    {1:5.4f}'.format(calendar_day, self.bathtub['total_thiswave'][i])+'\n')
        
        f.close()

    def write_background(self,background_file='background.txt', thisday=None):
        calendar = self.bkg_data['calendar']

        if thisday in calendar:
            thisday_index = np.where(thisday == calendar)[0][0]
        else:
            print("The input calendar day {}".format(thisday)+" is not available")
            return

        f = open(background_file,'w')
        header_text = ["# Output of JWST_backgrounds version " + str(__version__) + "\n",
                       "# background cache version " + str(self.cache_version) + '\n',
                       "\n"
                       "# for RA="+str(self.ra) + ", DEC=" + str(self.dec) + " On calendar day " + str(thisday) + "\n",
                       "# Columns: \n",
                       "# - Wavelength [micron] \n",
                       "# - Total background (MJy/sr)\n", 
                       "# - In-field zodiacal light (MJy/sr)\n",
                       "# - In-field galactic light (MJy/sr)\n",
                       "# - Stray light (MJy/sr)\n", 
                       "# - Thermal self-emission (MJy/sr)\n"]
        for line in header_text:
            f.write(line)               
        
        for i,wavelength in enumerate(self.bkg_data['wave_array']):
            f.write('{0:f}    {1:5.4f}    {2:5.4f}    {3:5.4f}    {4:5.4f}    {5:5.4f}'.format(wavelength, \
                    self.bkg_data['total_bg'][thisday_index][i],self.bkg_data['zodi_bg'][thisday_index][i],self.bkg_data['nonzodi_bg'][i], \
                    self.bkg_data['stray_light_bg'][thisday_index][i],self.bkg_data['thermal_bg'][i])+'\n')
        
        f.close()

            
def get_background(ra, dec, wavelength, thresh=1.1, plot_background=True, plot_bathtub=True, thisday=None,
                   showsubbkgs=False, write_background=True, write_bathtub=True, background_file='background.txt', 
                   bathtub_file='background_versus_day.txt'):
    """
    This is the main method, which serves as a wrapper to get the background data and create plots and outputs with one command.
    
    Parameters
    ----------
    ra: float
        Right ascension in decimal degrees
    dec: float 
        Declination in decimal degrees
    wavelength: float 
        Wavelength at which the bathtub curve is calculated, in micron
    thresh: float
        the background threshold, relative to the minimum.  Default=1.1, which corresponds to <5% above the minimum background noise.
        Note that the actual noise difference will be even smaller, as there are often other significant sources of noise than just the
        background (source shot noise, detector noise, etc.).
    plot_background: bool 
        whether to plot the background spectrum (and its components) for this day.
    thisday: int
        calendar day to use for plot_spec.  If not given, will use the average of visible calendar days.
    plot_days: bool
        whether to show the plot of background at wavelength_input versus calendar days
    showsubbkgs: bool
        whether to show the components of the background in the bathtub plot.
    write_bathtub: bool
        whether to print the background levels that are plotted in plot_days to an output file
    outfile:     output filename

    """    
    bkg = JWSTBackground(ra,dec,wavelength, thresh=thresh)
    calendar = bkg.bkg_data['calendar']
    
    print("These coordinates are observable by JWST", len(bkg.bkg_data['calendar']), "days per year.")
    print("For", bkg.bathtub['good_days'], "of those days, the background is <", thresh, "times the minimum, at wavelength", wavelength, "micron")
    
    # Figure out which day to use for the single-day background 
    if thisday not in calendar: 
        ndays = calendar.size
        if ndays>0:
            thisday_input = thisday
            thisday = calendar[int(ndays/2)] # plot the middle of the available dates in the calendar
            print("Warning: The input calendar day {}".format(thisday_input)+" is not available, assuming the middle day: {} instead".format(thisday))
        else:
            print("No valid days")
            return

    if plot_background:
        bkg.plot_background(thisday=thisday)
    
    if write_background:
        bkg.write_background(thisday=thisday, background_file=background_file)

    if plot_bathtub:
        bkg.plot_bathtub(showsubbkgs=showsubbkgs)
    
    if write_bathtub:
        bkg.write_bathtub(bathtub_file=bathtub_file)
        
        
