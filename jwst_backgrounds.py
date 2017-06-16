''' This is a simple program to predict the background levels for JWST observations, 
for use in proposal planning.
It accesses a precompiled background cache prepared by STScI, to do the following:
- Plot the background versus calendar day.
- Compute the number of days per year that a target is observable at low background,
  for a given wavelength and a selectable threshold.
This code was written by Jane Rigby, Jane.Rigby@nasa.gov
The background cache was prepared by Wane Kinzel at STScI.  
Software is provided as-is, with no warranty.

 Here is the schema for the precompiled background cache.
 (JRR verified the schema against the source code, generate_stray_light_with_threads.c)
 The cache uses a Healpix RING tesselation, with NSIDE=128.  Every point on the sky 
 (tesselated tile) corresponds to one binary file, whose name includes its healpix 
 pixel number, in a directory corresponding to the first 4 digits of the healpix number.  

Here's the schema of each binary file in the cache:
double RA
double DEC
double pos[3]
double nonzodi_bg[SL_NWAVE]
int[366] date_map  # This maps dates to indices.  **There are 366 days, not 365!**
for each day in FOR:
  double zodi_bg[SL_NWAVE]
  double stray_light_bg[SL_NWAVE]


Note: This is a subset of Rigby's longer code JRR_Code/read_JWST_bkg.py .  This is
just the part that a user needs to make a background-versus-calendar "bathtub" curve.
'''

#import jrr  # Moved these routines to this file, for portability.
import glob
import subprocess
from os.path import basename
import struct
import re
import healpy
import numpy as np
import pandas
from scipy.interpolate import interp1d
from matplotlib import pyplot as plt
#import seaborn as sns
from matplotlib.backends.backend_pdf import PdfPages

#### User should change these paths  ############################################
bkg_dir  = "/Volumes/Apps_and_Docs/MISSIONS/JWST/Zody_bathtubs/sl_cache.v1.0/"  # on satchmo
#################################################################################

#### User should not touch these #################################################
nside = 128  # Healpy parameter, from generate_backgroundmodel_cache.c .  
##################################################################################



# Routines to make this all work
###########  Copied over from jrr ####################################
def convert_RADEC_Galactic(RA_deg, DEC_deg) :    # Convert (decimal) RA, DEC to Galactic.  Can be LISTS []
    thisradec = SkyCoord(RA_deg, DEC_deg, unit=(units.deg, units.deg), frame='icrs')
    return(thisradec.galactic.l.value, thisradec.galactic.b.value)

def convert_RADEC_Ecliptic(RA_deg, DEC_deg) :    # Convert (decimal) RA, DEC to Ecliptic. Can be LISTS []
    thisradec = SkyCoord(RA_deg, DEC_deg, unit=(units.deg, units.deg), frame='icrs')
    return(thisradec.barycentrictrueecliptic.lon.value, thisradec.barycentrictrueecliptic.lat.value)

def convert_RADEC_GalEclip_df(df, colra='RA', coldec='DEC') :   # copied over from jrr.util
    # For a dataframe, compute Galactic & Ecliptic coords from RADEC
    (tempL, tempB) = convert_RADEC_Galactic(df[colra], df[coldec])
    (templon, templat) = convert_RADEC_Ecliptic(df[colra], df[coldec])
    df['Gal_lon'] = tempL
    df['Gal_lat'] = tempB
    df['Ecl_lon'] = templon
    df['Ecl_lat'] = templat
    return(0)  # acts on the dataframe

def put_header_on_file(infile, header_text, outfile) :
    ''' Pandas doesn't allow user to add explanatory header
    to output files.  Headers are useful.  So, writing wrapper to add one.'''
    tmp = "/tmp/header"
    with open(tmp, "w") as myfile:  myfile.write(header_text)
    subprocess.check_output("cat " + tmp + " " + infile + " > " + outfile, shell=True)
    return(0)
######################################################################

def find_nearest(array,value):
    idx = (np.abs(array-value)).argmin()
    return array[idx]

def rebin_spec_new(wave, specin, new_wave, fill=np.nan):
    f = interp1d(wave, specin, bounds_error=False, fill_value=fill)  # With these settings, writes NaN to extrapolated regions
    new_spec = f(new_wave)
    return(new_spec)

def read_JWST_precompiled_bkg(infile, bkg_dir, showplot=False, thisday=-99, verbose=False) :
    # Reads one JWST background file, in binary format.
    wave_file = bkg_dir + "updated_std_spectrum_wavelengths.txt"  # Standard wavelength array.  Should be SL_NWave=108 long
    wave_array = np.loadtxt(wave_file)    
    SL_NWAVE = len(wave_array)  # should be 108.  Size of wavelength array
    thermal_file = "thermal_curve_jwst_jrigby_1.1.csv"  # The constant (not time variable) thermal self-emission curve
    temp_thermal = np.genfromtxt(bkg_dir + thermal_file, delimiter=',')
    thermal = rebin_spec_new(temp_thermal[:, 0], temp_thermal[:,1],  wave_array, fill=0.0)  # rebin to same wavelength_array as others.
    
    sbet_file = open(bkg_dir + infile)
    sbet_data = sbet_file.read()
    # Unpack the constant first part
    if verbose: print "File has", len(sbet_data), "bytes, which is", len(sbet_data)/8., "doubles"
    size_calendar = struct.calcsize("366i") # bytes, not doubles
    partA = struct.unpack(str(5 + SL_NWAVE)+'d', sbet_data[0: (5 + SL_NWAVE)*8])
    RA = partA[0]
    DEC = partA[1]
    pos = partA[2:5]
    nonzodi_bg = np.array(partA[5:5+SL_NWAVE])

    # Unpack the calendar dates      # code goes from 0 to 365 days.
    date_map = np.array(struct.unpack('366i', sbet_data[(5 + SL_NWAVE)*8  : (5 + SL_NWAVE)*8 + size_calendar]))
    if verbose: print "Out of", len(date_map), "days, these many are legal:", np.sum(date_map >=0)
    #print "indices of days:", date_map[date_map>=0]
    calendar = np.where(date_map >=0)[0]
    #print "calendar date:", calendar
    # So, the index dd in zodi_bg[dd, : ]  corresponds to the calendar day lookup[dd]
    Ndays = len(calendar) 
    if verbose: print len(date_map), Ndays

    # Unpack part B, the time-variable part
    zodi_bg        = np.zeros((Ndays,SL_NWAVE))
    stray_light_bg = np.zeros((Ndays,SL_NWAVE))
    perday = SL_NWAVE*2
    partB= struct.unpack(str((len(calendar))*SL_NWAVE*2)+'d', sbet_data[perday*Ndays*-8 : ])

    for dd in range(0, int(Ndays)):
        br1 = dd*perday
        br2 = br1 + SL_NWAVE
        br3 = br2 + SL_NWAVE
        #print "Breaking at:", br1, br2, br3
        zodi_bg[dd, ]        = partB[br1 : br2]
        stray_light_bg[dd, ] = partB[br2 : br3]
    expand = np.ones((Ndays,SL_NWAVE))  # same shape as zodi_bg
    total = nonzodi_bg * expand + thermal * expand + stray_light_bg + zodi_bg

    if showplot :
        fontsize=16
        plt.clf()
        if thisday in calendar : pass
        else :   thisday = find_nearest(calendar, np.mean(calendar))  # plot the middle of the calendar
        print "Plotting spectrum for calendar day", thisday
        plt.plot(wave_array, nonzodi_bg, label="ISM")
        plt.plot(wave_array, zodi_bg[thisday, :], label="Zodi")
        plt.plot(wave_array, stray_light_bg[thisday, :], label="Stray light")
        plt.plot(wave_array, thermal, label = "Thermal")
        plt.plot(wave_array, total[thisday, :], label = "Total", color='black', lw=3)
        plt.xlim(0.6,30)
        plt.xlabel("wavelength (micron)", fontsize=fontsize)
        plt.ylabel("Equivalent in-field radiance (MJy/SR)", fontsize=fontsize)
        plt.legend()
        plt.yscale('log')
        plt.show()
    return((calendar, RA, DEC, pos, wave_array, nonzodi_bg, thermal, zodi_bg, stray_light_bg, total))  #pack it up as a tuple

def index_of_wavelength(wave_array, desired_wavelength) :  # look up index of wavelength array corresponding to desired wavelength
    the_index = np.where(wave_array == desired_wavelength)
    return(the_index[0][0])

def myfile_from_healpix(healpix) :
    return ( str(healpix)[0:4] + "/sl_pix_" + str(healpix) + ".bin")

def calc_the_healpix(df):  # For each row of a dataframe, take RA, DEC in degrees and calculate the healpix number
    df['healpix'] = df.apply(lambda row : format(healpy.pixelfunc.ang2pix(nside, row.RA_deg, row.DEC_deg, nest=False, lonlat=True), '06d'), axis=1)
    df['healpix'] = df['healpix'].astype('str')
    df['healpix_asindex'] = df['healpix'].str.lstrip('0').astype('int')
    return(0)

def make_bathtub(results, wavelength_desired, thresh, showthresh=True, showplot=False, showsubbkgs=False, showannotate=True, title=False, label=False) :
    # Once binary file was read w  read_JWST_precompiled_bkg(), compute the bathtub, and optionally, make a plot.
    # thresh is threshold above minimum background, to calculate number of good days
    (calendar, RA, DEC, pos, wave_array, nonzodi_bg, thermal, zodi_bg, stray_light_bg, total)  = results # show how to break it up
    the_index = index_of_wavelength(wave_array, wavelength_desired)
    total_thiswave = total[ :, the_index]
    stray_thiswave = stray_light_bg[ :, the_index]
    zodi_thiswave =  zodi_bg[ :, the_index]
    
    themin = np.min(total_thiswave)
    allgood =  int(np.sum(total_thiswave < themin * thresh)*1.0)
    if showplot:
        if showannotate:
            annotation = str(allgood) + " good days out of " + str(len(calendar)) + " days observable, for thresh " + str(thresh)
            #print str(wavelength_desired) + " " + annotation
            plt.annotate(annotation, (0.05,0.05), xycoords="axes fraction", fontsize=12)
        if not label : label="Total " + str(wavelength_desired) + " micron"
        plt.scatter(calendar, total_thiswave, s=20, label=label)
        if showsubbkgs :
            plt.scatter(calendar, zodi_thiswave, s=20, label="Zodiacal")
            plt.scatter(calendar, stray_thiswave, s=20, label="Stray light")
            plt.scatter(calendar, nonzodi_bg[the_index]*np.ones_like(zodi_thiswave), s=20, label="ISM+CIB")
            plt.scatter(calendar, thermal[the_index]*np.ones_like(zodi_thiswave), s=20, label="Thermal")
            plt.legend(fontsize=10, frameon=False, labelspacing=0)
            plt.grid()
            plt.locator_params(axis='x', nbins=10)
            plt.locator_params(axis='y', nbins=10)
        percentiles = (themin, themin*thresh)
        if showthresh : plt.hlines(percentiles, 0, 365, color='black')
        plt.xlabel("Day of the year", fontsize=fontsize)
        plt.xlim(0,366)
        if showannotate : plt.ylabel("bkg at " + str(wave_array[the_index]) + " um (MJy/SR)", fontsize=fontsize)
        else : plt.ylabel("bkg (MJy/SR)", fontsize=fontsize)
        if title : plt.title(title)
    return(allgood)  # Returns the number of days in the FOR with a background below 



###################################################
###  This is a wrapper to plot the background spectrum, and then plot the background versus calendar day
def plot_bkg(RA, DEC, wavelength_input, thresh=1.1, plot_spec=True, thisday=-99, plot_days=True, showsubbkgs=False, print2file=False, outfile='background_versus_day.csv') :
    # Required arguments are RA and DEC (in decimal degrees), and wavelength (in micron).
    # Optional arguments:
    #         thresh:      the background threshold, relative to the minimum.  Default=1.1, which is 10% above the min bkg.
    #         plot_spec:   [bool], whether to plot the spectrum for one day
    #         thisday:     calendar day to use for plot_spec.  If not given, will use the middle of the oberving window.
    #         plot_days:   [bool], whether to show the plot of background at wavelength_input versus calendar days
    #         showsubbkgs: [bool], whether to show the components of the background in the plot_days plot.
    #         print2file   [bool], whether to print the background levels that are plotted in plot_days to an outfile
    #         outfile:     output filename
    if not (RA and DEC and wavelength_input) : raise Exception("Required arguments for plot_bkg are RA, DEC, and wavelength.  Optional argument is threshold.")
    print "Plotting the background versus calendar day for RA, DEC, wave, thresh of:", RA, DEC, wavelength_input, thresh
    healpix = healpy.pixelfunc.ang2pix(nside, RA, DEC, nest=False, lonlat=True)  # old versions of healpy don't have lonlat
    myfile = myfile_from_healpix(healpix)   # Retrieve the name of the healpix file, including leading zero formatting

    if plot_spec:  print "Plotting the background spectrum for an example day.  Turn this off w plot_spec=False"
    results = read_JWST_precompiled_bkg(myfile, bkg_dir, showplot=plot_spec, thisday=thisday)  # Retrieve the bkg file
    (calendar, RA, DEC, pos, wave_array, nonzodi_bg, thermal, zodi_bg, stray_light_bg, total)  = results  # parse results
    wavelength_desired = find_nearest(wave_array, wavelength_input)  # Nearest neighbor interpolation of wavelength
    if wavelength_desired != wavelength_input :
        print "Using wave", wavelength_desired, "as the nearest neighbor to input", wavelength_input, "micron"
    allgood = make_bathtub(results, wavelength_desired, thresh, showplot=plot_days, showsubbkgs=showsubbkgs)  # Compute bathtub, plot it.
    if plot_days : print "Plotting background versus calendar day of the year"
    plt.show()
    print "RESULTS:  These coordinates are observable by JWST", len(calendar), "days per year."
    print "RESULTS:  For", allgood, "of those days, the background is <", thresh, "of the minimum, at wavelength", wavelength_desired, "micron"
    if print2file:  # print the total background levels
         the_index = index_of_wavelength(wave_array, wavelength_desired)
         bath = pandas.DataFrame({'calendar_day': calendar, 'total_bkg' : total[ :, the_index]}).set_index('calendar_day')
         bath.to_csv("tempfile")
         print "Just wrote total background versus calendar day to", outfile, "\nHere is a subset:"
         header_text = "#Output of JWST_backgrounds.py, plot_bkg().\n# Columns are Calendar day (Jan1=0) and total background (in MegaJanskies per sterradian)\n"
         header_text += "#  for  RA="+str(RA)+ "   DEC="+str(DEC)+"    wave="+str(wavelength_desired)+" micron\n"
         put_header_on_file("tempfile", header_text, outfile) 
         print bath.head(10)
    return(results)

## Here is an example of using this code:
#fontsize=16  ### User may want to change these
#results = plot_bkg(261.6833333, -73.33222222, 2.15, thresh=1.1, plot_spec=True, thisday=200, plot_days=True, showsubbkgs=False, print2file=True)
#(calendar, RA, DEC, pos, wave_array, nonzodi_bg, thermal, zodi_bg, stray_light_bg, total) = results  #if you want to dig into the results.
