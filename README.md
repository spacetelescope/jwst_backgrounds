jwst_backgrounds is a a simple program to predict the levels of background emission
in JWST observations, for use in proposal planning.

It accesses a precompiled background cache prepared by STScI.

For a given target (RA, DEC), and wavelength, it does the following:
- Plot the spectrum of the background for that target on a given calendar day.
- Plot the total background for that target versus calendar day.
- Compute the number of days per year that the target is observable at low background,
  for a given wavelength and a selectable threshold.
  
This code was written by Jane Rigby, Jane.Rigby@nasa.gov
The background cache was prepared by Wane Kinzel at STScI.  
Software is provided as-is, with no warranty.

  
INSTALLATION:
0) jwst_backgrounds has been tested on regular Anaconda, and AstroConda 0.0.1
1) Install healpy with pip install --user healpy
2) Download jwst_backgrounds 
3) Install jwst_bacgrounds with "python setup.py install --user"
Note) healpy (version >= 1.10) is required. The setup should install it automatically, but a 
bug in the healpy setup currently prevents this from working. 
   
RUNNING THE CODE:
>python			# Start python.
from jwst_backgrounds import bg_tools 	# Import the background module

Below is an example that plots a background curve for a given RA, DEC, wavelength, threshold
bg_tools.get_background(261.6833333, -73.33222222, 2.15, thresh=1.21, plot_background=True, plot_bathtub=True, write_bathtub=True) 


TROUBLESHOOTING:
If matplotlib does not display the images, then try editing your ~/.matplotlib/matplotlibrc file,
and choosing a different backend:  
backend: MacOSX
backend: TkAgg
backend: GTKCairo

