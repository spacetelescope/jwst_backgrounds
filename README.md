jwst_backgrounds is a a simple program to predict the levels of background emission
in JWST observations, for use in proposal planning.

It accesses a precompiled background cache prepared by Space Telescope Science Institute. The background cache is hosted by the 
Mikulski Archive for Space Telescopes (MAST), so you need internet access to run the tool with the remote cache. It is possible to
download the full background cache to your local machine. Instructions for downloading the background cache can be found at http://archive.stsci.edu/archive_news/2017/08-Aug/index.html#article1

For a given target (RA, DEC), and wavelength, jwst_backgrounds does the following:
- Plot the spectrum of the background for that target on a given calendar day.
- Plot the total background for that target versus calendar day.
- Compute the number of days per year that the target is observable at low background,
  for a given wavelength and a selectable threshold.
- Save the retrieved background data to file.
  
This code was written by Jane Rigby (GSFC, Jane.Rigby@nasa.gov) and Klaus Pontoppidan (STScI, pontoppi@stsci.edu)
The background cache was prepared by Wayne Kinzel at STScI, and is the same as used by the JWST Exposure Time Calculator.

This software is provided as-is, with no warranty.

  
INSTALLATION

Using pip:
----------
```
pip install jwst_backgrounds
```

Note: healpy (version >= 1.10) is a required dependency, so if you don't have it pip will install it automatically. 

Note: to upgrade the JBT with pip use `pip install jwst_background --upgrade`

Using Conda
-----------
First clone the repository

```
git clone git@github.com:spacetelescope/jwst_backgrounds.git
cd jwst_backgrounds
conda create --name <env> --file requirements.txt
```

where `<env>` is the name of the environment you wish to create and requirements is the `requirements.txt` in the package directory.

Manually
----------
Clone the repository from github and install using `easy_install`.

```
git clone git@github.com:spacetelescope/jwst_backgrounds.git
cd jwst_backgrounds
easy_install .
```

   
RUNNING THE CODE:
```
python			# Start python.
from jwst_backgrounds import jbt 	# Import the background module
```

Below is an example that plots a background curve for a given RA, DEC, wavelength, threshold
```
jbt.get_background(261.6833333, -73.33222222, 2.15, thresh=1.1, \
                        plot_background=True, plot_bathtub=True, write_bathtub=True) 
```

TROUBLESHOOTING:
-----------
If matplotlib does not display the images, then try editing your ~/.matplotlib/matplotlibrc file,
and choosing a different backend:  
```
backend: MacOSX
backend: TkAgg
backend: GTKCairo
```

