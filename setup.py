from setuptools import setup

from jwst_backgrounds.version import __version__

setup(name='jwst_backgrounds',
      version=__version__,
      description='Retrieve and plot JWST background information',
      author='Jane Rigby (GSFC) and Klaus Pontoppidan (STScI)',
      author_email='Jane.Rigby@nasa.gov',
      url='http://jwst.stsci.edu/',
      download_url = 'https://github.com/spacetelescope/jwst_backgrounds/',
      packages=['jwst_backgrounds'],
      package_data={'jwst_backgrounds': ['refdata/*.csv','refdata/*.txt','agreement/*.pdf']},
      install_requires=['healpy>=1.10',
                        'matplotlib==3.1.1',
                        'numpy==1.17.0',
                        'scipy==1.1.0'],
      entry_points = {'console_scripts': ['jwst_backgrounds=jwst_backgrounds.cli:main']}
      )

    
