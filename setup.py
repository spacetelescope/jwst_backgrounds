from setuptools import setup
 
setup(name='jwst_backgrounds',
      version='1.1',
      description='Retrieve and plot JWST background information',
      author='Jane Rigby (GSFC) and Klaus Pontoppidan (STScI)',
      author_email='Jane.Rigby@nasa.gov',
      url='http://jwst.stsci.edu/',
      download_url = 'https://github.com/spacetelescope/jwst_backgrounds/archive/jwst_backgrounds_1.0.tar.gz',
      packages=['jwst_backgrounds'],
      package_data={'jwst_backgrounds': ['refdata/*.csv','refdata/*.txt']},
      install_requires=['healpy>=1.10']
      )

    
