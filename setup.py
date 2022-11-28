from os import path
from setuptools import setup

from jwst_backgrounds.version import __version__

here = path.abspath(path.dirname(__file__))

# Get the long description from the README file
with open(path.join(here, 'README.md'), encoding='utf-8') as f:
    long_description = f.read()

setup(name='jwst_backgrounds',
      version=__version__,
      description='Retrieve and plot JWST background information',
      long_description=long_description,
      long_description_content_type='text/markdown',

      # The project's main homepage.
      url='https://github.com/spacetelescope/jwst_background',

      # Author details
      author='Jane Rigby (GSFC) and Klaus Pontoppidan (STScI)',
      author_email='Jane.Rigby@nasa.gov',

      # Choose your license
      license='BSD',
      download_url = 'https://github.com/spacetelescope/jwst_backgrounds/',
      packages=['jwst_backgrounds'],
      package_data={'jwst_backgrounds': ['refdata/*.csv','refdata/*.txt','agreement/*.pdf']},
      install_requires=['healpy>=1.10',
                        'matplotlib>=3.1.1',
                        'numpy>=1.17.0',
                        'scipy>=1.1.0'],
      entry_points = {'console_scripts': ['jwst_backgrounds=jwst_backgrounds.cli:main']}
      )
