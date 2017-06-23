from setuptools import setup
 
setup(name='jwst_backgrounds',
      version='1.0',
      description='Retrieve and plot JWST background information',
      author='Jane Rigby',
      author_email='Jane.Rigby@nasa.gov',
      url='http://jwst.stsci.edu/',
      packages=['jwst_backgrounds'],
      package_data={'jwst_backgrounds': ['refdata/*.csv','refdata/*.txt']},
      install_requires=['healpy>=1.10']
      )

    
