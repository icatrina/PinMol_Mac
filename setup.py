from setuptools import setup

setup(name='PinMol',
      version='1.0-beta',
      description='Find molecular beacons for live cell imaging of mRNAs',
      url='https://github.com/icatrina/PinMol_Mac',
      author='Irina E. Catrina',
      author_email='iecatrina@gmail.com',
      license='Hunter College',
      #packages=['PinMol'],
      install_requires=[
          'pandas<0.19', 'Biopython'],


      zip_safe=False,
      platforms=['any'])
