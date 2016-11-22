from setuptools import setup

setup(name='PinMol',
      version='0.1',
      description='Find molecular beacons for live cell imaging of mRNAs',
      url='https://github.com/icatrina/PinMol',
      author='Irina Catrina',
      author_email='icatrina@hunter.cuny.edu',
      license='Hunter College',
      #packages=['PinMol'],
      install_requires=[
          'pandas<0.19', 'Biopython'],


      zip_safe=False,
      platforms=['any'])
