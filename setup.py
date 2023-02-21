from setuptools import setup

setup(
    name='glacier_centerlines',
    version='0.5.3',
    description='Compute glacier centerlines using RGI data and a well known algorithm (see documentetion)',
    author='Francesc Roura Adserias @OGGM developpers',
    author_email='info@oggm.org',
    packages=['glacier_centerlines'],
    install_requires=[
    	'numpy==1.22.0',
	'scipy=1.7.3',
	'salem',
	'shapely==1.8.0',
	'rioxarray==0.9.1',
	'geopandas==0.10.2',
#	'skimage',
        'scikit-image',
	'oggm==1.5.3',
	'gdal'#osgeo is not accepted' 
#	'functools'
    ],
)
