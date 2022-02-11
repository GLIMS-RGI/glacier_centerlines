"""
#-----------------
# compute mask as in https://github.com/OGGM/oggm/blob/447a49d7f936dae4870453d7c65bf2c6f861d0d8/oggm/core/gis.py#L798
"""
import os
from rasterio.mask import mask as riomask
import rasterio
import shapely.geometry as shpg
import geopandas as gpd
import numpy as np
import rioxarray as rio
from scipy.ndimage import binary_erosion
import pandas as pd
from utils import (ncDataset, nicenumber)


import salem
#from salem.gis import transform_proj

# declare general paths
data_path = "/home/francesc/data/glacier_centerlines/"
out_path = "~/results/glacier_centerlines/"

# open the geotiff with rioxarray
tif_file = "Norway_DEM_sel.tif"
tif = rio.open_rasterio(os.path.join(data_path, tif_file))

# open shapefile 
shape_file = "Norway_Inventory_sel/Norway_Inventory_sel.shp"
crop_extent = gpd.read_file(os.path.join(data_path, shape_file))

#### modify original function:
"""Compute glacier masks based on much simpler rules than OGGM's default.
This is therefore more robust: we use this function to compute glacier
hypsometries.
Parameters
----------
gdir : :py:class:`oggm.GlacierDirectory`
    where to write the data
write_hypsometry : bool
    whether to write out the hypsometry file or not - it is used by e.g,
    rgitools
"""

# In case nominal, just raise
#if gdir.is_nominal:
#    raise GeometryError('{} is a nominal glacier.'.format(gdir.rgi_id))

#if not os.path.exists(gdir.get_filepath('gridded_data')):
#    # In a possible future, we might actually want to raise a
#    # deprecation warning here
#    process_dem(gdir)

# Geometries
#geometry = gdir.read_shapefile('outlines').geometry[0]
geometry = crop_extent.geometry[10] # warning: we are taking the first polygon only

# rio metadata
#with rasterio.open(gdir.get_filepath('dem'), 'r', driver='GTiff') as ds:
#    data = ds.read(1).astype(rasterio.float32)
#    profile = ds.profile
with rasterio.open(os.path.join(data_path, tif_file)) as ds:
    data = ds.read(1).astype(rasterio.float32)
    profile = ds.profile

# simple trick to correct invalid polys:
# http://stackoverflow.com/questions/20833344/
# fix-invalid-polygon-python-shapely
geometry = geometry.buffer(0)

if not geometry.is_valid:
# This is a OGGM error
#    raise InvalidDEMError('This glacier geometry is not valid.')
    raise ValueError('This glacier geometry is not valid.')

# Compute the glacier mask using rasterio
# Small detour as mask only accepts DataReader objects
profile['dtype'] = 'int16'
profile.pop('nodata', None)
with rasterio.io.MemoryFile() as memfile:
    with memfile.open(**profile) as dataset:
        dataset.write(data.astype(np.int16)[np.newaxis, ...])
    dem_data = rasterio.open(memfile.name)
    masked_dem, _ = riomask(dem_data, [shpg.mapping(geometry)],
                            filled=False)
glacier_mask = ~masked_dem[0, ...].mask

# Same without nunataks
with rasterio.io.MemoryFile() as memfile:
    with memfile.open(**profile) as dataset:
        dataset.write(data.astype(np.int16)[np.newaxis, ...])
    dem_data = rasterio.open(memfile.name)
    poly = shpg.mapping(shpg.Polygon(geometry.exterior))
    masked_dem, _ = riomask(dem_data, [poly],
                            filled=False)
glacier_mask_nonuna = ~masked_dem[0, ...].mask

# Glacier exterior excluding nunataks
erode = binary_erosion(glacier_mask_nonuna)
glacier_ext = glacier_mask_nonuna ^ erode
glacier_ext = np.where(glacier_mask_nonuna, glacier_ext, 0)

#dem = read_geotiff_dem(gdir)
with rasterio.open(os.path.join(data_path, tif_file)) as ds:
    topo = ds.read(1).astype(rasterio.float32)
    topo[topo <= -999.] = np.NaN
    topo[ds.read_masks(1) == 0] = np.NaN
    dem = topo


# Last sanity check based on the masked dem
tmp_max = np.nanmax(dem[glacier_mask])
tmp_min = np.nanmin(dem[glacier_mask])
if tmp_max < (tmp_min + 1):
    raise ValueError("'({}) min equal max in the masked DEM.'")
    #raise InvalidDEMError('({}) min equal max in the masked DEM.'
    #                      .format(gdir.rgi_id))
# until here ok #
#-------------------------------------------------------------------
# class GriddedNcdfFile(object):
#     """Creates or opens a gridded netcdf file template.
#     The other variables have to be created and filled by the calling
#     routine.
#     """
#     def __init__(self, gdir, basename='gridded_data', reset=False):
#         self.fpath = gdir.get_filepath(basename)
#         self.grid = gdir.grid
#         if reset and os.path.exists(self.fpath):
#             os.remove(self.fpath)

#     def __enter__(self):

#         if os.path.exists(self.fpath):
#             # Already there - just append
#             self.nc = ncDataset(self.fpath, 'a', format='NETCDF4')
#             return self.nc

# Create and fill
fpath = "netCDF_out_test.nc"
nc = ncDataset(fpath, 'w', format='NETCDF4')

# dummy data test
dx = 100
dy = 100
nx = 2685
ny = 2364
ulx = 4.633e+05
uly = 7.436e+06
x0y0 = (ulx+dx/2, uly-dx/2)  # To pixel center coordinates

# try the curent crs
utm_proj = salem.check_crs(crop_extent.crs)

grid = salem.Grid(proj=utm_proj, nxny=(nx, ny), dxdy=(dx, -dx), x0y0=x0y0)
# try with dummy grid nx and grid ny
nc.dimensions['x'] = grid.nx
nc.dimensions['y'] = grid.ny
#nc.createDimension('x', grid.nx)
#nc.createDimension('y', grid.ny)

nc.author = 'OGGM'
nc.author_info = 'Open Global Glacier Model'
nc.pyproj_srs = grid.proj.srs

x = grid.x0 + np.arange(grid.nx) * grid.dx
y = grid.y0 + np.arange(grid.ny) * grid.dy

# i get problems with the "'int' object has no attribute '_dimid'", maybe because i am  not
# working with an objectof class GriddedNcdfFile 
# v = nc.createVariable('xx', 'f4', ('x',), zlib=True)
# v.units = 'm'
# v.long_name = 'x coordinate of projection'
# v.standard_name = 'projection_x_coordinate'
# v[:] = x

# v = nc.createVariable('yy', 'f4', ('y',), zlib=True)
# v.units = 'm'
# v.long_name = 'y coordinate of projection'
# v.standard_name = 'projection_y_coordinate'
# v[:] = y

    #     self.nc = nc
    #     return nc

    # def __exit__(self, exc_type, exc_value, exc_traceback):
    #     self.nc.close()

# write out the grids in the netcdf file
#with GriddedNcdfFile(gdir) as nc:

#if 'glacier_mask' not in nc.variables:
#    v = nc.createVariable('glacier_mask', 'i1', ('y', 'x', ),
#                      zlib=True)
#    v.units = '-'
#    v.long_name = 'Glacier mask'
#else:
#    v = nc.variables['glacier_mask']
#v = glacier_mask
#
#if 'glacier_ext' not in nc.variables:
#    v = nc.createVariable('glacier_ext', 'i1', ('y', 'x', ),
#                          zlib=True)
#    v.units = '-'
#    v.long_name = 'Glacier external boundaries'
#else:
#    v = nc.variables['glacier_ext']
#v = glacier_ext

## Log DEM that needed processing within the glacier mask
#if 'topo_valid_mask' not in nc.variables:
#    msg = ('You seem to be running from old preprocessed directories. '
#           'See https://github.com/OGGM/oggm/issues/1095 for a fix.')
#    raise ValueError(msg)
#    #raise InvalidWorkflowError(msg)
#valid_mask = nc.variables['topo_valid_mask'][:]
#if gdir.get_diagnostics().get('dem_needed_interpolation', False):
#    pnan = (valid_mask == 0) & glacier_mask
#    gdir.add_to_diagnostics('dem_invalid_perc_in_mask',
#                            np.sum(pnan) / np.sum(glacier_mask))

# add some meta stats and close
nc.max_h_dem = np.nanmax(dem)
nc.min_h_dem = np.nanmin(dem)
dem_on_g = dem[np.where(glacier_mask)]
nc.max_h_glacier = np.nanmax(dem_on_g)
nc.min_h_glacier = np.nanmin(dem_on_g)

# Last sanity check
if nc.max_h_glacier < (nc.min_h_glacier + 1):
#        raise InvalidDEMError('({}) min equal max in the masked DEM.'
#                              .format(gdir.rgi_id))
    raise ValueError('({}) min equal max in the masked DEM.')


# hypsometry if asked for
#if not write_hypsometry:
#    return

bsize = 50.
dem_on_ice = dem[glacier_mask]
bins = np.arange(nicenumber(dem_on_ice.min(), bsize, lower=True),
                 nicenumber(dem_on_ice.max(), bsize) + 0.01, bsize)

h, _ = np.histogram(dem_on_ice, bins)
h = h / np.sum(h) * 1000  # in permil

# We want to convert the bins to ints but preserve their sum to 1000
# Start with everything rounded down, then round up the numbers with the
# highest fractional parts until the desired sum is reached.
hi = np.floor(h).astype(int)
hup = np.ceil(h).astype(int)
aso = np.argsort(hup - h)
for i in aso:
    hi[i] = hup[i]
    if np.sum(hi) == 1000:
        break

# slope
#sy, sx = np.gradient(dem, gdir.grid.dx)
sy, sx = np.gradient(dem, grid.dx)
aspect = np.arctan2(np.mean(-sx[glacier_mask]), np.mean(sy[glacier_mask]))
aspect = np.rad2deg(aspect)
if aspect < 0:
    aspect += 360
slope = np.arctan(np.sqrt(sx ** 2 + sy ** 2))
avg_slope = np.rad2deg(np.mean(slope[glacier_mask]))

    # write
df = pd.DataFrame()
#df['RGIId'] = [gdir.rgi_id]
#df['GLIMSId'] = [gdir.glims_id]
df['Zmin'] = [dem_on_ice.min()]
df['Zmax'] = [dem_on_ice.max()]
df['Zmed'] = [np.median(dem_on_ice)]
#df['Area'] = [gdir.rgi_area_km2]
df['Slope'] = [avg_slope]
df['Aspect'] = [aspect]
for b, bs in zip(hi, (bins[1:] + bins[:-1])/2):
    df['{}'.format(np.round(bs).astype(int))] = [b]
df.to_csv("out_test.csv", index=False)