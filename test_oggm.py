""" Workflow:
This scripit does cut and plot a DEM correcponding to each polygon (glacier) 
given a shapefile and a DEM.

Using the same methodology as in 
https://github.com/OGGM/oggm/blob/447a49d7f936dae4870453d7c65bf2c6f861d0d8/oggm/core/gis.py#L798
"""
# import libraries --------------------
import rioxarray as rio
import rasterio
import geopandas as gpd
import shapely.geometry as shpg
import os
import numpy as np
import matplotlib.pyplot as plt
from rasterio.mask import mask as riomask
from scipy.ndimage import binary_erosion
from utils import (ncDataset, nicenumber)
import salem
import pandas as pd
#--------------------------------------

# declare general paths
data_path = "/home/francesc/data/glacier_centerlines/"
out_path = "~/results/glacier_centerlines/"

# open the geotiff with rioxarray
dem_file = "Norway_DEM_sel.tif"
dem_0 = rio.open_rasterio(os.path.join(data_path, dem_file))

# open shapefile 
shape_file = "Norway_Inventory_sel/Norway_Inventory_sel.shp"
crop_extent = gpd.read_file(os.path.join(data_path, shape_file))

# view all polygons in shapefile:
crop_extent.plot()

# First check that projections are equal:
print('DEM crs: ', dem_0.rio.crs)
print('crop extent crs: ', crop_extent.crs)

if dem_0.rio.crs == crop_extent.crs:
    print("DEM and shapefile are in the same projection")
else:
    raise ValueError('Projections do not match.')

# Geometries
geometry = crop_extent.geometry[10] # warning: we are taking the first polygon only

# rio metadata
with rasterio.open(os.path.join(data_path, dem_file)) as ds:
    data = ds.read(1).astype(rasterio.float32)
    profile = ds.profile
    #what is "profile??"
    
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

#plot mask
plt.imshow(glacier_mask, cmap='hot')
plt.show()    

# Same without nunataks
with rasterio.io.MemoryFile() as memfile:
    with memfile.open(**profile) as dataset:
        dataset.write(data.astype(np.int16)[np.newaxis, ...])
    dem_data = rasterio.open(memfile.name)
    poly = shpg.mapping(shpg.Polygon(geometry.exterior))
    masked_dem, _ = riomask(dem_data, [poly],
                            filled=False)
glacier_mask_nonuna = ~masked_dem[0, ...].mask
#plt.imshow(glacier_mask_nonuna, cmap='hot')
#plt.show() 

# Glacier exterior excluding nunataks
erode = binary_erosion(glacier_mask_nonuna)
glacier_ext = glacier_mask_nonuna ^ erode
glacier_ext = np.where(glacier_mask_nonuna, glacier_ext, 0)
    

with rasterio.open(os.path.join(data_path, dem_file)) as ds:
    topo = ds.read(1).astype(rasterio.float32)
    topo[topo <= -999.] = np.NaN
    topo[ds.read_masks(1) == 0] = np.NaN
    dem = topo
    
#plt.imshow(dem, cmap='hot')
#plt.show()

# Last sanity check based on the masked dem
tmp_max = np.nanmax(dem[glacier_mask])
tmp_min = np.nanmin(dem[glacier_mask])
if tmp_max < (tmp_min + 1):
    raise ValueError("'({}) min equal max in the masked DEM.'")

# Create and fill netCDF
fpath = "netCDF_out_test.nc"
nc = ncDataset(fpath, 'w', format='NETCDF4')

# create grid for nc file
nx = glacier_mask.shape[0]
ny = glacier_mask.shape[1]
dx = abs(max(dem_0.x)-min(dem_0.x))/nx
dy = abs(max(dem_0.y)-min(dem_0.y))/ny
ulx = min(dem_0.x)
uly = max(dem_0.y)
x0y0 = (ulx+dx/2, uly-dx/2)  # To pixel center coordinates

# try the curent crs
utm_proj = salem.check_crs(crop_extent.crs)
    
grid = salem.Grid(proj=utm_proj, nxny=(nx, ny), dxdy=(dx, -dx), x0y0=x0y0)

nc.dimensions['x'] = grid.nx
nc.dimensions['y'] = grid.ny

#nc.createDimension('x', grid.nx)
#nc.createDimension('y', grid.ny)

nc.author = 'OGGM'
nc.author_info = 'Open Global Glacier Model'
nc.pyproj_srs = grid.proj.srs

x = grid.x0 + np.arange(grid.nx) * grid.dx
y = grid.y0 + np.arange(grid.ny) * grid.dy

#v = nc.createVariable('glacier_mask', 'i1', ('y', 'x', ), zlib=True)
#v.units = '-'
#v.long_name = 'Glacier mask'


# add some meta stats and close
nc.max_h_dem = np.nanmax(dem)
nc.min_h_dem = np.nanmin(dem)
dem_on_g = dem[np.where(glacier_mask)]
nc.max_h_glacier = np.nanmax(dem_on_g)
nc.min_h_glacier = np.nanmin(dem_on_g)

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

nc.close()
