""" Workflow:

assume that the DEM and the outlines are in a cartesian projection (units: m)
open the geotiff with rioxarray
start with one outline, and crop the DEM to the outline + a few grid point
compute a mask of the glacier as in: https://github.com/OGGM/oggm/blob/447a49d7f936dae4870453d7c65bf2c6f861d0d8/oggm/core/gis.py#L798
use the OGGM code to compute the heads, terminus, but: do not simplify geometries as done in OGGM. I would really try to see if its possible to work and compute glacier heads and terminus in the native geometry resolution. OGGM code: https://github.com/OGGM/oggm/blob/master/oggm/core/centerlines.py
I don't think there is a need for the OGGM Centerline object. All calculations should be doable with shapely only.

The tools you will need:

rioxarray to read geotiff
geopandas to read and write geometries
shapely to do the geometrical stuff (as OGGM does)
scipy for the routing algorithm (as OGGM does)
"""
# import libraries --------------------
import rioxarray as rio
import geopandas as gpd
import shapely.geometry as shpg
import scipy as sc
import os
import numpy as np
import matplotlib.pyplot as plt
#--------------------------------------

# declare general paths
data_path = "/home/francesc/data/glacier_centerlines/"
out_path = "~/results/glacier_centerlines/"

# open the geotiff with rioxarray
tif_file = "Norway_DEM_sel.tif"
tif = rio.open_rasterio(os.path.join(data_path, tif_file))

# open shapefile 
shape_file = "Norway_Inventory_sel/Norway_Inventory_sel.shp"
crop_extent = gpd.read_file(os.path.join(data_path, shape_file))

# view shapefile:
crop_extent.plot()

# start with one outline, and crop the DEM to the outline + a few grid point
## First check that projections are equal:
print('DEM crs: ', tif.rio.crs)
print('crop extent crs: ', crop_extent.crs)

if tif.rio.crs == crop_extent.crs:
    print("DEM and shapefile are in the same projection")
else:
    raise ValueError('Projections do not match.')
    

## crop with buffer (I am not sure about the units of the buffer, i assume meters)
tif_clipped = tif.rio.clip(crop_extent.geometry.buffer(200).apply(shpg.mapping),
                                      crop_extent.crs)
f, ax = plt.subplots(figsize=(8, 10))

# assign some value to outside crop: e.g. 0 (default number is too large)
tif_clipped.values[0][tif_clipped.values[0] > 1500] = 0
tif_clipped.plot(ax=ax)
ax.set(title="Raster Layer Cropped to Geodataframe Extent")
ax.set_axis_off()
plt.show()

#-----------------
# compute mask as in https://github.com/OGGM/oggm/blob/447a49d7f936dae4870453d7c65bf2c6f861d0d8/oggm/core/gis.py#L798
from rasterio.mask import mask as riomask
import rasterio

# rename objects to fit oggm names 
dem_data = rasterio.open(os.path.join(data_path, tif_file))

geometry = crop_extent.geometry

masked_dem, _ = riomask(dem_data, [shpg.mapping(geometry)], filled=False)
glacier_mask = ~masked_dem[0, ...].mask

# --> it seems its all mask points are 0 :( 
# even though dem_data.crs = geometry.crs
