""" Workflow:
This scripit does cut and plot a DEM correcponding to each polygon (glacier) 
given a shapefile and a DEM.
"""
# import libraries --------------------
import rioxarray as rio
import geopandas as gpd
import shapely.geometry as shpg
import os
import numpy as np
import matplotlib.pyplot as plt
#--------------------------------------

# declare general paths
data_path = "/home/francesc/data/glacier_centerlines/"
out_path = "~/results/glacier_centerlines/"

# open the geotiff with rioxarray
dem_file = "Norway_DEM_sel.tif"
dem = rio.open_rasterio(os.path.join(data_path, dem_file))

# open shapefile 
shape_file = "Norway_Inventory_sel/Norway_Inventory_sel.shp"
crop_extent = gpd.read_file(os.path.join(data_path, shape_file))

# view all polygons in shapefile:
crop_extent.plot()

# First check that projections are equal:
print('DEM crs: ', dem.rio.crs)
print('crop extent crs: ', crop_extent.crs)

if dem.rio.crs == crop_extent.crs:
    print("DEM and shapefile are in the same projection")
else:
    raise ValueError('Projections do not match.')

#loop over all geometries
for i in np.arange(len(crop_extent)): 
    # start with one outline, and crop the DEM to the outline + a few grid point
    crp1 = crop_extent.iloc[[i]]
    
    ## crop with buffer (I am not sure about the units of the buffer, i assume meters)
    dem_clipped = dem.rio.clip(crp1.buffer(20).apply(shpg.mapping),
                               crop_extent.crs)
    f, ax = plt.subplots(figsize=(8, 10))

    # assign some value to outside crop: e.g. 0 (default number is too large)
    dummy_val = dem_clipped.values[0][dem_clipped.values[0] < 1500].mean()
    dem_clipped.values[0][dem_clipped.values[0] > 1500] = dummy_val
    dem_clipped.plot(ax=ax)
    ax.set(title="Raster Layer Cropped to Geodataframe Extent")
    plt.show()
    
