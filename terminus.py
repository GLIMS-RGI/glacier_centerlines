""" 
Workflow:
This scripit does cut and plot a DEM correcponding to each polygon (glacier) 
given a shapefile and a DEM.
Then the terminus is found following 
## https://github.com/OGGM/oggm/blob/master/oggm/core/centerlines.py
and https://tc.copernicus.org/articles/8/503/2014/
https://github.com/OGGM/oggm/blob/447a49d7f936dae4870453d7c65bf2c6f861d0d8/oggm/core/gis.py#L798
"""
# import libraries --------------------
import rioxarray as rio
import geopandas as gpd
import shapely.geometry as shpg
import os
import numpy as np
import matplotlib.pyplot as plt
#--------------------------------------
# import functions
from functions import get_terminus_coord, profile, coordinate_change

plot = True #True

# declare general paths
data_path = "/home/francesc/data/glacier_centerlines/"
out_path = "~/results/glacier_centerlines/"

# open the geotiff (DEM) with rioxarray
dem_file = "Norway_DEM_sel.tif"
dem_path = os.path.join(data_path, dem_file)
dem = rio.open_rasterio(dem_path)

# open shapefile 
shape_file = "Norway_Inventory_sel/Norway_Inventory_sel.shp"
crop_extent = gpd.read_file(os.path.join(data_path, shape_file))

# view all polygons in shapefile:
if plot:
    crop_extent.plot()

# Sanity check. Check that projections are equal:
print('DEM crs: ', dem.rio.crs)
print('crop extent crs: ', crop_extent.crs)

if dem.rio.crs == crop_extent.crs:
    print("DEM and shapefile are in the same projection")
else:
    raise ValueError('Projections do not match.')

################################################################
# todo: smooth and fiter DEM
# I don't know exactly what do they mean with that in the paper.
#
#
#
################################################################

# Compute heads and tails
## start with the weighting function:
## https://github.com/OGGM/oggm/blob/master/oggm/core/centerlines.py

# get altitude and pixel info: 
# altitude, (xorigin, yorigin, pixelH, pixelW)
data, pix_params = coordinate_change(dem_path)

#loop over all geometries, plot DEM and zoutline
for i in np.arange(len(crop_extent)): 
    # start with one outline, and crop the DEM to the outline + a few grid point
    crp1 = crop_extent.iloc[[i]]
    # crop with buffer. Buffer in meters
    dem_clipped = dem.rio.clip(crp1.buffer(20).apply(shpg.mapping),
                               crop_extent.crs)

    # assign some value to outside crop: e.g. 0 (default number is too large)
    dummy_val = dem_clipped.values[0][dem_clipped.values[0] < 1500].mean()
    dem_clipped.values[0][dem_clipped.values[0] > 1500] = 0 #dummy_val
   
    if plot:
        f, ax = plt.subplots(figsize=(8, 10))
        dem_clipped.plot(ax=ax)
        ax.set(title="Raster Layer Cropped to Geodataframe Extent")
        plt.show()

#plot zoutline
    area = crop_extent.geometry[i].area
    points_yx=crop_extent.geometry.exterior[i].coords #list of X,Y coordinates
        
    if points_yx[0] == points_yx[-1]:
        points_yx = points_yx[:-1]
    
    #get profile under glacier outline: distance (arbitrary units) - altitude (m)  
    prof = profile(points_yx, data, pix_params)
    
    # get terminus coordinates and position in or data (index)
    xyterm, ind_term = get_terminus_coord(points_yx, prof[1]) #yx, zoutline

    if plot:
        #plt.plot(prof[0], prof[1]) 
        plt.plot(prof[0], prof[1], 'o-') #horizontal distance vs altitude
        plt.plot(prof[0][ind_term], prof[1][ind_term], 'r*') #terminus
    
        plt.show()    
    