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
from osgeo import gdal
#--------------------------------------
plot = True #True

# declare general paths
data_path = "/home/francesc/data/glacier_centerlines/"
out_path = "~/results/glacier_centerlines/"

# open the geotiff (DEM) with rioxarray
dem_file = "Norway_DEM_sel.tif"
dem = rio.open_rasterio(os.path.join(data_path, dem_file))

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

dataset = gdal.Open(os.path.join(data_path, dem_file))
band = dataset.GetRasterBand(1)

cols = dataset.RasterXSize
rows = dataset.RasterYSize

#map pixel/line coordinates into georeferenced space
transform = dataset.GetGeoTransform()

xOrigin = transform[0] 
yOrigin = transform[3]
pixelWidth = transform[1]
pixelHeight = -transform[5]

data = band.ReadAsArray(0, 0, cols, rows)

def profile(points_yx):
    """
    Parameters
    ----------
    points_list : list with lat, lon.

    Returns
    -------
    profile distance (arbitrary units) - altitude (m)
    """
    
    # initialize vectors
    alt = np.zeros(1)
    dist = np.zeros(1)
    dumdist=0
    
    # altitude
    for point in points_yx:
        col = int((point[0] - xOrigin) / pixelWidth)
        row = int((yOrigin - point[1] ) / pixelHeight)
    
        alt = np.append(alt, data[row][col])   
    
    #remove dummy 0 in the beginning 
    alt = alt[1:]
    
    # distance along line
    # Distance between  2 points
 
    #repeat the first point at the end
    #np.append(points_list, points_list[0])

    for i in np.arange(len(points_yx)):
        i=int(i)
        a=shpg.Point(points_yx[i])
        #last point
        if i == len(points_yx)-1:
            d = a.distance(shpg.Point(points_yx[0]))
        else:
            d = a.distance(shpg.Point(points_yx[i+1]))
        dumdist = dumdist + d
        dist = np.append(dist, dumdist)   
      
    #remove the dummy 0 ini point
    dist = dist[1:]
     
    return dist, alt

# Find for local maxima on the outline
#x, y = tuple2int(poly_pix.exterior.xy)
#ext_yx = tuple(reversed(poly_pix.exterior.xy))

#zoutline = topo[y[:-1], x[:-1]]  # last point is first point
    


def get_terminus_coord(ext_yx, zoutline):
    """This finds the terminus coordinate of the glacier.
    There is a special case for marine terminating glaciers/
    """

    perc = 10 #cfg.PARAMS['terminus_search_percentile']
    deltah = 20 #20 problem #50m (?) #cfg.PARAMS['terminus_search_altitude_range']

    #if gdir.is_tidewater and (perc > 0):
    if perc > 0:

        plow = np.percentile(zoutline, perc).astype(np.int64)

        # the minimum altitude in the glacier outline
        mini = np.min(zoutline)

        # indices of where in the outline the altitude is lower than the qth
        # percentile and lower than $delatah meters higher, than the minimum altitude
        ind = np.where((zoutline < plow) & (zoutline < (mini + deltah)))[0]

        # We take the middle of this area --> is that good? when we have several minima this does not hold...
        # We take the middle of this area
        try:
            ind_term = ind[np.round(len(ind) / 2.).astype(int)]
        except IndexError:
            # Sometimes the default perc is not large enough
            try:
                # Repeat
                perc *= 2
                plow = np.percentile(zoutline, perc).astype(np.int64)
                mini = np.min(zoutline)
                ind = np.where((zoutline < plow) &
                               (zoutline < (mini + deltah)))[0]
                ind_term = ind[np.round(len(ind) / 2.).astype(int)]
            except IndexError:
                # Last resort
                ind_term = np.argmin(zoutline)
    else:
        # easy: just the minimum
        ind_term = np.argmin(zoutline)
        # find coordinated from ind_term
    xterm = ext_yx[ind_term][0]
    yterm = ext_yx[ind_term][1]
        
    xyterm = shpg.Point(xterm, yterm)
        
    return xyterm, ind_term

#loop over all geometries, plot DEM
for i in np.arange(len(crop_extent)): 
    # start with one outline, and crop the DEM to the outline + a few grid point
    crp1 = crop_extent.iloc[[i]]
    ## crop with buffer (I am not sure about the units of the buffer, i assume meters)
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

#loop over all geometries, plot zoutline
    area = crop_extent.geometry[i].area
    points_yx=crop_extent.geometry.exterior[i].coords #list of X,Y coordinates
        
    if points_yx[0] == points_yx[-1]:
        points_yx = points_yx[:-1]
        
    prof = profile(points_yx)

    zoutline = prof[1]
    ext_yx=points_yx[:]

    xyterm, ind_term = get_terminus_coord(ext_yx, zoutline)

    if plot:
        #plt.plot(prof[0], prof[1]) 
        plt.plot(prof[0], zoutline, 'o-') #horizontal distance vs altitude
        plt.plot(prof[0][ind_term], zoutline[ind_term], 'r*') #terminus
    
        plt.show()    
    