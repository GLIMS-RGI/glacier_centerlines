""" Workflow:
This scripit does cut and plot a DEM correcponding to each polygon (glacier) 
given a shapefile and a DEM.
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
    dem_clipped.values[0][dem_clipped.values[0] > 1500] = 0 #dummy_val
    dem_clipped.plot(ax=ax)
    ax.set(title="Raster Layer Cropped to Geodataframe Extent")
    plt.show()
 
# Compute heads and tails
## start with the weighting function:
## https://github.com/OGGM/oggm/blob/master/oggm/core/centerlines.py
crop_extent.geometry.exterior[0].coords[0]

from osgeo import gdal

driver = gdal.GetDriverByName('GTiff')
#filename = "/home/zeito/pyqgis_data/aleatorio.tif" #path to raster
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

def profile(points_list):
    """

    Parameters
    ----------
    points_list : list (?) with latlon.

    Returns
    -------
    profile distance-altitude

    """
    alt = np.zeros(1)
    dist = np.zeros(1)
    dumdist=0
    
    
    for point in points_list:
        col = int((point[0] - xOrigin) / pixelWidth)
        row = int((yOrigin - point[1] ) / pixelHeight)
    
        alt = np.append(alt, data[row][col])   
    
    # distance along line
    # Distance between  2 points
    #a=shpg.Point(points_list[0])
    #a.distance(shpg.Point(points_list[1]))
    
    np.append(points_list, points_list[0])
    for i in np.arange(len(points_list)):
        i=int(i)
        a=shpg.Point(points_list[i])
        if i == len(points_list)-1:
            d = a.distance(shpg.Point(points_list[0]))
        else:
            d = a.distance(shpg.Point(points_list[i+1]))
        dumdist = dumdist + d
        dist = np.append(dist, dumdist)   
    
    #remove dummy 0 in the begining
    alt = alt[1:]
    
    #remove repeated last value 
    #alt = alt[1:]
    
    dist = dist[1:]
     
    return dist, alt

fig = list()
for i in np.arange(len(crop_extent)):
    points_list=crop_extent.geometry.exterior[i].coords #list of X,Y coordinates
    prof = profile(points_list)

    plt.plot(prof[0], prof[1])
    plt.show()
    
    
# Find for local maxima on the outline
#x, y = tuple2int(poly_pix.exterior.xy)
#ext_yx = tuple(reversed(poly_pix.exterior.xy))

#zoutline = topo[y[:-1], x[:-1]]  # last point is first point
    
zoutline = dem_clipped.copy()

def _get_terminus_coord(gdir, ext_yx, zoutline):
    """This finds the terminus coordinate of the glacier.
    There is a special case for marine terminating glaciers/
    """

    perc = 20 #cfg.PARAMS['terminus_search_percentile']
    deltah = 50 #50m (?) #cfg.PARAMS['terminus_search_altitude_range']

    if 0 == 0: #gdir.is_tidewater and (perc > 0):
        # There is calving

        # find the lowest percentile
        inglac = zoutline.values[zoutline.values != 0]
        inglacx = np.repeat(zoutline.x, len(zoutline.y)) 
        inglacy = np.repeat(zoutline.y, len(zoutline.x))

        plow = np.percentile(inglac, perc).astype(np.int64)

        # the minimum altitude in the glacier
        mini = np.min(inglac)

        # indices of where in the outline the altitude is lower than the qth
        # percentile and lower than $delatah meters higher, than the minimum altitude
        ind = np.where((inglac < plow) & (inglac < (mini + deltah)))[0]

        # We take the middle of this area --> is that good? when we have several minima this does not hold...
        #try:
        ind_term = ind[np.round(len(ind) / 2.).astype(int)]
        #except IndexError:
        #x = np.array(zoutline.x)
        #y = np.array(zoutline.y)
        # find coordinated from ind_term
        xterm = inglacx.values[ind_term]
        yterm = inglacy.values[ind_term]
        xyterm = shpg.Point(xterm, yterm)
        return xyterm
    
    
#        except IndexError:
#            # Sometimes the default perc is not large enough
#            try:
#                # Repeat
#                perc *= 2
#                plow = np.percentile(zoutline, perc).astype(np.int64)
#                mini = np.min(zoutline)
#                ind = np.where((zoutline < plow) &
#                               (zoutline < (mini + deltah)))[0]
#                ind_term = ind[np.round(len(ind) / 2.).astype(int)]
#            except IndexError:
#                # Last resort
#                ind_term = np.argmin(zoutline)
#    else:
#        # easy: just the minimum
#        ind_term = np.argmin(zoutline)
    
#    return np.asarray(ext_yx)[:, ind_term].astype(np.int64)
    
    