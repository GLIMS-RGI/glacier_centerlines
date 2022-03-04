"""
Functions to be called in the main centerline script
by froura, Mar 2022
"""
# import libraries --------------------
import shapely.geometry as shpg
import numpy as np
from osgeo import gdal


def coordinate_change(tif_path):
    dataset = gdal.Open(tif_path)
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
    pix_params= [xOrigin,yOrigin,pixelHeight,pixelWidth]
    
    return data, pix_params

def profile(points_yx, data, pix_params):
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
    
    xOrigin, yOrigin, pixelHeight, pixelWidth = pix_params
    
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

