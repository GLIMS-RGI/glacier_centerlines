"""
Functions to be called in the main centerline script
by froura, Mar 2022
"""
# import libraries --------------------
import shapely.geometry as shpg
import shapely
import numpy as np
from osgeo import gdal
from scipy.ndimage.morphology import distance_transform_edt

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
    # NOTE: possible problems in tidewater because of not constant distance between points
    
    perc = 10 # from oggm #cfg.PARAMS['terminus_search_percentile']
    deltah = 20 #20 problem #50m (?) #cfg.PARAMS['terminus_search_altitude_range']

    #if gdir.is_tidewater and (perc > 0):
    if min(zoutline) == 0 and perc > 0:

        plow = np.percentile(zoutline, perc).astype(np.int64)

        # the minimum altitude in the glacier outline
        mini = np.min(zoutline)

        # indices of where in the outline the altitude is lower than the qth
        # percentile and lower than $delatah meters higher, than the minimum altitude
        ind = np.where((zoutline < plow) & (zoutline < (mini + deltah)))[0]

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

def _make_costgrid(mask, ext, z):
    """Computes a costgrid following Kienholz et al. (2014) Eq. (2)
    Parameters
    ----------
    mask : numpy.array
        The glacier mask.
    ext : numpy.array
        The glacier boundaries' mask.
    z : numpy.array
        The terrain height.
    Returns
    -------
    numpy.array of the costgrid
    """
    # Kienholz et al eq (2)
    f1 = 1000.
    f2 = 3000.
    a = 4.25
    b = 3.7  
    
    dis = np.where(mask, distance_transform_edt(mask), np.NaN)
    z = np.where(mask, z, np.NaN)

    dmax = np.nanmax(dis)
    zmax = np.nanmax(z)
    zmin = np.nanmin(z)
    cost = ((dmax - dis) / dmax * f1) ** a + \
           ((z - zmin) / (zmax - zmin) * f2) ** b
#    cost = ((dmax - dis) / dmax * cfg.PARAMS['f1']) ** cfg.PARAMS['a'] + \
#           ((z - zmin) / (zmax - zmin) * cfg.PARAMS['f2']) ** cfg.PARAMS['b']

    # This is new: we make the cost to go over boundaries
    # arbitrary high to avoid the lines to jump over adjacent boundaries
#    cost[np.where(ext)] = np.nanmax(cost[np.where(ext)]) * 50
    cost[0][np.where(ext)] = np.nanmax(cost[0][np.where(ext)]) * 50
    

    return np.where(mask, cost, np.Inf)

def _polygon_to_pix(polygon):
    """Transforms polygon coordinates to integer pixel coordinates. It makes
    the geometry easier to handle and reduces the number of points.
    Parameters
    ----------
    polygon: the shapely.geometry.Polygon instance to transform.
    Returns
    -------
    a shapely.geometry.Polygon class instance.
    """

    def project(x, y):
        return np.rint(x).astype(np.int64), np.rint(y).astype(np.int64)

    def project_coarse(x, y, c=2):
        return ((np.rint(x/c)*c).astype(np.int64),
                (np.rint(y/c)*c).astype(np.int64))

    poly_pix = shapely.ops.transform(project, polygon)

    # simple trick to correct invalid polys:
    tmp = poly_pix.buffer(0)

    # try to deal with a bug in buffer where the corrected poly would be null
    c = 2
    while tmp.length == 0 and c < 7:
        project = partial(project_coarse, c=c)
        poly_pix = shapely.ops.transform(project_coarse, polygon)
        tmp = poly_pix.buffer(0)
        c += 1

    # We tried all we could
    if tmp.length == 0:
        raise InvalidGeometryError('This glacier geometry is not valid for '
                                   'OGGM.')

    # sometimes the glacier gets cut out in parts
    if tmp.type == 'MultiPolygon':
        # If only small arms are cut out, remove them
        area = np.array([_tmp.area for _tmp in tmp.geoms])
        _tokeep = np.argmax(area).item()
        tmp = tmp.geoms[_tokeep]

        # check that the other parts really are small,
        # otherwise replace tmp with something better
        area = area / area[_tokeep]
        for _a in area:
            if _a != 1 and _a > 0.05:
                # these are extremely thin glaciers
                # eg. RGI40-11.01381 RGI40-11.01697 params.d1 = 5. and d2 = 8.
                # make them bigger until its ok
                for b in np.arange(0., 1., 0.01):
                    tmp = shapely.ops.transform(project, polygon.buffer(b))
                    tmp = tmp.buffer(0)
                    if tmp.type == 'MultiPolygon':
                        continue
                    if tmp.is_valid:
                        break
                if b == 0.99:
                    raise InvalidGeometryError('This glacier geometry is not '
                                               'valid for OGGM.')

    if not tmp.is_valid:
        raise InvalidGeometryError('This glacier geometry is not valid '
                                   'for OGGM.')

    return tmp