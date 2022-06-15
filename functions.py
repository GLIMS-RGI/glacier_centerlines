"""
Functions to be called in the main centerline script
by froura, Mar 2022
"""
# import libraries --------------------
import shapely.geometry as shpg
import shapely
import numpy as np
from osgeo import gdal
import scipy
from scipy.ndimage.morphology import distance_transform_edt
from scipy.interpolate import RegularGridInterpolator
import copy
from scipy.ndimage.filters import gaussian_filter1d
from functools import partial
import geopandas as gpd
#from main import Centerline

def coordinate_change(tif_path):
    """
    Parameters
    ----------
    tif_path :str
        path to raster file

    Returns
    -------
    Raster values and raster parameters (xOrigin, yOrigin, pixelHeight,
                                         pixelWidth)
    """
    dataset = gdal.Open(tif_path)
    band = dataset.GetRasterBand(1)

    cols = dataset.RasterXSize
    rows = dataset.RasterYSize

    # map pixel/line coordinates into georeferenced space
    transform = dataset.GetGeoTransform()

    xOrigin = transform[0]
    yOrigin = transform[3]
    pixelWidth = transform[1]
    pixelHeight = -transform[5]

    data = band.ReadAsArray(0, 0, cols, rows)
    pix_params = [xOrigin,yOrigin,pixelHeight,pixelWidth]

    return data, pix_params


def profile(points_yx, data, pix_params):
    """
    Parameters
    ----------
    points_list :
        list with lat, lon.
    data :
        np.ndarray, altitude (topography) of each pixel
    pix_params :
        list  with (xorigin, yorigin, pixelH, pixelW)

    Returns
    -------
    tuple: profile distance (arbitrary units) - altitude (m)
    """

    # initialize vectors
    alt = np.zeros(1)
    dist = np.zeros(1)
    dumdist = 0

    xOrigin, yOrigin, pixelHeight, pixelWidth = pix_params

    # altitude
    for point in points_yx:
        col = int((point[0] - xOrigin) / pixelWidth)
        row = int((yOrigin - point[1]) / pixelHeight)

        alt = np.append(alt, data[row][col])

    # remove dummy 0 in the beginning
    alt = alt[1:]

    # distance along line
    # Distance between  2 points

    # repeat the first point at the end
    # np.append(points_list, points_list[0])

    for i in np.arange(len(points_yx)):
        i = int(i)
        a = shpg.Point(points_yx[i])
        # last point
        if i == len(points_yx)-1:
            d = a.distance(shpg.Point(points_yx[0]))
        else:
            d = a.distance(shpg.Point(points_yx[i+1]))
        dumdist = dumdist + d
        dist = np.append(dist, dumdist)

    # remove the dummy 0 ini point
    dist = dist[1:]

    return dist, alt


def get_terminus_coord(ext_yx, zoutline):
    """This finds the terminus coordinate of the glacier.
    There is a special case for marine terminating glaciers/
    Parameters
    ----------
    ext_yx : list
        list with the coordinates (y,x) from the points on the glacier outline
    zoutline : np.ndarray
        altitude of the outline points
    Returns
    -------
    xy - coordinates (shapely.geometry.point.Point) of the glacier terminus.
    index: integer, index of the terminus in the input list.
    """
    # NOTE: possible problems in tidewater because of not constant distance
    # between points

    perc = 10  # from oggm #cfg.PARAMS['terminus_search_percentile']
    deltah = 20  # problem #50m(?)#cfg.PARAMS['terminus_search_altitude_range']

    # if gdir.is_tidewater and (perc > 0):
    if min(zoutline) == 0 and perc > 0:

        plow = np.percentile(zoutline, perc).astype(np.int64)

        # the minimum altitude in the glacier outline
        mini = np.min(zoutline)

        # indices of where in the outline the altitude is lower than the qth
        # percentile and lower than $delatah meters higher,
        # than the minimum altitude
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
    cost[np.where(ext)] = np.nanmax(cost[np.where(ext)]) * 50
    # this works but makes the costgrid plot ugly + I dont get why are we
    # doing it
    # cost[np.where(ext)] = np.nanmax(cost[np.where(ext)]) * 50

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
        raise ValueError('This glacier geometry is not valid for OGGM.')

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
                    raise ValueError('This glacier geometry is not \
                                     valid for OGGM.')
    # TODO incorporate InvalidGeometryError

    if not tmp.is_valid:
        raise ValueError('This glacier geometry is not valid for OGGM.')
    # TODO incorporate InvalidGeometryError

    return tmp


def _filter_lines(lines, heads, k, r):
    """Filter the centerline candidates by length.
    Kienholz et al. (2014), Ch. 4.3.1
    Parameters
    ----------
    lines : list of shapely.geometry.LineString instances
        The lines to filter out (in raster coordinates).
    heads :  list of shapely.geometry.Point instances
        The heads corresponding to the lines. (also in raster coordinates)
    k : float
        A buffer (in raster coordinates) to cut around the selected lines
    r : float
        The lines shorter than r will be filtered out.
    Returns
    -------
    (lines, heads) a list of the new lines and corresponding heads
    """

    olines = []
    oheads = []
    ilines = copy.copy(lines)

    lastline = None
    while len(ilines) > 0:  # loop as long as we haven't filtered all lines
        if len(olines) > 0:  # enter this after the first step only

            toremove = lastline.buffer(k)  # buffer centerlines the last line
            tokeep = []
            for l in ilines:
                # loop over all remaining lines and compute their diff
                # to the last longest line
                diff = l.difference(toremove)
                if diff.type == 'MultiLineString':
                    # Remove the lines that have no head
                    diff = list(diff.geoms)
                    for il in diff:
                        hashead = False
                        for h in heads:
                            if il.intersects(h):
                                hashead = True
                                diff = il
                                break
                        if hashead:
                            break
                        else:
                            diff = None
                # keep this head line only if it's long enough
                if diff is not None and diff.length > r:
                    # Fun fact. The heads can be cut by the buffer too
                    diff = shpg.LineString(l.coords[0:2] + diff.coords[2:])
                    tokeep.append(diff)
            ilines = tokeep

        # it could happen that we're done at this point
        if len(ilines) == 0:
            break

        # Otherwise keep the longest one and continue
        lengths = np.array([])
        for l in ilines:
            lengths = np.append(lengths, l.length)
        ll = ilines[np.argmax(lengths)]

        ilines.remove(ll)
        if len(olines) > 0:
            # the cut line's last point is not guaranteed
            # to on straight coordinates. Remove it
            olines.append(shpg.LineString(np.asarray(ll.xy)[:, 0:-1].T))
        else:
            olines.append(ll)
        lastline = ll

    # add the corresponding head to each line
    for l in olines:
        for h in heads:
            if l.intersects(h):
                oheads.append(h)
                break

    return olines, oheads


def _filter_lines_slope(lines, heads, topo, gdir, min_slope):
    """Filter the centerline candidates by slope: if they go up, remove
    Kienholz et al. (2014), Ch. 4.3.1
    Parameters
    ----------
    lines : list of shapely.geometry.LineString instances
        The lines to filter out (in raster coordinates).
    topo : the glacier topography
    gdir : the glacier directory for simplicity
    min_slope: rad
    Returns
    -------
    (lines, heads) a list of the new lines and corresponding heads
    """
    import params
    dx_cls = params.flowline_dx
    lid = params.flowline_junction_pix
    sw = params.flowline_height_smooth

    # Bilinear interpolation
    # Geometries coordinates are in "pixel centered" convention, i.e
    # (0, 0) is also located in the center of the pixel
    xy = (np.arange(0, gdir.grid.ny-0.1, 1),
          np.arange(0, gdir.grid.nx-0.1, 1))
    interpolator = RegularGridInterpolator(xy, topo)

    olines = [lines[0]]
    oheads = [heads[0]]
    for line, head in zip(lines[1:], heads[1:]):

        # The code below mimics what initialize_flowlines will do
        # this is a bit smelly but necessary
        points = line_interpol(line, dx_cls)

        # For tributaries, remove the tail
        #points = points[0:-lid]

        new_line = shpg.LineString(points)

        # Interpolate heights
        x, y = new_line.xy
        hgts = interpolator((y, x))

        # If smoothing, this is the moment
        hgts = gaussian_filter1d(hgts, sw)

        # Finally slope
        slope = np.arctan(-np.gradient(hgts, dx_cls*gdir.grid.dx))

        # And altitude range
        z_range = np.max(hgts) - np.min(hgts)

        # arbitrary threshold with which we filter the lines, otherwise bye bye
        if np.sum(slope >= min_slope) >= 5 and z_range > 10:
            olines.append(line)
            oheads.append(head)

    return olines, oheads


def _normalize(n):
    """Computes the normals of a vector n.
    Returns
    -------
    the two normals (n1, n2)
    """
    nn = n / np.sqrt(np.sum(n*n))
    n1 = np.array([-nn[1], nn[0]])
    n2 = np.array([nn[1], -nn[0]])
    return n1, n2


def _projection_point(centerline, point):
    """Projects a point on a line and returns the closest integer point
    guaranteed to be on the line, and guaranteed to be far enough from the
    head and tail.
    Parameters
    ----------
    centerline : Centerline instance
    point : Shapely Point geometry
    Returns
    -------
    (flow_point, ind_closest): Shapely Point and indice in the line
    """
    prdis = centerline.line.project(point, normalized=False)
    ind_closest = np.argmin(np.abs(centerline.dis_on_line - prdis)).item()
    flow_point = shpg.Point(centerline.line.coords[int(ind_closest)])
    return flow_point


def line_order(line):
    """Recursive search for the line's hydrological level.
    Parameters
    ----------
    line: a Centerline instance
    Returns
    -------
    The line's order
    """

    if len(line.inflows) == 0:
        return 0
    else:
        levels = [line_order(s) for s in line.inflows]
        return np.max(levels) + 1


def line_interpol(line, dx):
    """Interpolates a shapely LineString to a regularly spaced one.
    Shapely's interpolate function does not guaranty equally
    spaced points in space. This is what this function is for.
    We construct new points on the line but at constant distance from the
    preceding one.
    Parameters
    ----------
    line: a shapely.geometry.LineString instance
    dx: the spacing
    Returns
    -------
    a list of equally distanced points
    """

    # First point is easy
    points = [line.interpolate(dx / 2.)]

    # Continue as long as line is not finished
    while True:
        pref = points[-1]
        pbs = pref.buffer(dx).boundary.intersection(line)
        if pbs.type == 'Point':
            pbs = [pbs]
        elif pbs.type == 'LineString':
            # This is rare
            pbs = [shpg.Point(c) for c in pbs.coords]
            assert len(pbs) == 2
        elif pbs.type == 'GeometryCollection':
            # This is rare
            opbs = []
            for p in pbs.geoms:
                if p.type == 'Point':
                    opbs.append(p)
                elif p.type == 'LineString':
                    opbs.extend([shpg.Point(c) for c in p.coords])
            pbs = opbs
        else:
            if pbs.type != 'MultiPoint':
                raise RuntimeError('line_interpol: we expect a MultiPoint '
                                   'but got a {}.'.format(pbs.type))

        try:
            # Shapely v2 compat
            pbs = pbs.geoms
        except AttributeError:
            pass

        # Out of the point(s) that we get, take the one farthest from the top
        refdis = line.project(pref)
        tdis = np.array([line.project(pb) for pb in pbs])
        p = np.where(tdis > refdis)[0]
        if len(p) == 0:
            break
        points.append(pbs[int(p[0])])

    return points


def gaussian_blur(in_array, size):
    """Applies a Gaussian filter to a 2d array.
    Parameters
    ----------
    in_array : numpy.array
        The array to smooth.
    size : int
        The half size of the smoothing window.
    Returns
    -------
    a smoothed numpy.array
    """

    # expand in_array to fit edge of kernel
    padded_array = np.pad(in_array, size, 'symmetric')

    # build kernel
    x, y = np.mgrid[-size:size + 1, -size:size + 1]
    g = np.exp(-(x**2 / float(size) + y**2 / float(size)))
    g = (g / g.sum()).astype(np.float)  # I had to change that

    # do the Gaussian blur
    return scipy.signal.fftconvolve(padded_array, g, mode='valid')


def _chaikins_corner_cutting(line, refinements=5):
    """Some magic here.
    https://stackoverflow.com/questions/47068504/where-to-find-python-
    implementation-of-chaikins-corner-cutting-algorithm
    """
    coords = np.array(line.coords)

    for _ in range(refinements):
        L = coords.repeat(2, axis=0)
        R = np.empty_like(L)
        R[0] = L[0]
        R[2::2] = L[1:-1:2]
        R[1:-1:2] = L[2::2]
        R[-1] = L[-1]
        coords = L * 0.75 + R * 0.25

    return shpg.LineString(coords)


def save_lines(lines, file_name, crs):
    """
    Given lines and a name of file, centerlines are saved as multiline.
    
    Parameters
    ----------
    lines : list of shapely.geometry.linestring.LineString
        Centerlines coming from the centerline algorithm.
    file_name : string
        Name of file to be saved to.
    crs : pyproj.crs.crs.CRS
        projection to be saved to.
        
    Returns
    -------
    Saves file as 'file_name'. 

    """
    # save as shapefile
    multilinec = shpg.MultiLineString(lines)  # lines merged in a multiline

    # transform to geodataframe and save as shapefile
    gpdmultil = gpd.GeoSeries(multilinec)
    
    # add crs
    gpdmultil.crs = crs
    
    # save object
    gpdmultil.to_file(driver='ESRI Shapefile',
                      filename=file_name)
