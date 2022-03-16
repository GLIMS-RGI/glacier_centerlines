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
import copy
from scipy.interpolate import RegularGridInterpolator
from utils import line_interpol
from scipy.ndimage.filters import gaussian_filter1d

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
    #this works but makes the costgrid plot ugly
    #cost[np.where(ext)] = np.nanmax(cost[np.where(ext)]) * 50
    

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

    dx_cls = flowline_dx = 2 #= cfg.PARAMS['flowline_dx']
    lid = flowline_junction_pix = int(3) #int(cfg.PARAMS['flowline_junction_pix'])
    sw = flowline_height_smooth = 1 #cfg.PARAMS['flowline_height_smooth']

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
        points = points[0:-lid]

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

def _join_lines(lines, heads):
    """Re-joins the lines that have been cut by _filter_lines
     Compute the rooting scheme.
    Parameters
    ----------
    lines: list of shapely lines instances
    Returns
    -------
    Centerline instances, updated with flow routing properties
     """

    olines = [Centerline(l, orig_head=h) for l, h
              in zip(lines[::-1], heads[::-1])]
    nl = len(olines)
    if nl == 1:
        return olines

    # per construction the line cannot flow in a line placed before in the list
    for i, l in enumerate(olines):

        last_point = shpg.Point(*l.line.coords[-1])

        totest = olines[i+1:]
        dis = [last_point.distance(t.line) for t in totest]
        flow_to = totest[np.argmin(dis)]

        flow_point = _projection_point(flow_to, last_point)

        # Interpolate to finish the line, bute force:
        # we interpolate 20 points, round them, remove consecutive duplicates
        endline = shpg.LineString([last_point, flow_point])
        endline = shpg.LineString([endline.interpolate(x, normalized=True)
                                   for x in np.linspace(0., 1., num=20)])
        # we keep all coords without the first AND the last
        grouped = groupby(map(tuple, np.rint(endline.coords)))
        endline = [x[0] for x in grouped][1:-1]

        # We're done
        l.set_line(shpg.LineString(l.line.coords[:] + endline))
        l.set_flows_to(flow_to, check_tail=False)

        # The last one has nowhere to flow
        if i+2 == nl:
            break

    return olines[::-1]