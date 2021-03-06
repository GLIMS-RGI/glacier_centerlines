"""
Francesc Roura Adserias @franra9 from OGGM/core/centerlines.py
## https://github.com/OGGM/oggm/blob/master/oggm/core/centerlines.py
Compute glacier centerlines (flow lines) as in OGGM:
and https://tc.copernicus.org/articles/8/503/2014/
https://github.com/OGGM/oggm/blob/447a49d7f936dae4870453d7c65bf2c6f861d0d8/oggm/core/gis.py#L798
"""
import numpy as np
import scipy
import utils as utils
from utils import lazy_property
from utils import _filter_heads
import salem
import shapely.geometry as shpg
import matplotlib.pyplot as plt
import os
import rioxarray as rio
import geopandas as gpd
from shapely.ops import transform as shp_trafo
import shapely.affinity as shpa
import shapely.ops as shpo

try:
    import skimage.draw as skdraw
except ImportError:
    pass
from functools import partial
try:
    from skimage.graph import route_through_array
except ImportError:
    pass
from itertools import groupby
try:
    from scipy.signal.windows import gaussian
except AttributeError:
    # Old scipy
    from scipy.signal import gaussian
# ------------------------ Import functions ---------
from functions import (get_terminus_coord, profile, coordinate_change,
                       _make_costgrid, _polygon_to_pix, _filter_lines,
                       _filter_lines_slope, _normalize,
                       _projection_point, line_order, gaussian_blur,
                       _chaikins_corner_cutting, save_lines)

# TODO: ? from  Kienholtz paper, depression filling is not currently done.

# load parameters to be used (more info in params.py)
import params

# load class
from utils import SuperclassMeta

############################################################
# TODO: ? parse params as in OGGM cfg.PARAMS('parameter-name')

q1 = params.q1
q2 = params.q2  # m
rmax = params.rmax  # m
localmax_window = params.localmax_window  # In units of dx
flowline_dx = params.flowline_dx
flowline_junction_pix = params.flowline_junction_pix = 1
kbuffer = params.kbuffer
min_slope_flowline_filter = params.min_slope_flowline_filter
filter_min_slope = params.filter_min_slope
flowline_height_smooth = params.flowline_height_smooth
is_first_call = params.is_first_call
plot = params.plot
data_path = params.data_path
out_path = params.out_path
dem_file = params.dem_file
shape_file = params.shape_file
smooth_window = params.smooth_window
single_fl = params.single_fl

############################################################

# this is a variable used in the centerline class, but I don't know what
# is it used for
GAUSSIAN_KERNEL = dict()
for ks in [5, 7, 9]:
    kernel = gaussian(ks, 1)
    GAUSSIAN_KERNEL[ks] = kernel / kernel.sum()


##################################################################
# Define classes to reuse some functions from original OGGM code
class Centerline(object, metaclass=SuperclassMeta):
    """Geometry (line and widths) and flow rooting properties, but no thickness
    """

    def __init__(self, line, dx=None, surface_h=None, orig_head=None,
                 rgi_id=None, map_dx=None):
        """ Initialize a Centerline
        Parameters
        ----------
        line : :py:class:`shapely.geometry.LineString`
            The geometrically calculated centerline
        dx : float
            Grid spacing of the initialised flowline in pixel coordinates
        surface_h :  :py:class:`numpy.ndarray`
            elevation [m] of the points on ``line``
        orig_head : :py:class:`shapely.geometry.Point`
            geometric point of the lines head
        rgi_id : str
            The glacier's RGI identifier
        map_dx : float
            the map's grid resolution. Centerline.dx_meter = dx * map_dx
        """

        self.line = None  # Shapely LineString
        self.head = None  # Shapely Point
        self.tail = None  # Shapely Point
        self.dis_on_line = None
        self.nx = None
        if line is not None:
            self.set_line(line)  # Init all previous properties
        else:
            self.nx = len(surface_h)
            self.dis_on_line = np.arange(self.nx) * dx

        self.order = None  # Hydrological flow level (~ Strahler number)

        # These are computed at run time by compute_centerlines
        self.flows_to = None  # pointer to a Centerline object (when available)
        self.flows_to_point = None  # point of the junction in flows_to
        self.inflows = []  # list of Centerline instances (when available)
        self.inflow_points = []  # junction points

        # Optional attrs
        self.dx = dx  # dx in pixels (assumes the line is on constant dx
        self.map_dx = map_dx  # the pixel spacing
        try:
            self.dx_meter = self.dx * self.map_dx
        except TypeError:
            # For backwards compatibility we allow this for now
            self.dx_meter = None
        self._surface_h = surface_h
        self._widths = None
        self.is_rectangular = None
        self.is_trapezoid = None
        self.orig_head = orig_head  # Useful for debugging and for filtering
        self.geometrical_widths = None  # these are kept for plotting and such
        self.apparent_mb = None  # Apparent MB, NOT weighted by width.
        self.mu_star = None  # the mu* associated with the apparent mb
        self.mu_star_is_valid = False  # if mu* leeds to good flux, keep it
        self.flux = None  # Flux (kg m-2)
        self.flux_needs_correction = False  # whether this branch was baaad
        self.rgi_id = rgi_id  # Useful if line is used with another glacier

    def set_flows_to(self, other, check_tail=True, to_head=False):
        """Find the closest point in "other" and sets all the corresponding
        attributes. Btw, it modifies the state of "other" too.
        Parameters
        ----------
        other : :py:class:`oggm.Centerline`
            another flowline where self should flow to
        """

        self.flows_to = other

        if check_tail:
            # Project the point and Check that its not too close
            prdis = other.line.project(self.tail, normalized=False)
            ind_closest = np.argmin(np.abs(other.dis_on_line - prdis)).item()
            n = len(other.dis_on_line)
            if n >= 9:
                ind_closest = utils.clip_scalar(ind_closest, 4, n-5)
            elif n >= 7:
                ind_closest = utils.clip_scalar(ind_closest, 3, n-4)
            elif n >= 5:
                ind_closest = utils.clip_scalar(ind_closest, 2, n-3)
            p = shpg.Point(other.line.coords[int(ind_closest)])
            self.flows_to_point = p
        elif to_head:
            self.flows_to_point = other.head
        else:
            # just the closest
            self.flows_to_point = _projection_point(other, self.tail)
        other.inflow_points.append(self.flows_to_point)
        other.inflows.append(self)

    def set_line(self, line):
        """Update the Shapely LineString coordinate.
        Parameters
        ----------
        line : :py:class`shapely.geometry.LineString`
        """

        self.nx = len(line.coords)
        self.line = line
        dis = [line.project(shpg.Point(co)) for co in line.coords]
        self.dis_on_line = np.array(dis)
        xx, yy = line.xy
        self.head = shpg.Point(xx[0], yy[0])
        self.tail = shpg.Point(xx[-1], yy[-1])

    @lazy_property
    def flows_to_indice(self):
        """Indices instead of geometry"""

        ind = []
        tofind = self.flows_to_point.coords[0]
        for i, p in enumerate(self.flows_to.line.coords):
            if p == tofind:
                ind.append(i)
        assert len(ind) == 1, 'We expect exactly one point to be found here.'
        return ind[0]

    @lazy_property
    def inflow_indices(self):
        """Indices instead of geometries"""

        inds = []
        for p in self.inflow_points:
            ind = [i for (i, pi) in enumerate(self.line.coords)
                   if (p.coords[0] == pi)]
            inds.append(ind[0])
        assert len(inds) == len(self.inflow_points), ('For every inflow point '
                                                      'there should be '
                                                      'exactly one inflow '
                                                      'indice')
        return inds

    @lazy_property
    def normals(self):
        """List of (n1, n2) normal vectors at each point.
        We use second order derivatives for smoother widths.
        """

        pcoords = np.array(self.line.coords)

        normals = []
        # First
        normal = np.array(pcoords[1, :] - pcoords[0, :])
        normals.append(_normalize(normal))
        # Second
        normal = np.array(pcoords[2, :] - pcoords[0, :])
        normals.append(_normalize(normal))
        # Others
        for (bbef, bef, cur, aft, aaft) in zip(pcoords[:-4, :],
                                               pcoords[1:-3, :],
                                               pcoords[2:-2, :],
                                               pcoords[3:-1, :],
                                               pcoords[4:, :]):
            normal = np.array(aaft + 2*aft - 2*bef - bbef)
            normals.append(_normalize(normal))
        # One before last
        normal = np.array(pcoords[-1, :] - pcoords[-3, :])
        normals.append(_normalize(normal))
        # Last
        normal = np.array(pcoords[-1, :] - pcoords[-2, :])
        normals.append(_normalize(normal))

        return normals

    @property
    def widths(self):
        """Needed for overriding later"""
        return self._widths

    @property
    def widths_m(self):
        return self.widths * self.map_dx

    @widths.setter
    def widths(self, value):
        self._widths = value

    @property
    def surface_h(self):
        """Needed for overriding later"""
        return self._surface_h

    @surface_h.setter
    def surface_h(self, value):
        self._surface_h = value

    def set_apparent_mb(self, mb, mu_star=None):
        """Set the apparent mb and flux for the flowline.
        MB is expected in kg m-2 yr-1 (= mm w.e. yr-1)
        This should happen in line order, otherwise it will be wrong.
        Parameters
        ----------
        mu_star : float
            if appropriate, the mu_star associated with this apparent mb
        """

        self.apparent_mb = mb
        self.mu_star = mu_star

        # Add MB to current flux and sum
        # no more changes should happen after that
        flux_needs_correction = False
        flux = np.cumsum(self.flux + mb * self.widths * self.dx)

        # We filter only lines with two negative grid points,
        # the rest we can cope with
        if flux[-2] < 0:
            flux_needs_correction = True

        self.flux = flux
        self.flux_needs_correction = flux_needs_correction

        # Add to outflow. That's why it should happen in order
        if self.flows_to is not None:
            n = len(self.flows_to.line.coords)
            ide = self.flows_to_indice
            if n >= 9:
                gk = GAUSSIAN_KERNEL[9]
                self.flows_to.flux[ide-4:ide+5] += gk * flux[-1]
            elif n >= 7:
                gk = GAUSSIAN_KERNEL[7]
                self.flows_to.flux[ide-3:ide+4] += gk * flux[-1]
            elif n >= 5:
                gk = GAUSSIAN_KERNEL[5]
                self.flows_to.flux[ide-2:ide+3] += gk * flux[-1]


# class defined to be able to use some functions from OGGM "as they are".
class glacier_dir(object):
    def __init__(self, grid):
        self.grid = grid


class grid_inf(object):
    def __init__(self, grid):
        self.grid = grid

     
#####################################################################
# FUNCTION: It may have to go to functions.py but it gave some cross
# import problems so at the moment I keep it here
# def _join_lines(lines, heads):
#     """Re-joins the lines that have been cut by _filter_lines
#      Compute the rooting scheme.
#     Parameters
#     ----------
#     lines: list of shapely lines instances
#     Returns
#     -------
#     Centerline instances, updated with flow routing properties
#      """

#     olines = [Centerline(l, orig_head=h) for l, h
#               in zip(lines[::-1], heads[::-1])]
#     nl = len(olines)
#     if nl == 1:
#         return olines

#     # per construction the line cannot flow in a line placed before in the list
#     for i, l in enumerate(olines):

#         last_point = shpg.Point(*l.line.coords[-1])

#         totest = olines[i+1:]
#         dis = [last_point.distance(t.line) for t in totest]
#         flow_to = totest[np.argmin(dis)]

#         flow_point = _projection_point(flow_to, last_point)

#         # Interpolate to finish the line, bute force:
#         # we interpolate 20 points, round them, remove consecutive duplicates
#         endline = shpg.LineString([last_point, flow_point])
#         endline = shpg.LineString([endline.interpolate(x, normalized=True)
#                                    for x in np.linspace(0., 1., num=20)])
#         # we keep all coords without the first AND the last
#         #grouped = groupby(map(tuple, np.rint(endline.coords)))
#         grouped = groupby(map(tuple, endline.coords))
#         endline = [x[0] for x in grouped][1:-1]
#         #endline = [x[0] for x in grouped][:] # take all points

#         # We're done
#         l.set_line(shpg.LineString(l.line.coords[:] + endline))
#         l.set_flows_to(flow_to, check_tail=False)

#         # The last one has nowhere to flow
#         if i+2 == nl:
#             break

#     return olines[::-1]

def cls_intersec_outline(cls, outline):
    """
    fine tunning: flowlines must finish onto an outline

    Parameters
    ----------
    lines : list of centerline instances
        flowlines
    outline : shapely.geometry.polygon.Polygon in raster coordinates
        glacier boundaries

    Returns
    -------
    modified lines, in form of shapely.lines list (?)

    """
    lis = []
    for li in cls:
        bbb = shpo.split(li.line, outline)
        print(len(bbb.geoms))
        print(li.line.length)
        print(len(bbb.geoms))
        if len(bbb.geoms) != 1:
            
            if len(bbb.geoms) == 2: #intersect head only
                ind=np.argmax([bbb.geoms[0].length, bbb.geoms[1].length])
                print(list(bbb.geoms))
                lis.append(bbb.geoms[ind])
            
            elif len(bbb.geoms) == 3: # intersect in the head and tail both
                 ind=np.argmax([bbb.geoms[0].length, 
                                bbb.geoms[1].length,
                                bbb.geoms[2].length])
                 lis.append(bbb.geoms[ind])

            else:
                 print("Some centerline crosses the outline \
                                 more than 2 times")
        else:
            lis.append(li.line)
            
    return lis

# summ=0
# for l in cls:
#     summ =+ l.line.length
    
# print(summ)
 
#type(aaa.geoms[0])
# for i, feat_line in enumerate(feats_lines):
#     for j, feat_polygon in enumerate(feats_polygons):
#         if feat_polygon.geometry().intersects(feat_line.geometry()):
#             geom = feat_polygon.geometry().intersection(feat_line.geometry())
#             print i, j, (geom.length()/feat_line.geometry().length())*100
#             feat = QgsFeature()
#             feat.setAttributes([i])
#             feat.setGeometry(geom)
#             prov.addFeatures( [feat] )

    
##################
# move this to utils:
    
def geoline_to_cls(lines):
    """
    list of shapely.geometry.lines to be converted to Centerline list 

    Parameters
    ----------
    lines : list of shapely.geometry.lines
        list of shapely.geometry.lines

    Returns
    -------
    list of centerlines instances.

    """
    clss = []
    for li in lines:
        clss.append(Centerline(li))
    cls = clss
    return cls


################################################################
# read the geotiff (DEM) with rioxarray.
dem_path = os.path.join(data_path, dem_file)
dem = rio.open_rasterio(dem_path)

# smooth DEM using gaussian smoothing.
dx = abs(dem.x[0] - dem.x[1])
gsize = int(np.rint(smooth_window / dx))
smoothed_dem = gaussian_blur(np.array(dem.values[0]), int(gsize))
dem.values[0] = smoothed_dem

# read shapefile (glacier outlines)
crop_extent = gpd.read_file(os.path.join(data_path, shape_file))

# Check that projection is in metre
import warnings
try:
    assert crop_extent.crs.axis_info[0].unit_name == 'metre'
    # TODO: dem crs check in meters
    #assert crop_extent.crs.axis_info[0].unit_name == 'metre'        
except Exception:
    warnings.warn('Projection from input data is not in meters.')

# view all glaciers at once:
if plot:
    crop_extent.plot()

# Get altitude and pixel info: 
# altitude, (xorigin, yorigin, pixelH, pixelW)
data, pix_params = coordinate_change(dem_path)

# loop over all glaciers
for i in np.arange(len(crop_extent)):
    print(i)
    # select i-th glacier
    crp1 = crop_extent.iloc[[i]]

    # crop the DEM to the outline + a few buffer points. Buffer in meters
    dem_clipped = dem.rio.clip(crp1.buffer(20).apply(shpg.mapping),
                               crop_extent.crs)

    # assign some value to outside crop: e.g. -1 (default number is too large)
    dem_clipped.values[0][dem_clipped.values[0] < 0 ] = -1
    dem_clipped.values[0][dem_clipped.values[0] > 8848 ] = -1                           

    # Determine heads and tails #
    area = crop_extent.geometry[i].area

    # list of outline X,Y coordinates
    points_xy = crop_extent.geometry.exterior[i].coords

    # Circular outline: if the first element is repeated at the end,
    # delete the last one
    if points_xy[0] == points_xy[-1]:
        points_xy = points_xy[:-1]

    # get profile under glacier outline: distance
    # (arbitrary units (?)) - altitude (m)
    prof = profile(points_xy, data, pix_params)

    # get terminus coordinates and position (index) in our data
    xyterm, ind_term = get_terminus_coord(points_xy, prof[1])

    zoutline = prof[1]

    # here heads start
    # create grid
    nx, ny = len(dem_clipped.x), len(dem_clipped.y)

    ulx, uly, dx, dy = pix_params

    # Initial topmost-leftmost point
    # To pixel center coordinates
    x0y0 = (ulx + dx/2, uly - dy/2)

    # try the shapefile curent crs for the raster grid
    utm_proj = salem.check_crs(crop_extent.crs)

    # build raster grid properties
    grid = salem.Grid(proj=utm_proj, nxny=(nx, ny), dxdy=(dx, -dy), x0y0=x0y0)

    # fill gdir class with grid data
    gdir = glacier_dir(grid)

    # Size of the half window to use to look for local maxima
    maxorder = np.rint(localmax_window/dx)

    # Order: number of points to take at each side.
    # minumum 5, maximum = maxorder
    maxorder = utils.clip_scalar(maxorder, 5., np.rint((len(zoutline) / 5.)))
    heads_idx = scipy.signal.argrelmax(zoutline, mode='wrap',
                                       order=int(maxorder))
    if single_fl or len(heads_idx[0]) <= 1:
        # small glaciers with one or less heads: take the absolute max
        heads_idx = (np.atleast_1d(np.argmax(zoutline)),)

    # TODO: if is tidewater:
    #    do something

    # Remove the heads that are too low
    zglacier = dem_clipped.values[dem_clipped.values != -1]
    head_threshold = np.percentile(zglacier, (1./3.)*100)
    _heads_idx = heads_idx[0][np.where(zoutline[heads_idx] > head_threshold)]

    if len(_heads_idx) == 0:
        # this is for bad ice caps where the outline is far off in altitude
        _heads_idx = [heads_idx[0][np.argmax(zoutline[heads_idx])]]
    heads_idx = _heads_idx

    # TODO: undo the loop, do something like this
    # headsxx = gpd.GeoDataFrame(points_xy)[0]
    # headsyy = gpd.GeoDataFrame(points_x)[1]

    headsx = headsy = np.array(0)

    for j in heads_idx:
        headsy = np.append(headsy, points_xy[int(j)][1])
        headsx = np.append(headsx, points_xy[int(j)][0])

    headsx, headsy = headsx[1:], headsy[1:]
    heads_z = zoutline[heads_idx]

    # careful, the coords are in y, x order!
    heads = [shpg.Point(x, y) for y, x in zip(headsy,
                                              headsx)]

    # radius for descarting heads
    radius = q1 * area + q2
    radius = utils.clip_scalar(radius, 0, rmax)
    radius /= grid.dx   # radius in raster coordinates
    # (warning: this is assuming |xd| = |dy|)

    # parameters taken from default OGGM values
    # Plus our criteria, quite useful to remove short lines:
    radius += flowline_junction_pix * flowline_dx

    # OK. Filter and see.
    glacier_poly_hr = crop_extent.geometry[i]

    # heads' coordinates and altitude
    heads, heads_z = _filter_heads(heads, list(heads_z), float(radius),
                                   glacier_poly_hr)

    # after head filtering, refill heads xy coordinates:
    headsx = headsy = np.zeros(1)

    for khead in heads:
        headsy = np.append(headsy, khead.y)
        headsx = np.append(headsx, khead.x)

    headsx, headsy = headsx[1:], headsy[1:]

    # topography on glacier
    z = dem_clipped.values[0]

# Rounded nearest pix
    glacier_poly_pix = _polygon_to_pix(glacier_poly_hr)

    tuple2int = partial(np.array, dtype=np.int64)

    # Compute the glacier mask (currently: center pixels + touched)
    glacier_mask = np.zeros((ny, nx), dtype=np.uint8)
    glacier_ext = np.zeros((ny, nx), dtype=np.uint8)
    (x, y) = glacier_poly_pix.exterior.xy

    # transform coordinates to pixels and assign to 1 inside, 0 otherwise
    xx, yy = grid.transform(x, y, crs=utm_proj)

    # I have added a np.clip because some errors when regridding data
    xx = np.clip(xx, 0, glacier_mask.shape[1] - 1)
    yy = np.clip(yy, 0, glacier_mask.shape[0] - 1)
    glacier_mask[skdraw.polygon(np.array(yy), np.array(xx))] = 1 #1 inside plygon

    for gint in glacier_poly_pix.interiors:
        x, y = tuple2int(gint.xy)
        xx, yy = grid.transform(x, y, crs=utm_proj)
        xx, yy = np.round(xx), np.round(yy)
        xx, yy = xx.astype(int), yy.astype(int)
        glacier_mask[skdraw.polygon(yy, xx)] = 0  # inside the nunataks
        glacier_mask[yy, xx] = 0  # onto nunatacks boundaries

    x, y = tuple2int(glacier_poly_pix.exterior.xy)

    # project xy to our local (raster) grid
    xx, yy = grid.transform(x, y, crs=utm_proj)
    xx, yy = np.round(xx), np.round(yy)
    xx, yy = xx.astype(int), yy.astype(int)

    # I have added a np.clip because some errors when regridding data
    xx = np.clip(xx, 0, glacier_mask.shape[1]-1)
    yy = np.clip(yy, 0, glacier_mask.shape[0]-1)

    glacier_mask[yy, xx] = 1  # glacier bundaries belong to the glacier
    glacier_ext[yy, xx] = 1

############################################################
# Compute centerlines
    # Cost array
    costgrid = _make_costgrid(glacier_mask, glacier_ext, z)

    # Terminus
    t_coord = xyterm

    # Compute the least cost routes for all possible head-terminus trajectories
    lines = []
    heads_pix = []
    count = 0
    for h in heads:
        try:
            h_coord_pix = np.array(grid.transform(h.x, h.y, crs=utm_proj))
            h_coord_pix = np.round(h_coord_pix).astype(int)
            t_coord_pix = np.array(grid.transform(t_coord.x,
                                                  t_coord.y, crs=utm_proj))
            t_coord_pix = np.round(t_coord_pix).astype(int)
            heads_pix.append(shpg.Point(h_coord_pix))

            indices, _ = route_through_array(costgrid, np.roll(h_coord_pix, 1),
                                             np.roll(t_coord_pix, 1))
            lines.append(shpg.LineString(np.array(indices)[:, [1, 0]]))
        except:
            print("There is a problem when computing the least-cost route")
            # raise Exception("There is a problem when computing
            # the least-cost route")
            count += 1
        #finally:
        #    print('')
    print(str(count) + " centerlines out of " + str(len(heads)) +
          " were not able to be computed")

###############################
 # TODO: line postprocessing according to https://github.com/OGGM/oggm/blob/e7c90e9cb90c0a8a1767e4111fb2d18cc01ac007/oggm/utils/_workflow.py#L722-L748

    flowline_junction_pix
    # Filter the shortest lines out
    dx_cls = flowline_dx
    radius = flowline_junction_pix * dx_cls
    radius += 6 * dx_cls

###############################

    # Smooth centerlines  
    tra_func = partial(gdir.grid.transform, crs=crop_extent.crs)
    exterior = shpg.Polygon(shp_trafo(tra_func, crop_extent.geometry[0].exterior))
    #exterior = shpg.Polygon(crop_extent.geometry[0].exterior)
    
    liness = []
    for j, line in enumerate(lines):
        mm = 1 if j == (len(lines)-1) else 0
    
        ensure_exterior_match = True
        if ensure_exterior_match:
            # Extend line at the start by 10
            fs = shpg.LineString(line.coords[:2])
            # First check if this is necessary - this segment should
            # be within the geometry or it's already good to go
            if fs.within(exterior):
                fs = shpa.scale(fs, xfact=3, yfact=3, origin=fs.boundary.geoms[1])
                line = shpg.LineString([*fs.coords, *line.coords[2:]])
            # If last also extend at the end
            if mm == 1:  # mm means main
                ls = shpg.LineString(line.coords[-2:])
                if ls.within(exterior):
                    ls = shpa.scale(ls, xfact=3, yfact=3, origin=ls.boundary.geoms[0])
                    line = shpg.LineString([*line.coords[:-2], *ls.coords])
        
            # Simplify and smooth?
            simplify_line = True
            if simplify_line:
                line = line.simplify(simplify_line)
            corner_cutting = True
            if corner_cutting:
                line = _chaikins_corner_cutting(line, refinements=5)
        
            # Intersect with exterior geom
            line = line.intersection(exterior)
            if line.type == 'MultiLineString':
                # Take the longest
                lens = [il.length for il in line.geoms]
                line = line.geoms[np.argmax(lens)]            
        liness.append(line)
    
    lines = liness

    # heads in raster coordinates:
    olines, oheads = _filter_lines(lines, heads_pix, kbuffer, radius)
    

    # Filter the lines which are going up instead of down
    min_slope = np.deg2rad(min_slope_flowline_filter)

    if filter_min_slope: # this kills most of the branches! --> too many killed!
        topo = z
        olines, oheads = _filter_lines_slope(olines, oheads, topo,
                                             gdir, min_slope)
# tmp: convert to cls
    olines = geoline_to_cls(olines)
#
    # Adds the line level
    for cl in olines:
        cl.order = line_order(cl)

    # And sort them per order !!! several downstream tasks rely on this
    cls = []
    for k in np.argsort([cl.order for cl in olines]):
        cls.append(olines[k])


###########################################################################

    lines = cls_intersec_outline(cls, exterior)

    cls = geoline_to_cls(lines)
    
    ###### the error about the non matching centerlines-outlines is here!!!
    # these 2 lines match
    #save_lines(utils.cls_to_geoline(cls),'centerlines_tmp.shp', crop_extent.crs)
    #save_lines(exterior.boundary,'boundary.shp', crop_extent.crs)
    
    cls_xy = []
    for li in cls:
        # ij_to_xy
        ii = np.array(li.line.xy[0])
        jj = np.array(li.line.xy[1])
        x = ulx + ii * dx + 0.5 * dx
        y = uly - jj * dy - 0.5 * dy

        xy = np.zeros((len(x), 2))

        xy[:, 0] = x
        xy[:, 1] = y

        lxy = shpg.LineString(xy)

        cls_xy.append(lxy)
   
    save_lines(cls_xy, out_path + "/11_rgi60_central_europe.shp", crop_extent.crs)


###############################################################
    if plot:  # (this is provisional, for checking urposes only)
        # plot profile
        # plot profile + terminus + heads:
        # NOTE: in this plot, removed heads are displayed anyway.
        plt.plot(prof[0], zoutline, '-+')   # horizontal distance vs altitude
        plt.plot(prof[0][ind_term], zoutline[ind_term],
                 'r*', label="terminus")  # terminus
        plt.plot(prof[0][heads_idx], zoutline[heads_idx],
                 'g*', label="head")  # head
        plt.xlabel("Distance along outline (a.u.)")
        plt.ylabel("Altitude (m)")
        plt.legend()
        plt.show()

        # plot altitude raster
        f, ax = plt.subplots(figsize=(8, 10))
        dem_clipped.plot(ax=ax)
        ax.set(title="Raster Layer Cropped to Geodataframe Extent")
        plt.scatter(headsx, headsy, marker="*", s=1000, c="g")
        plt.scatter(xyterm.x, xyterm.y, marker="*", s=1000, c="r")
        crp1.boundary.plot(ax=ax)
        plt.show()

        # plot costgrid
        plt.imshow(costgrid)
        plt.colorbar()

        plt.scatter(heads_pix[0].xy[0], heads_pix[0].xy[1],
                    marker="*", s=100, c="g")
        if len(lines) > 1:
            plt.scatter(heads_pix[1].xy[0], heads_pix[1].xy[1],
                        marker="*", s=100, c="g")
        plt.scatter(t_coord_pix[0], t_coord_pix[1], marker="*",
                    s=100, c="r")

        plt.scatter(lines[0].xy[0], lines[0].xy[1], marker="o",
                    s=5000/(nx*ny), c="y")

        for lin in np.arange(len(lines)):
            plt.scatter(lines[lin].xy[0], lines[lin].xy[1], marker="o",
                        s=5000/(nx*ny), c="y")
        plt.show()
##############################################################
# END #
