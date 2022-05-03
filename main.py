"""
Compute glacier centerlines (flow lines) as in OGGM:
## https://github.com/OGGM/oggm/blob/master/oggm/core/centerlines.py
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
#------------------------ Import functions ---------
from functions import (get_terminus_coord, profile, coordinate_change, 
                       _make_costgrid, _polygon_to_pix, _filter_lines,
                       _filter_lines_slope, _normalize, 
                       _projection_point, line_order, gaussian_blur)

# libraries to save output
import gzip
import pickle

## Comment: from  Kienholtz paper, depression filling is not currently done.
 
# TODO:  
# - error at i=12, see what's happening

# load parameters to be used (more info in params.py)
import params

q1 = params.q1
q2 = params.q2 #m
rmax = params.rmax #m
localmax_window = params.localmax_window #In units of dx
flowline_dx = params.flowline_dx
flowline_junction_pix = params.flowline_junction_pix
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

#

# this is a variable used in the centerline class, but I don't know what 
# is it used for
GAUSSIAN_KERNEL = dict()
for ks in [5, 7, 9]:
    kernel = gaussian(ks, 1)
    GAUSSIAN_KERNEL[ks] = kernel / kernel.sum()

### Define classes to reuse some functions from original OGGM code
from utils import SuperclassMeta

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

### FUNCTION: It may have to go to functions.py but it gave some cross 
# import problems so at the moment I keep it here 
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

# view all glaciers at once:
if plot:
    crop_extent.plot()

# Get altitude and pixel info: 
# altitude, (xorigin, yorigin, pixelH, pixelW)
data, pix_params = coordinate_change(dem_path)

# skip glacier i=12 because here is some problem determining the minimum cost trajectory
#dum = np.arange(12)
#dum = np.append(dum, np.arange(13, len(crop_extent))) 

# loop over all glaciers
for i in np.arange(len(crop_extent)): 
#for i in dum: 
    print(i)
    # select i-th glacier
    crp1 = crop_extent.iloc[[i]]
    
    # crop the DEM to the outline + a few buffer points. Buffer in meters
    dem_clipped = dem.rio.clip(crp1.buffer(20).apply(shpg.mapping),
                               crop_extent.crs)

    # assign some value to outside crop: e.g. 0 (default number is too large)
    #dem_clipped.values[0][dem_clipped.values[0] > 32000] = 0
    dem_clipped.values[0][dem_clipped.values[0] < -32000] = 0
    
    ### Determine heads and tails ###
    area = crop_extent.geometry[i].area
    points_yx=crop_extent.geometry.exterior[i].coords #list of outline X,Y coordinates
    
    # Circular outline: if the first element is repeated at the end, delete the last one
    if points_yx[0] == points_yx[-1]:
        points_yx = points_yx[:-1]
    
    # get profile under glacier outline: distance (arbitrary units) - altitude (m)  
    prof = profile(points_yx, data, pix_params)
    
    # get terminus coordinates and position (index) in our data 
    xyterm, ind_term = get_terminus_coord(points_yx, prof[1]) # get_terminus_coord(yx, zoutline)

    zoutline = prof[1]

    ## here heads start
    # create grid
    nx = len(dem_clipped.x)
    ny = len(dem_clipped.y)
    dx = abs(max(dem_clipped.x) - min(dem_clipped.x))/nx
    dy = abs(max(dem_clipped.y) - min(dem_clipped.y))/ny
    ulx = min(dem_clipped.x)
    uly = max(dem_clipped.y)
    
    ### I dont undersatnd this dx instead of dy ????
    x0y0 = (ulx + dx/2, uly - dx/2) # To pixel center coordinates
    
    # try the curent crs
    utm_proj = salem.check_crs(crop_extent.crs)
        
    ### I dont undersatnd this -dx instead of - dy ????
    grid = salem.Grid(proj=utm_proj, nxny=(nx, ny), dxdy=(dx, -dx), x0y0=x0y0) 
     
    # fill gdir class with grid data
    gdir = glacier_dir(grid)

    # Size of the half window to use to look for local maxima
    maxorder = np.rint(localmax_window/dx) #np.rint(cfg.PARAMS['localmax_window'] / gdir.grid.dx)
    
    #order: number of points to take at each side. minumum 5, maximum = maxorder
    maxorder = utils.clip_scalar(maxorder, 5., np.rint((len(zoutline) / 5.)))
    heads_idx = scipy.signal.argrelmax(zoutline, mode='wrap',
                                       order=int(maxorder))
    if single_fl or len(heads_idx[0]) <= 1:
        # small glaciers with one or less heads: take the absolute max
        heads_idx = (np.atleast_1d(np.argmax(zoutline)),)
    
    # Remove the heads that are too low
    zglacier = dem_clipped.values[dem_clipped.values != 0] 
    head_threshold = np.percentile(zglacier, (1./3.)*100)
    _heads_idx = heads_idx[0][np.where(zoutline[heads_idx] > head_threshold)]
    
    if len(_heads_idx) == 0:
        # this is for bad ice caps where the outline is far off in altitude
        _heads_idx = [heads_idx[0][np.argmax(zoutline[heads_idx])]]
    heads_idx = _heads_idx
    
    headsx = np.array(0)
    headsy = np.array(0)

    for j in heads_idx:
        headsy = np.append(headsy, points_yx[int(j)][1])
        headsx = np.append(headsx, points_yx[int(j)][0])
 
    headsx = headsx[1:]
    headsy = headsy[1:]
    heads_z = zoutline[heads_idx]
    
    # careful, the coords are in y, x order!
    heads = [shpg.Point(x, y) for y, x in zip(headsy,
                                              headsx)]
    
    #radius for descarting heads
    radius = q1 * area + q2 
    radius = utils.clip_scalar(radius, 0, rmax) 
    radius /= grid.dx #in raster coordinates

    ## params taken from default OGGM values  
    # Plus our criteria, quite useful to remove short lines:
    radius += flowline_junction_pix * flowline_dx #cfg.PARAMS['flowline_junction_pix'] * cfg.PARAMS['flowline_dx']
    
    # OK. Filter and see.
    poly_pix = crop_extent.geometry[i]
    
    # heads' coordinates and altitude
    heads, heads_z = _filter_heads(heads, list(heads_z), float(radius), poly_pix)
    
    #if some head removed:
    headsx=np.zeros(1)
    headsy=np.zeros(1)
    for k in np.arange(len(heads)):
        headsy = np.append(headsy, heads[k].y)
        headsx = np.append(headsx, heads[k].x)
         
    headsx = headsx[1:]
    headsy = headsy[1:]
    
    # topography on glacier
    z = dem_clipped.values[0]
    
    # glacier (polygon) (for some reason i have to keep the index 'i' and i cannot use just [0])
    glacier_poly_hr = crp1.geometry[i]

# Rounded nearest pix
    glacier_poly_pix = _polygon_to_pix(glacier_poly_hr)

    tuple2int = partial(np.array, dtype=np.int64)

    # Compute the glacier mask (currently: center pixels + touched)
    #glacier_mask = mask.values[0]
    glacier_mask = np.zeros((ny, nx), dtype=np.uint8)
    glacier_ext = np.zeros((ny, nx), dtype=np.uint8)
    (x, y) = glacier_poly_pix.exterior.xy
    
    #transform coordinates to pixels and assign to 1 inside, 0 otherwise
    xx, yy = grid.transform(x,y,crs=utm_proj)
    
    # I have added a np.clip because some errors when regridding data
    xx = np.clip(xx, 0, glacier_mask.shape[1]-1)
    yy = np.clip(yy, 0, glacier_mask.shape[0]-1)
    glacier_mask[skdraw.polygon(np.array(yy), np.array(xx))] = 1
    
    for gint in glacier_poly_pix.interiors:
         x, y = tuple2int(gint.xy)
         xx, yy = grid.transform(x,y,crs=utm_proj)
         xx, yy = np.round(xx), np.round(yy)
         xx, yy = xx.astype(int), yy.astype(int)
         glacier_mask[skdraw.polygon(yy, xx)] = 0
         glacier_mask[yy, xx] = 0  # on the nunataks
    
    x, y = tuple2int(glacier_poly_pix.exterior.xy)

    #project xy to our shapefile (raster) grid
    xx, yy = grid.transform(x,y,crs=utm_proj)
    xx, yy = np.round(xx), np.round(yy)
    xx, yy = xx.astype(int), yy.astype(int)
    
    # I have added a np.clip because some errors when regridding data
    xx = np.clip(xx, 0, glacier_mask.shape[1]-1)
    yy = np.clip(yy, 0, glacier_mask.shape[0]-1)
    
    glacier_mask[yy, xx] = 1
    glacier_ext[yy, xx] = 1  
    
#####------------------------- Compute centerlines --------------------
    # Cost array
    costgrid = _make_costgrid(glacier_mask, glacier_ext, z)

    # Terminus
    t_coord = xyterm
    
    # Compute the least cost routes
    lines = []
    heads_pix = []
    count = 0
    for h in heads:
        try:
            h_coord_pix = np.array(grid.transform(h.x, h.y, crs=utm_proj))
            h_coord_pix = np.round(h_coord_pix).astype(int)
            t_coord_pix = np.array(grid.transform(t_coord.x, t_coord.y, crs=utm_proj))
            t_coord_pix = np.round(t_coord_pix).astype(int)
            heads_pix.append(shpg.Point(h_coord_pix))
         
            indices, _ = route_through_array(costgrid, np.roll(h_coord_pix,1), np.roll(t_coord_pix,1))
            lines.append(shpg.LineString(np.array(indices)[:, [1, 0]]))
        except:
            print("There is a problem when computing the least-cost route")
            #raise Exception("There is a problem when computing the least-cost route")
            count += 1
        finally:
            print(".")
    print(str(count) + " centerlines out of " + str(len(heads)) + " were not able to be computed" )
    
    # Filter the shortest lines out
    dx_cls = flowline_dx #cfg.PARAMS['flowline_dx']
    radius = flowline_junction_pix * dx_cls #cfg.PARAMS['flowline_junction_pix'] * dx_cls
    radius += 6 * dx_cls
    
    #heads in raster coordinates:   
    olines, oheads = _filter_lines(lines, heads_pix, kbuffer, radius) #cfg.PARAMS['kbuffer'], radius)

    # Filter the lines which are going up instead of down
    min_slope = np.deg2rad(min_slope_flowline_filter)
    do_filter_slope = filter_min_slope

    if do_filter_slope:
        topo = z
        olines, oheads = _filter_lines_slope(olines, oheads, topo,
                                             gdir, min_slope)

    # And rejoin the cut tails
    olines = _join_lines(olines, oheads)

    # Adds the line level
    for cl in olines:
        cl.order = line_order(cl)

    # And sort them per order !!! several downstream tasks  rely on this
    cls = []
    for k in np.argsort([cl.order for cl in olines]):
        cls.append(olines[k])

    # Write centerlines
    # TODO: save it in geographic and not raster coordinates
    #    
    #
    
    use_comp = True
    _open = gzip.open if use_comp else open
    fp =  out_path + "centerline_glacier_" + str(i) + ".pkl"
    with _open(fp, 'wb') as f:
        pickle.dump(cls, f, protocol = 4)

    if plot: #(this is provisional, for checking urposes only)
        # plot profile
        #plot profile + terminus + heads:
        ## NOTE: in this plot, removed heads are displayed anyway.
        #plt.plot(prof[0], zoutline, '-+') #horizontal distance vs altitude
        #plt.plot(prof[0][ind_term], zoutline[ind_term], 'r*', label="terminus") #terminus
        #plt.plot(prof[0][heads_idx], zoutline[heads_idx], 'g*', label="head") #head
        #plt.xlabel("Distance along outline (a.u.)")
        #plt.ylabel("Altitude (m)")
        #plt.legend()
        #plt.show()
            
        #plot altitude raster
        f, ax = plt.subplots(figsize=(8, 10))
#        dem_clipped.plot(ax=ax)
#        ax.set(title="Raster Layer Cropped to Geodataframe Extent")
#        plt.scatter(headsx,headsy, marker="*",s=1000, c="g")
#        plt.scatter(xyterm.x,xyterm.y, marker="*", s=1000, c="r")  
#        crp1.boundary.plot(ax=ax)
#        plt.show()
        
       
        # plot costgrid
        plt.imshow(costgrid)
        plt.colorbar()

        plt.scatter(heads_pix[0].xy[0], heads_pix[0].xy[1], marker="*",s=100, c="g")
        if len(lines) > 1 :
            plt.scatter(heads_pix[1].xy[0], heads_pix[1].xy[1], marker="*",s=100, c="g")
        plt.scatter(t_coord_pix[0], t_coord_pix[1], marker="*", s=100, c="r") 
        
        plt.scatter(lines[0].xy[0], lines[0].xy[1], marker="o", s=5000/(nx*ny), c="y")
        #if len(lines) > 1 :
        #    plt.scatter(lines[1].xy[0], lines.xy[1], marker="o", s=5000/(nx*ny), c="y") 


        for lin in np.arange(len(lines)):
            plt.scatter(lines[lin].xy[0], lines[lin].xy[1], marker="o", s=5000/(nx*ny), c="y") 
        plt.show() 
        


        