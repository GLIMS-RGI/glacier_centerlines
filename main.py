"""
Francesc Roura Adserias @franra9 from OGGM/core/centerlines.py
## https://github.com/OGGM/oggm/blob/master/oggm/core/centerlines.py

Compute glacier centerlines (flow lines) as in here:
https://tc.copernicus.org/articles/8/503/2014/
"""
import numpy as np
import scipy
import utils as utils
import salem
import shapely.geometry as shpg
import os
import rioxarray as rio
import geopandas as gpd
from shapely.ops import transform as shp_trafo
import shapely.affinity as shpa

try:
    import skimage.draw as skdraw
except ImportError:
    pass
from functools import partial
try:
    from skimage.graph import route_through_array
except ImportError:
    pass
try:
    from scipy.signal.windows import gaussian
except AttributeError:
    # Old scipy
    from scipy.signal import gaussian
# ------------------------ Import functions ---------
from functions import (get_terminus_coord, profile, coordinate_change,
                       _make_costgrid,_filter_lines,_filter_lines_slope,
                       gaussian_blur) # TODO: make costgrid has some issue when importing it from OGGM because I dont have the same params.py file\
    #idem for _filter_lines_slope
                       
from oggm.core.gis import (_polygon_to_pix)#,gaussian_blur) --> problem with oggm blur
from oggm.core.centerlines import (line_order)
from oggm.utils._workflow import _chaikins_corner_cutting

# load parameters to be used (more info in params.py)
import params

############################################################
# TODO: ? parse params as in OGGM cfg.PARAMS('parameter-name')
import oggm.cfg as cfg 

q1 = params.q1 #cfg.PARAMS['q1']
q2 = params.q2  # cfg.PARAMS['q2']
rmax = params.rmax  # cfg.PARAMS['rmax']
localmax_window = params.localmax_window #cfg.PARAMS['localmax_window']  # In units of dx
flowline_dx = params.flowline_dx #cfg.PARAMS['flowline_dx']
flowline_junction_pix = params.flowline_junction_pix = 1 #cfg.PARAMS['flowline_junction_pix']
kbuffer = params.kbuffer #cfg.PARAMS['kbuffer']
min_slope_flowline_filter = params.min_slope_flowline_filter #cfg.PARAMS['slope_flowline_filter']
filter_min_slope = params.filter_min_slope #cfg.PARAMS['filter_min_slope']
flowline_height_smooth = params.flowline_height_smooth #cfg.PARAMS['flowline_height_smooth']
is_first_call = params.is_first_call #cfg.PARAMS['is_first_call']
plot = params.plot #cfg.PARAMS['plot']
data_path = params.data_path #cfg.PARAMS['data_path']
out_path = params.out_path #cfg.PARAMS['out_path']
dem_file = params.dem_file #cfg.PARAMS['dem_file']
shape_file = params.shape_file #cfg.PARAMS['shape_file']
smooth_window = params.smooth_window #cfg.PARAMS['smooth_window']
single_fl = params.single_fl #cfg.PARAMS['single_fl']

############################################################

# # this is a variable used in the centerline class.
# GAUSSIAN_KERNEL = dict()
# for ks in [5, 7, 9]:
#     kernel = gaussian(ks, 1)
#     GAUSSIAN_KERNEL[ks] = kernel / kernel.sum()


# class defined to be able to use some functions from OGGM "as they are".
class glacier_dir(object):
    def __init__(self, grid):
        self.grid = grid
        self.sourceDem = os.path.join( f'{data_path}', f'{dem_file}')
        self.sourceOutline = os.path.join( f'{data_path}', f'{shape_file}')

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
        clss.append(utils.Centerline(li))
    cls = clss
    return cls


################################################################
# read shapefile (glacier outlines)

def get_centerlines_rgi(dem_path, outline_path):
    """
    Given outline and DEM get glacier centerlines.

    Parameters
    ----------
    dem_path : str
        Gopography. Path to .tif file
    outline_path : str
        Glacier outlines. Path to shapefile.

    Returns
    -------
    None.

    """

    dem = rio.open_rasterio(dem_path)

    # smooth DEM using gaussian smoothing.
    dx = abs(dem.x[0] - dem.x[1])
    gsize = int(np.rint(smooth_window / dx))
    smoothed_dem = gaussian_blur(np.array(dem.values[0]), int(gsize))
    dem.values[0] = smoothed_dem
    
    crop_extent = gpd.read_file(os.path.join(data_path, shape_file))

    crop_extent = crop_extent.to_crs(dem.rio.crs)
    
    # Check that projection is in metre
    try:
        assert crop_extent.crs.axis_info[0].unit_name == 'metre'
              
    except Exception:
        raise Exception('Projection from input shapefile data is not in meters.')
    
    # # Check that projection is in metre
    # try:
    #     assert dem.rio.crs.data['units'] == 'm'
             
    # except Exception:
    #     raise Exception('Projection from input DEM data is not in meters.')
        
    # view all glaciers at once:
    #if plot:
    #    crop_extent.plot()
    
    # Get altitude and pixel info: 
    # altitude, (xorigin, yorigin, pixelH, pixelW)
    data, pix_params = coordinate_change(dem_path)
    
    # select i-th glacier
    crp1 = crop_extent.iloc[[0]]
    
    # crop the DEM to the outline + a few buffer points. Buffer in meters
    dem_clipped = dem.rio.clip(crp1.buffer(20).apply(shpg.mapping), 
                               crop_extent.crs, drop=False)
    
    #dem=dem_clipped
    # get pix paramaters and topography of the new clipped DEM
    #dem_clipped.rio.to_raster('.tmp.tif')
    #data, pix_params = coordinate_change('.tmp.tif')
    
    #del dum_data
    
    # assign some value to outside crop: e.g. -1 (default number is too large)
    dem_clipped.values[0][dem_clipped.values[0] < 0 ] = -1
    dem_clipped.values[0][dem_clipped.values[0] > 8848 ] = -1                           
    
    # Determine heads and tails #
    area = crop_extent.geometry[0].area
    
    
    # list of outline X,Y coordinates
    points_xy = crop_extent.geometry.exterior[0].coords
    
    
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
    #nx, ny = len(dem.x), len(dem.y)
    
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
    glacier_poly_hr = crop_extent.geometry[0]
    
    
    
    # heads' coordinates and altitude
    heads, heads_z = utils._filter_heads(heads, list(heads_z), float(radius),
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
    
    #####
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
       
    #Convert to multiline
    if len(cls_xy) > 1:
        cls_xy = shpg.MultiLineString(cls_xy)
    
    # convert to geopandas datagrame to save afterwards
    cls_xy_gpd = gpd.GeoDataFrame(cls_xy)
    cls_xy_gpd['geometry'] = cls_xy_gpd[0]; del cls_xy_gpd[0]
    
    # add some metadata
    cls_xy_gpd['src_files'] = os.path.join(f'{gdir.sourceDem}',f'{gdir.sourceOutline}')
    
    if single_fl:
        fileout = 'main_flowline' 
    else:
        fileout = 'flowlines' 
     
    # save
    cls_xy_gpd.to_file(os.path.join(f"{out_path}", f"{fileout}.shp"))

    return None
##############################################################
# END #

dem_path = os.path.join(data_path, dem_file)
outline_path = os.path.join(data_path, shape_file)

get_centerlines_rgi(dem_path, outline_path)