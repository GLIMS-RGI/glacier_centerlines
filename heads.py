"""
Compute glacier heads following 
## https://github.com/OGGM/oggm/blob/master/oggm/core/centerlines.py
and https://tc.copernicus.org/articles/8/503/2014/
https://github.com/OGGM/oggm/blob/447a49d7f936dae4870453d7c65bf2c6f861d0d8/oggm/core/gis.py#L798
"""
import numpy as np
import scipy
from utils import *
import utils
from utils import _filter_heads
import salem
#import shapely
#import copy
import shapely.geometry as shpg
import matplotlib.pyplot as plt
import os
import rioxarray as rio
import geopandas as gpd
#import logging
import copy
try:
    import skimage.draw as skdraw
except ImportError:
    pass

try:
    from skimage import measure
    from skimage.graph import route_through_array
except ImportError:
    pass
#------------------------ Import functions ---------
from functions import (get_terminus_coord, profile, coordinate_change, 
                       _make_costgrid, _polygon_to_pix, _filter_lines,
                       _filter_lines_slope, _join_lines)

#some parameters:
# get radius of the buffer according to Kienholz eq. (1)
q1 = 2/10**6 # 1/m
q2 = 500 #m
rmax = 1000 #m
localmax_window = 500 #In units of dx
flowline_dx = 2
flowline_junction_pix = 3
kbuffer = 2.5
min_slope_flowline_filter = 0
filter_min_slope = True
flowline_height_smooth = 1
#


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

#need to run terminus .py

#toadd: single_fl, for now, single_fl = 0
#todo: discard heads, now no head is discarded

# Module logger
#log = logging.getLogger(__name__)
#logging.basicConfig(filename="test.log", level=logging.my_DEBUG)

single_fl = False

# get altitude and pixel info: 
# altitude, (xorigin, yorigin, pixelH, pixelW)
data, pix_params = coordinate_change(dem_path)

#loop over all geometries
for i in np.arange(len(crop_extent)): 
    # start with one outline, and crop the DEM to the outline + a few grid point
    crp1 = crop_extent.iloc[[i]]
    # crop with buffer. Buffer in meters
    dem_clipped = dem.rio.clip(crp1.buffer(20).apply(shpg.mapping),
                               crop_extent.crs)
    dem_clipped0 = dem_clipped.copy()

    # assign some value to outside crop: e.g. 0 (default number is too large)
    dummy_val = dem_clipped.values[0][dem_clipped.values[0] < 32000].mean()
    dem_clipped.values[0][dem_clipped.values[0] > 32000] = dummy_val

#plot zoutline
    area = crop_extent.geometry[i].area
    points_yx=crop_extent.geometry.exterior[i].coords #list of X,Y coordinates
        
    if points_yx[0] == points_yx[-1]:
        points_yx = points_yx[:-1]
    
    #get profile under glacier outline: distance (arbitrary units) - altitude (m)  
    prof = profile(points_yx, data, pix_params)
    
    # get terminus coordinates and position in or data (index)
    xyterm, ind_term = get_terminus_coord(points_yx, prof[1]) #yx, zoutline

    zoutline = prof[1]
    ext_yx = points_yx
    
    ## here heads star
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
        
    # Size of the half window to use to look for local maxima
    maxorder = np.rint(localmax_window/dx) #np.rint(cfg.PARAMS['localmax_window'] / gdir.grid.dx)
    
    #order: number of points to take at eack side. minumum 5, maximum = maxorder
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
        headsy = np.append(headsy, ext_yx[int(j)][1])
        headsx = np.append(headsx, ext_yx[int(j)][0])
 
    headsx = headsx[1:]
    headsy = headsy[1:]
    heads_z = zoutline[heads_idx]
    
    # careful, the coords are in y, x order!
    heads = [shpg.Point(x, y) for y, x in zip(headsy,
                                              headsx)]
    
    radius = q1 * area + q2 # cfg.PARAMS['q1'] * geom['polygon_area'] + cfg.PARAMS['q2']
    radius = utils.clip_scalar(radius, 0, rmax) #cfg.PARAMS['rmax'])
    radius /= grid.dx #gdir.grid.dx  # in raster coordinates

    # params from default OGGM
    # Grid spacing of a flowline in pixel coordinates
    flowline_dx = 2
    # Number of pixels to arbitrarily remove at junctions
    flowline_junction_pix = 3
    
    # Plus our criteria, quite useful to remove short lines:
    radius += flowline_junction_pix * flowline_dx #cfg.PARAMS['flowline_junction_pix'] * cfg.PARAMS['flowline_dx']
    
    # OK. Filter and see.
    poly_pix = crop_extent.geometry[i]
    
    #for r=500 it does the job in the last glacier
    heads, heads_z = _filter_heads(heads, list(heads_z), float(radius), poly_pix)
    
    #if some head removed:
    headsx=np.zeros(1)
    headsy=np.zeros(1)
    for k in np.arange(len(heads)):
        headsy = np.append(headsy, heads[k].y)
        headsx = np.append(headsx, heads[k].x)
         
    headsx = headsx[1:]
    headsy = headsy[1:]
    #heads_z = zoutline[heads_idx]
    
    
    #TODO:
    # unify dem_clipped0 and dem_clipped
    #
    
    # assign some value to outside crop: e.g. 0 (default number is too large)
    dem_clipped0.values[0][dem_clipped0.values[0] < 1500] = 1
    dem_clipped0.values[0][dem_clipped0.values[0] != 1] = 0
    mask = dem_clipped0
    ext = ext_yx
    z = dem_clipped.values[0]
    
#points = dem_clipped0[0, [y_ind], [x_ind]]
    glacier_poly_hr = crp1.geometry[i]
    proj=utm_proj

# Rounded nearest pix
    glacier_poly_pix = _polygon_to_pix(glacier_poly_hr)

    from functools import partial, wraps

    tuple2int = partial(np.array, dtype=np.int64)

    # Compute the glacier mask (currently: center pixels + touched)
    glacier_mask = mask.values[0]
    glacier_ext = np.zeros((ny, nx), dtype=np.uint8)
    (x, y) = glacier_poly_pix.exterior.xy
    
    for gint in glacier_poly_pix.interiors:
         x, y = tuple2int(gint.xy)
         xx, yy = grid.transform(x,y,crs=utm_proj)
         xx, yy = np.round(xx), np.round(yy)
         xx, yy = xx.astype(int), yy.astype(int)
         glacier_mask[yy, xx] = 0  # on the nunataks
    
    x, y = tuple2int(glacier_poly_pix.exterior.xy)
    
    #project xy to our shapefile grid
    xx, yy = grid.transform(x,y,crs=utm_proj)
    xx, yy = np.round(xx), np.round(yy)
    xx, yy = xx.astype(int), yy.astype(int)
    #glacier_mask[y, x] = 1
    glacier_ext[yy, xx] = 1  
    ext = glacier_ext
    
    #mask = mask.values[0]
    #ext = dem_clipped0.values[0]
    #TODO: do mask in the oggm way
    #
    #
    
    mask = glacier_mask
    costgrid = _make_costgrid(mask, ext, z)

    
#####------------------------- Compute centerlines --------------------

    # Cost array
    #costgrid = _make_costgrid(glacier_mask, glacier_ext, topo)

    # Terminus
    #t_coord = _get_terminus_coord(gdir, ext_yx, zoutline)
    t_coord =  xyterm
    
    # Compute the routes
    lines = []
    heads_pix = []
    for h in heads:
        h_coord_pix = np.array(grid.transform(h.x, h.y, crs=utm_proj))
        h_coord_pix = np.round(h_coord_pix).astype(int)
        t_coord_pix = np.array(grid.transform(t_coord.x, t_coord.y, crs=utm_proj))
        t_coord_pix = np.round(t_coord_pix).astype(int)
        heads_pix.append(shpg.Point(h_coord_pix))
        indices, _ = route_through_array(costgrid, np.roll(h_coord_pix,1), np.roll(t_coord_pix,1))
        lines.append(shpg.LineString(np.array(indices)[:, [1, 0]]))
        
    # Filter the shortest lines out
    dx_cls = flowline_dx #cfg.PARAMS['flowline_dx']
    radius = flowline_junction_pix * dx_cls #cfg.PARAMS['flowline_junction_pix'] * dx_cls
    radius += 6 * dx_cls
    
    #heads in raster coordinates:
    
    
    olines, oheads = _filter_lines(lines, heads_pix, kbuffer, radius) #cfg.PARAMS['kbuffer'], radius)
    #log.debug('(%s) number of heads after lines filter: %d',
    #          gdir.rgi_id, len(olines))

    # Filter the lines which are going up instead of down
    min_slope = np.deg2rad(min_slope_flowline_filter)
    do_filter_slope = filter_min_slope


    ##########
    #TODO: bring this up in the code
    #
    #
    
#    class gdir(object):
#        def __init__(self, grid):
#            self.grid = grid
#    
#    class grid_inf(object):
#        def __init__(self,dx,dy,nx,ny):
#            self.dx = dx
#            self.dy = dy
#            self.nx = nx
#            self.ny = ny
            
    class glacier_dir(object):
        def __init__(self, grid):
            self.grid = grid
    
    class grid_inf(object):
        def __init__(self, grid):
            self.grid = grid

    grid = salem.Grid(proj=utm_proj, nxny=(nx, ny), dxdy=(dx, -dx), x0y0=x0y0) 
    
    gdir=glacier_dir(grid)
    
    if do_filter_slope:
        topo = z
        olines, oheads = _filter_lines_slope(olines, oheads, topo,
                                             gdir, min_slope)
       # log.debug('(%s) number of heads after slope filter: %d',
       #           gdir.rgi_id, len(olines))

    #TODO: problem here with the Centerline class
    #
    #
    # And rejoin the cut tails
    #olines = _join_lines(olines, oheads)

    # Adds the line level
    # for cl in olines:
    #     cl.order = line_order(cl)

    # # And sort them per order !!! several downstream tasks  rely on this
    # cls = []
    # for i in np.argsort([cl.order for cl in olines]):
    #     cls.append(olines[i])

    # # Final check
    # if len(cls) == 0:
    #     raise GeometryError('({}) no valid centerline could be '
    #                         'found!'.format(gdir.rgi_id))

    # # Write the data
    # gdir.write_pickle(cls, 'centerlines')

    # if is_first_call:
    #     # For diagnostics of filtered centerlines
    #     gdir.add_to_diagnostics('n_orig_centerlines', len(cls))
        
    #plot profile + terminus + heads:
    if plot:
        #profile
        # NOTE: in this plot, removed heads are displayed anyway.
        # plt.plot(prof[0], zoutline, '-+') #horizontal distance vs altitude
        # plt.plot(prof[0][ind_term], zoutline[ind_term], 'r*', label="terminus") #terminus
        # plt.plot(prof[0][heads_idx], zoutline[heads_idx], 'g*', label="head") #head
        # plt.xlabel("Distance along outline (a.u.)")
        # plt.ylabel("Altitude (m)")
        # plt.legend()
        # plt.show()
            
        #raster
        f, ax = plt.subplots(figsize=(8, 10))
        dem_clipped.plot(ax=ax)
        ax.set(title="Raster Layer Cropped to Geodataframe Extent")
        plt.scatter(headsx,headsy, marker="*",s=1000, c="g")
        plt.scatter(xyterm.x,xyterm.y, marker="*", s=1000, c="r")  
        crp1.boundary.plot(ax=ax)
        plt.show()
        
        #costgrid
        plt.imshow(costgrid)
        plt.colorbar()

        plt.scatter(heads_pix[0].xy[0], heads_pix[0].xy[1], marker="*",s=100, c="g")
        if len(lines) > 1 :
            plt.scatter(heads_pix[1].xy[0], heads_pix[1].xy[1], marker="*",s=100, c="g")
        plt.scatter(t_coord_pix[0], t_coord_pix[1], marker="*", s=100, c="r") 
        
        #plt.Line2D(lines[0].xy[1], lines[0].xy[0])
        plt.scatter(lines[0].xy[0], lines[0].xy[1], marker="o", s=10, c="y")
        if len(lines) > 1 :
            plt.scatter(lines[1].xy[0], lines[1].xy[1], marker="o", s=10, c="y") 

        #crp1.boundary.plot(ax=ax)
        plt.show()
        
        
#        costgrid.plot(ax=ax)
#        ax.set(title="Raster Layer Cropped to Geodataframe Extent")
#        plt.scatter(headsx,headsy, marker="*",s=1000, c="g")
#        plt.scatter(xyterm.x,xyterm.y, marker="*", s=1000, c="r")  
        