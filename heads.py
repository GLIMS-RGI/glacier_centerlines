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
import shapely
import copy
import shapely.geometry as shpg
import matplotlib.pyplot as plt

#need to run terminus .py

#toadd: single_fl, for now, single_fl = 0
#todo: discard heads, now no head is discarded

single_fl = False

#loop over all geometries
for i in np.arange(len(crop_extent)): 
    # Crop the DEM to the outline + a few grid point
    crp1 = crop_extent.iloc[[i]]
    ## crop with buffer (in meters)
    dem_clipped = dem.rio.clip(crp1.buffer(20).apply(shpg.mapping),
                               crop_extent.crs)

    # assign some value to outside crop: e.g. 0 (default number is too large)
    dummy_val = dem_clipped.values[0][dem_clipped.values[0] < 1500].mean()
    dem_clipped.values[0][dem_clipped.values[0] > 1500] = 0 #dummy_val
    
    # plot DEM for each glacier
    if plot:
        f, ax = plt.subplots(figsize=(8, 10))
        dem_clipped.plot(ax=ax)
        ax.set(title="Raster Layer Cropped to Geodataframe Extent")
        plt.show()

    ## work with outlines
    area = crop_extent.geometry[i].area
    points_yx=crop_extent.geometry.exterior[i].coords #list of X,Y coordinates
    prof = profile(points_yx)

    zoutline = prof[1][:-1]
    ext_yx=points_yx[:-1]

    xyterm, ind_term = get_terminus_coord(ext_yx, zoutline)
    
    # plot profile (zoutline) for each glacier
    if plot:
        #plt.plot(prof[0], prof[1]) 
        plt.plot(prof[0][:-1], zoutline, '-') #horizontal distance vs altitude
        plt.plot(prof[0][:-1][ind_term], zoutline[ind_term], 'r*') #terminus
    
        plt.show()  
    
    ## here heads start
    # create grid
    nx = len(dem_clipped.x)
    ny = len(dem_clipped.y)
    dx = abs(max(dem_clipped.x)-min(dem_clipped.x))/nx
    dy = abs(max(dem_clipped.y)-min(dem_clipped.y))/ny
    ulx = min(dem_clipped.x)
    uly = max(dem_clipped.y)
    x0y0 = (ulx+dx/2, uly-dx/2)  # To pixel center coordinates
    
    # try the curent crs
    utm_proj = salem.check_crs(crop_extent.crs)
        
    grid = salem.Grid(proj=utm_proj, nxny=(nx, ny), dxdy=(dx, -dx), x0y0=x0y0)
    
    localmax_window = 10 #I don't know units. It has to be in units of dx
    
    # Size of the half window to use to look for local maximas
    maxorder = np.rint(localmax_window/dx) #np.rint(cfg.PARAMS['localmax_window'] / gdir.grid.dx)
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
    for i in heads_idx:
        headsx = np.append(headsx, ext_yx[int(i)][1])
        headsy = np.append(headsy, ext_yx[int(i)][0])
    headsx = headsx[1:]
    headsy = headsy[1:]
    heads_z = zoutline[heads_idx]
    
    # careful, the coords are in y, x order!
    heads = [shpg.Point(x, y) for y, x in zip(headsy,
                                              headsx)]
    
    # get radius of the buffer according to Kienholz eq. (1)
    q1 = 2/10**6 # 1/m
    q2=500 #m
    rmax = 1000 #m
    
    radius = q1*area + q2# cfg.PARAMS['q1'] * geom['polygon_area'] + cfg.PARAMS['q2']
    radius = utils.clip_scalar(radius, 0, rmax)#cfg.PARAMS['rmax'])
    radius /= dx #gdir.grid.dx  # in raster coordinates
    
    # Plus our criteria, quite useful to remove short lines:
    #radius += cfg.PARAMS['flowline_junction_pix'] * cfg.PARAMS['flowline_dx']
    #log.debug('(%s) radius in raster coordinates: %.2f',
    #          gdir.rgi_id, radius)
    
    # OK. Filter and see.
    #log.debug('(%s) number of heads before radius filter: %d',
    #          gdir.rgi_id, len(heads))
    #poly_pix = crop_extent.geometry[i]
    #heads, heads_z = _filter_heads(heads, list(heads_z), radius, poly_pix)
    #log.debug('(%s) number of heads after radius filter: %d',
    #          gdir.rgi_id, len(heads))

    # plot profile + terminus + heads:
    plt.plot(prof[0][:-1], zoutline, '-') #horizontal distance vs altitude
    plt.plot(prof[0][:-1][ind_term], zoutline[ind_term], 'r*', label="terminus") #terminus
    plt.plot(prof[0][:-1][heads_idx], zoutline[heads_idx], 'g*', label="head") #head
    plt.xlabel("Distance along outline (n.u.)")
    plt.ylabel("Altitude (m)")
    plt.legend()
    plt.show()
    
        
