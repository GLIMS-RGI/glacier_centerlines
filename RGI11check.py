#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 21 12:15:39 2022

@author: francesc
"""
# get DEM for RGIglaciers
# https://oggm.org/tutorials/master/notebooks/dem_sources.html#use-pre-processed-dems-from-rgi-topo
# use NASADEM, the default DEM for low and mid-latitudes in OGGM, you can also change to e.g. 'COPDEM'
import rioxarray as rioxr
import matplotlib.pyplot as plt
import oggm
from oggm import cfg, utils, workflow, tasks, graphics
from oggm.core import gis
import numpy as np
import geopandas as gpd
import os
cfg.initialize(logging_level='WARNING')
cfg.PARAMS['border'] = 10

from oggm.shop import rgitopo
cfg.PATHS['working_dir'] = utils.gettempdir('rgitopo', reset=True)


################################
# current oggm_flowlines
################################
shape_file = '/home/francesc/data/glacier_centerlines/RGI11_centerlines/RGI11_centerlines.shp'

# find length glaciers region RGI11:
oggm_flowlines = gpd.read_file(os.path.join("", shape_file))
RGIId = oggm_flowlines.RGIID

length = oggm_flowlines.geometry.length
new_len = []

for ii in RGIId:
    if sum(oggm_flowlines['RGIID'] == ii) !=1: #glacier with multiple branches, we take the longest:
        #print(sum(crop_extent['RGIID'] == ii)) #number of branches
        ilen = np.max(length[oggm_flowlines['RGIID'] == ii])
        new_len.append(ilen)
    else:
        #print(sum(crop_extent['RGIID'] == ii)) #if only one branch, just its length: 
        new_len.append(float(length[oggm_flowlines['RGIID'] == ii]))

final_len = list(set(new_len)) # glacier length over all rgi11
ini_len = sum(length) # glacier length of all lines (main flowlines + tributaries)

# Length oggm centerlines:
print(f'the length og OGGM centerlines is {sum(new_len)}') 


######################################
# run glacier_centerlines.main for all RGI11:
######################################   
#tif_path = os.path.join(data_path,'RGI60-11.00.tar/RGI60-11.00999.tar.gz/NADASEM/dem.tif')
#'/home/francesc/data/OGGM/rgi/RGIV60/11_rgi60_CentralEurope/RGI60-11.00/RGI60-11.00001/NASADEM/dem.tif'

# DEM stored in 
# /home/francesc/data/OGGM/download_cache/cluster.klima.uni-bremen.de/data/gdirs/dems_v1/default/RGI62/b_010/L1/RGI60-11
for i in np.arange(3920, 3927+1):
    idrgi='RGI60-11.0' + str(i)
    
