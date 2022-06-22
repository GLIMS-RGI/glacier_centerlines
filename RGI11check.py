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
cfg.initialize(logging_level='WARNING')
cfg.PARAMS['border'] = 10

from oggm.shop import rgitopo
cfg.PATHS['working_dir'] = utils.gettempdir('rgitopo', reset=True)

# DEM stored in 
# /home/francesc/data/OGGM/download_cache/cluster.klima.uni-bremen.de/data/gdirs/dems_v1/default/RGI62/b_010/L1/RGI60-11
for i in np.arange(3920, 3927+1):
    idrgi='RGI60-11.0' + str(i)
    
shape_file = '/home/francesc/data/glacier_centerlines/RGI11_centerlines/RGI11_centerlines.shp'
# find length glaciers region RGI11:
crop_extent = gpd.read_file(os.path.join(data_path, shape_file))
RGIId = crop_extent.RGIID

length=0
for ii in RGIId[0:3]:
    a=0
    # lines_per_glacier = sum(crop_extent['RGIID'] == ii)
    # iloc = np.argmax(crop_extent['RGIID'] == ii)

    # if lines_per_glacier != 1:
    #     # take the longest
    #     dumlen=[]
    #     lengths = [ dumlen.append(crop_extent.iloc[k].geometry.length)\
    #                for k in np.arange(lines_per_glacier)]
    #     #length = np.max(lengths)
    #     print(lengths)
    # else:
    #     crop_extent.iloc[iloc].geometry.length
    #     length += crop_extent.iloc[iloc].geometry.length
    #     print(iloc)
#print(length)

length = crop_extent.geometry.length
length = length[0:5]
new_len = []

for ii in RGIId[0:5]:
    if sum(crop_extent['RGIID'] == ii) !=1:
        ilen = np.max(length[crop_extent['RGIID'] == ii])
        new_len.append(float(ilen))
    else:
        new_len.append(length[crop_extent['RGIID'] == ii])

sum(list(set(new_len)))
sum(length)

    
length = crop_extent.geometry.length

length = []
for ii in RGIId:
    length.append(crop_extent.geometry.length[crop_extent['RGIID'] == ii])    
