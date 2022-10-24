#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 21 16:48:31 2022

Testing functions that are not part of OGGM or modified from it.

@author: francesc
"""
import numpy as np
from osgeo import gdal
from functions import coordinate_change, profile

def test_coordinate_change():
    values = coordinate_change('./test_data/dem.tif')[0][1]
    pix_params = coordinate_change('./test_data/dem.tif')[1]
    assert all(values == [2157, 2157, 2156, 2156, 2156, 2157, 2159, 2160, 2161, 2162, 2164,
           2167, 2169, 2171, 2172, 2174, 2175, 2177, 2180, 2184, 2187, 2190,
           2194, 2197, 2202, 2206, 2209, 2213, 2216, 2219, 2221, 2223, 2224,
           2228, 2234, 2241, 2249, 2259, 2270, 2286, 2305, 2325, 2345, 2364,
           2381, 2398, 2416, 2437, 2461, 2480, 2494])
    assert pix_params == [-400.26259317302896, 5260483.663235244, 15.0, 15.0]
    
def test_profile():
    points_xy = [(0, 1),
             (2, 1),
             (2, 0),
             (0, 0)]
    data=np.array([[1,2,3],[5,4,3]])

    pix_params = [0, 0, 1, 1]
    prof = profile(points_xy, data, pix_params)
    
    distance = prof[0]
    altitude = prof[1]
    
    assert all(distance == [2., 3., 5., 6.])
    assert all(altitude == [5., 3., 3., 1.])