#!/usr/bin/env python
# coding: utf-8
# %%

# %%


import numpy as np
import matplotlib.pyplot as plt
import sys
import h5py
sys.path.append("/Users/inesulrich/Documents/PestoSeis/pestoseis")
import ttimerays as tr
import time
import unittest


# %%


class TestTraveltime(unittest.TestCase):
    
    def test_numerical(self):
        # number of layers
        Nlay = 120
        # depth of layers -- includes both top and bottom (Nlay+1)
        laydepth = np.linspace(0.0,2000.0,Nlay+1)[1:]
        # velocity
        velmod = np.linspace(2000.0,3000.0,Nlay)
        # origin of ray
        xystart = np.array([0.0, 0.0])
        # take off angle
        takeoffangle = 45.0

        # trace a single ray
        raypath,tt,dist = tr.tracerayhorlay(laydepth, velmod, xystart, takeoffangle)
        
        
        raypath_ref = np.load("reference_solutions/raypath_horizLayeredModel.npy") 
        tt_ref = np.load("reference_solutions/tt_horizLayeredModel.npy") 
        dist_ref = np.load("reference_solutions/dist_horizLayeredModel.npy")
        
        atol=1e-06
        rtol=1e-05
        if np.allclose(raypath,raypath_ref, rtol=rtol, atol=atol) != True:
            print("Test failed. Difference in raypath too large. Please check the code.")
        if np.allclose(tt,tt_ref, rtol=rtol, atol=atol) != True:
            print("Test failed. Difference in traveltime too large. Please check the code.")
        if np.allclose(dist,dist_ref, rtol=rtol, atol=atol) != True:
            print("Test failed. Difference in distance too large. Please check the code.")

        
        
        
        





