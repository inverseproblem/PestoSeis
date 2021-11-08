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
        # create a 2D grid with 50x30 cells
        nx,ny = 100,60
        # size of the cells
        dh = 2.5
        # origin of grid axes
        xinit,yinit = 0.0,0.0

        # create the dictionary containing the grid parameters
        gridpar = tr.setupgrid(nx, ny, dh, xinit, yinit)

        # define a velocity model
        velmod = np.zeros((nx,ny)) + 1500

        # define the position of sources and receivers, e.g.,
        recs = np.array([[30.4, 22.3],
                         [10.1,  20.0],
                         [12.4,  9.5]])
        srcs = np.array([[ 3.4,  2.3],
                         [42.4, 15.5]])
        ## calculate all traveltimes
        ttpick,ttime = tr.traveltime(velmod,gridpar,srcs,recs)

        ## now trace rays (ttime contains a set of 2D traveltime arrays)
        rays = tr.traceallrays(gridpar,srcs,recs,ttime)

        rays_xy = np.zeros([rays.shape[0],rays[0][0]["xy"].shape[0],rays[0][0]["xy"].shape[1]])
        rays_ij = np.zeros([rays.shape[0],rays[0][0]["ij"].shape[0],rays[0][0]["ij"].shape[1]])
        rays_segment_len = np.zeros([rays.shape[0],rays[0][0]["segment_len"].shape[0],rays[0][0]["segment_len"].shape[1]])
        for i in range(rays.shape[0]):
            rays_xy[i,:] = rays[0][0]["xy"]
            rays_ij[i,:] = rays[0][0]["ij"] 
            rays_segment_len[i,:] = rays[0][0]["segment_len"] 
        # np.save("reference_solutions/rays_homogeneous.npy",rays,allow_pickle=True)

        # load reference solution
        ref_rays = np.load("reference_solutions/rays_homogeneous.npy",allow_pickle=True)
        ref_rays_xy = np.zeros([ref_rays.shape[0],ref_rays[0][0]["xy"].shape[0],ref_rays[0][0]["xy"].shape[1]])
        ref_rays_ij = np.zeros([ref_rays.shape[0],ref_rays[0][0]["ij"].shape[0],ref_rays[0][0]["ij"].shape[1]])
        ref_rays_segment_len = np.zeros([ref_rays.shape[0],ref_rays[0][0]["segment_len"].shape[0],ref_rays[0][0]["segment_len"].shape[1]])
        for i in range(rays.shape[0]):
            ref_rays_xy[i,:] = ref_rays[0][0]["xy"]
            ref_rays_ij[i,:] = ref_rays[0][0]["ij"]
            ref_rays_segment_len[i,:] = ref_rays[0][0]["segment_len"]

        if np.allclose(rays_xy,ref_rays_xy, rtol=rtol, atol=atol) != True:
            print("Test failed. Difference between coordinates of rays too large. Please check the code.")
        if np.allclose(rays_ij,ref_rays_ij, rtol=rtol, atol=atol) != True:
            print("Test failed. Difference between indices too large. Please check the code.")
        if np.allclose(rays_segment_len,ref_rays_segment_len, rtol=rtol, atol=atol) != True:
            print("Test failed. Difference between segement lengths too large. Please check the code.")

        atol=1e-03
        rtol=1e-05
        self.assertTrue(np.allclose(rays_xy,ref_rays_xy, rtol=rtol, atol=atol))
        self.assertTrue(np.allclose(rays_ij,ref_rays_ij, rtol=rtol, atol=atol))
        self.assertTrue(np.allclose(rays_segment_len,ref_rays_segment_len, rtol=rtol, atol=atol))
    
        
        
        





