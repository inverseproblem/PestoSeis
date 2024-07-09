#!/usr/bin/env python
# coding: utf-8
# %%

# %%

import os
import numpy as np
import matplotlib.pyplot as plt
import sys
import h5py
import pestoseis.ttimerays as tr
import time
import unittest


# %%


class TestTraveltime(unittest.TestCase):

    def test_numerical(self):
        def analyticalsollingrad2D(gridpar,xsrcpos,ysrcpos):
            ##################################
            ## linear gradient of velocity

            ## source must be *on* the top surface
            #@show ysrcpos
            assert (ysrcpos == 0.0)

            # position of the grid nodes
            xgridpos = np.array([i * gridpar["dh"] + gridpar["xinit"] for i in range(gridpar["nx"])])
            ygridpos = np.array([i * gridpar["dh"] + gridpar["yinit"] for i in range(gridpar["ny"])])
            vel0 = 2.0
            # gradient of velociy
            gr = 0.05
            ## construct the 2D velocity model
            velmod = np.zeros([gridpar["nx"],gridpar["ny"]]) 
            for i in range(gridpar["nx"]):
                velmod[i,:] = vel0 + gr * (ygridpos - ysrcpos)

            # https://pubs.geoscienceworld.org/books/book/1011/lessons-in-seismic-computing
            # Lesson No. 41: Linear Distribution of Velocity—V. The Wave-Fronts 

            ## Analytic solution
            ansol2d  = np.zeros([gridpar["nx"],gridpar["ny"]])
            for j in range(gridpar["ny"]):
                for i in range(gridpar["nx"]):
                    x = xgridpos[i]-xsrcpos
                    h = ygridpos[j]-ysrcpos
                    ansol2d[i,j] = 1 / gr * np.arccosh( 1 + (gr ** 2 * (x ** 2 + h ** 2)) / (2 * vel0 * velmod[i,j]))
            return ansol2d,velmod
        
        nx = 100
        ny = 100
        h = 100
        xinit = 0.0
        yinit = 0.0
        hgrid = nx / h
        gridpar = tr.setupgrid(nx,ny,hgrid,xinit,yinit)
        nrec = 100
        x_rec = np.linspace(0.0,gridpar["nx"],nrec)
        y_rec = np.zeros(nrec)
        receivers = np.vstack((x_rec,y_rec))

        # one source on top of domain
        nsrc = 1
        xsrcpos = 20.0
        ysrcpos = 0.0
        srccoo = np.zeros([nsrc,2])
        srccoo[:,0] = xsrcpos;
        srccoo[:,1] = ysrcpos;

        ansol2d,velmod_analytical = analyticalsollingrad2D(gridpar, xsrcpos, ysrcpos)
        ttpick,ttime = tr.traveltime(velmod_analytical, gridpar, srccoo, receivers)
        
        # compute difference
        diff = np.abs(ansol2d - ttime[0][:-1,:-1]) 
        # check that difference is within tolerance of reference solution
        current_dir = os.path.dirname(os.path.abspath(__file__))
        reference_solutions_dir = os.path.join(current_dir, "reference_solutions")
        ref_diff_path = os.path.join(reference_solutions_dir, "diff_analytical_numerical_linearGradVelocity.npy")
        ref_diff = np.load(ref_diff_path, allow_pickle=True)
        if np.allclose(diff, ref_diff, rtol=1e-01, atol=0) != True:
            print("Test failed. Difference between analytical traveltimes and computed traveltimes greater than threshold w.r.t. reference difference. Please check the code")
        self.assertTrue(np.allclose(diff, ref_diff, rtol=1e-01, atol=0))
        
        
        
        
        





