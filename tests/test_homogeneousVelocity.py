#!/usr/bin/env python
# coding: utf-8

# In[4]:


import numpy as np
import matplotlib.pyplot as plt
import sys
import h5py
sys.path.append("/Users/inesulrich/Documents/PestoSeis/pestoseis")
import ttimerays as tr
import time
import unittest


# In[5]:


class TestTraveltime(unittest.TestCase):

    def test_numerical(self):
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

        # one source inside domain
        nsrc = 1
        xsrcpos = 50.0
        ysrcpos = 50.0
        srccoo = np.zeros([nsrc,2])
        srccoo[:,0] = xsrcpos;
        srccoo[:,1] = ysrcpos;

        # define a velocity model
        velmod = np.zeros((nx,ny)) + 1500

        ## Analytic solution
        # position of the grid nodes
        xgridpos = np.array([i * gridpar["dh"] + gridpar["xinit"] for i in range(gridpar["nx"])])
        ygridpos = np.array([i * gridpar["dh"] + gridpar["yinit"] for i in range(gridpar["ny"])])
        ansol2d  = np.zeros([gridpar["nx"],gridpar["ny"]])
        for j in range(gridpar["ny"]):
            for i in range(gridpar["nx"]):
                x = xgridpos[i]-xsrcpos
                h = ygridpos[j]-ysrcpos
                ansol2d[i,j] = np.sqrt(x ** 2 + h ** 2) / 1500
        ttpick,ttime = tr.traveltime(velmod, gridpar, srccoo, receivers)
        
        atol=1e-03
        rtol=1e-05
        if np.allclose(ansol2d,ttime[0][:-1,:-1], rtol=rtol, atol=atol) != True:
            print("Test failed. Difference between analytical traveltimes and computed traveltimes greater than threshold.")
        self.assertTrue(np.allclose(ansol2d,ttime[0][:-1,:-1],rtol=rtol, atol=atol))
        
        
        
        





