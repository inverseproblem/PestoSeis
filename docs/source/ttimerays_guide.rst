

.. _ttimerays_guide:

*******************************************
Traveltimes and rays -- using ``ttimerays``
*******************************************

====================================
Setup of grid parameters and models
====================================

In order to compute traveltimes and rays, a 2D grid must be set up by defining its dimensions in terms of number of cells in the two ``x`` and ``y`` directions and cell size.
See :func:`pestoseis.ttimerays.setupgrid` for how to create the *dictionary* holding the grid parameters needed for subsequent computations.

Example::
  
  import sys
  # add location of ttimerays package to path
  sys.path.append('pyteachseismo/')
  import pestoseis.ttimerays as TR
  # import the plotting library
  import matplotlib.pyplot as PL

  # create a 2D grid with 50x30 cells
  nx,ny = 50,30
  # size of the cells
  dh = 5.0
  # origin of grid axes
  xinit,yinit = 0.0,0.0

  # create the dictionary containing the grid parameters
  gridpar = TR.setupgrid(nx, ny, dh, xinit, yinit)

==================
 Visualization
==================

Some functions for different kinds of plots are provided (click on the function
name to get the docstring):

* plot the grid: :func:`pestoseis.ttimerays.plotgrid` 
* plot the traveltime array: :func:`pestoseis.ttimerays.plotttimemod` 
* plot the velocity model: :func:`pestoseis.ttimerays.plotvelmod`
* plot (previously traced) rays: :func:`pestoseis.ttimerays.plotrays`  


Example of plotting the grid::

  import sys
  # add location of ttimerays package to path
  sys.path.append('pyteachseismo/')
  import pestoseis.ttimerays as TR
  # import the plotting library
  import matplotlib.pyplot as PL

  # create a 2D grid with 50x30 cells
  nx,ny = 50,30
  # size of the cells
  dh = 5.5
  # origin of grid axes
  xinit,yinit = 0.0,0.0

  # create the dictionary containing the grid parameters
  gridpar = TR.setupgrid(nx, ny, dh, xinit, yinit)
  
  # plot the grid
  PL.figure()
  TR.plotgrid(gridpar)
  PL.show()

which produces the following image

.. figure::  images/gridpl.png
   :align:   center
   :width: 400px


==================
Traveltimes
==================

Traveltime calculation given a velocity model and one or more sources and related receivers can be performed using the function :func:`pestoseis.ttimerays.traveltime`. By default the function returns both the traveltimes at the receivers and also the entire 2D traveltime array(s) for subsequent ray tracing.
Example::

  [...]
  # define a grid
  [...]
  # define a velocity model
  velmod = 3.0*ones(nx,ny)
  # define the position of sources and receivers, e.g.,
  recs = NP.array([[30.4, 22.3],
                   [12.4,  9.5]])
  srcs = NP.array([[ 3.4,  2.3],
                   [42.4, 15.5]])
  ## calculate all traveltimes
  ttpick,ttime = TR.traveltime(velmod,gridpar,sources,receivers)


==================  
Rays 
==================

-------------------------------------------
Trace rays in a 2D heterogeneus model
-------------------------------------------

In order to trace rays in a 2D model, the traveltime arrays are needed first. So, after computing traveltimes, it is possible to trace (approximatively) the rays using the function :func:`pestoseis.ttimerays.traceallrays`

Example::

  [...]
  ## compute traveltimes
  ttpick,ttime = TR.traveltime(velmod,gridpar,sources,receivers)
  ## now trace rays (ttime contains a set of 2D traveltime arrays)
  rays = TR.traceallrays(gridpar,sources,receivers,ttime)

-------------------------------------------  
Trace *straight* rays 
-------------------------------------------

The function :func:`pestoseis.ttimerays.traceall_straight_rays` traces rays simply as a straight lines between source and receiver.

Example::

  [...]
  ## now trace straight rays 
  rays = TR.traceall_straight_rays(gridpar,sources,receivers)

  
-----------------------------------------------
Trace rays in a *horizontally layered* medium
-----------------------------------------------

The function :func:`pestoseis.ttimerays.tracerayhorlay` provides a way to compute ray paths, traveltime and distance covered in a horizontally layered medium, provided depths of layers and their velocity. The geometrical setup is as following::

    #
    # v1,theta1
    #  
    # --------xpt1------------
    #           \
    # v2,theta2  \
    #             \  
    # ------------xpt2--------
    #        
    #  
    #   |\
    #   | \   theta
    #   |  \
    #   |-->\
    #   |    *
    #
    
The angle *theta* is the take off angle, measured anti-clockwise from the vertical.
Example::

  [...]
  import numpy as NP
  # number of layers
  Nlay = 120
  # depth of layers -- includes both top and bottom (Nlay+1)
  laydepth = NP.linspace(0.0,2000.0,Nlay+1)[1:]
  # velocity
  vp = NP.linspace(2000.0,3000.0,Nlay)
  # origin of ray
  xystart = NP.array([0.0, 0.0])
  # take off angle
  takeoffangle = 45.0
  
  # trace a single ray
  TR.tracerayhorlay(laydep, vel, xystart, takeoffangle)

-----------------
Ray tomography
-----------------

A function to perform simple linear inversion under Gaussian assumptions (least squares approach) is provided in :func:`pestoseis.ttimerays.lininv`. In order to run the inversion the `tomography matrix` (containing the length of the rays in each cell), the prior mean model and covariances for observed data and model parameters are needed. The rays can be calcutated using :func:`ttimerays.traceallrays` and subsequently the tomography matrix can be built using :func:`pestoseis.ttimerays.buildtomomat`. This kind of inversion is quite primitive and therefore often unstable. The result are the posterior mean model and covariance matrix (we are under Gaussian assumptions).

.. math::
   
   S( \mathbf{m}) = \frac{1}{2} ( \mathbf{G} \mathbf{m} - \mathbf{d}_{\sf{obs}} )^{\sf{T}}
   \mathbf{C}^{-1}_{\rm{D}} ( \mathbf{G} \mathbf{m} - \mathbf{d}_{\sf{obs}} ) 
   + \frac{1}{2} ( \mathbf{m} - \mathbf{m}_{\sf{prior}} )^{\sf{T}} \mathbf{C}^{-1}_{\rm{M}}
     ( \mathbf{m} - \mathbf{m}_{\sf{prior}} ) 
 

The posterior covariance matrix is given by 

.. math::
   
   \mathbf{\widetilde{C}}_{\rm{M}} =  \left( \mathbf{G}^{\sf{T}} \,
   \mathbf{C}^{-1}_{\rm{D}} \, \mathbf{G} + \mathbf{C}^{-1}_{\rm{M}} \right)^{-1}

and the center of posterior Gaussian (the mean model) is 

.. math::  

   \mathbf{\widetilde{m}}  
   = \mathbf{m}_{\rm{prior}}+ \mathbf{\widetilde{C}}_{\rm{M}} \, \mathbf{G}^{\sf{T}} \, \mathbf{C}^{-1}_{\rm{D}} \left(\mathbf{d}_{\rm{obs}} - \mathbf{G} \mathbf{m}_{\rm{prior}} \right) .



Example::

  [...]

  # trace rays
  rays = TR.traceallrays(gridpar,sources,receivers,bkgttimegrd)
  # build the tomography matrix 
  tomomat,residualsvector = TR.buildtomomat(gridpar, rays, residuals)
  
  # Perform the actual inversion using a "least-squares" approach
  postm,postC_m = TR.lininv(tomomat,cov_m,cov_d,mprior,residualsvector)



======================================
API ``ttimerays``, list of functions
======================================

.. automodule:: pestoseis.ttimerays
   :members:
   :imported-members:  
.. Need to use :imported-members: since ttimerays imports
   the functions in __init__.py, so has no own members.
