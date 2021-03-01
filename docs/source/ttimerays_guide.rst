.. role:: raw-math(raw)
    :format: latex html
.. _ttimerays_guide:

*******************************************
Traveltimes and rays -- using ``ttimerays``
*******************************************

One possible way to model wave propagation in a medium is to assume that waves can be approximated by rays of infinite frequency along the path between a source and a receiver. If we consider a specific ray :math:`i` along a path :math:`\Gamma_i`, then we can obtain the travel time :math:`t_i` belonging to that ray by solving the line integral

.. math::

   t_i=\int_{\Gamma_i(s(\mathbf{x}))}s(\mathbf{x}(l))dl,
    
where :math:`s=s(\mathbf{x})` is the slowness map of the medium and is related to the speed of sound by :math:`s(\mathbf{x})=\frac{1}{c(\mathbf{x})}`, :math:`dl` is an infinitesimal line segment on the path and :math:`\mathbf{x(}l)` is the parametrization of the spatial variable in terms of :math:`l`. To solve the continuous line integral, one commonly discretizes the domain by approximating the continuous slowness map with a finite grid in which the slowness within each cell is constant. The travel times of an entire ray :math:`i` between two points is then given by the line connecting them and is calculated by the sum of the cell segments as

.. math::   

   t_i=\sum_{j=1}^nl_{ij}s_j,
	
where :math:`l_{ij}` is the line of ray :math:`i` in cell :math:`j` and :math:`n` is the total number of cells. In general, the path itself depends on the slowness structure of the medium, hence the equation is non-linear in the path. 

One common approach to simplify the problem is to linearize about some reference model and to assume that variations in the speed of sound are small, which means that the speed-of-sound map can be approximated by :math:`s(\mathbf{x})\approx 1/c_0`. With this, rays are fixed to straight lines between a source and a receiver and the path becomes independent of the slowness distribution of the medium. This simplifies the line integral to

.. math::
   
   t_i=\int_{\Gamma_i}s(\mathbf{x}(l))dl.

With the straight-ray approximation we model wave propagation by performing the simplest ray tracing approach for which Snell’s law at all cell boundaries is ignored so that waves traveling from sources to receivers are approximated by straight lines. If we imagine to perform an experiment using a grid of :math:`n` cells and a total number of :math:`m` source-receiver pairs, resulting in :math:`m` rays, then we can leverage the benefit of the linear forward problem by building a linear system of equations 

.. math::

   \begin{eqnarray}
      \begin{gathered}
         t_1=l_{11}s_1+\dots+l_{1j}s_{j}+\dots+l_{1n}s_n \\ 
         \vdots \\
         t_i=l_{i1}s_1+\dots+l_{ij}s_{j}+\dots+l_{in}s_n \\ 
         \vdots \\
         t_m=l_{m1}s_1+\dots+l_{mj}s_{j}+\dots+l_{mn}s_n. 
      \end{gathered}
   \end{eqnarray}

This can be condensed to matrix vector notation by introducing the forward modelling matrix :math:`\mathbf{F}` of dimensions :math:`m\cross n` that collects all line segments :math:`l_{ij}` for every source-receiver pair as

.. math::

      \mathbf{t}=\mathbf{F}\mathbf{s},

where :math:`\mathbf{t}` is the vector of travel times from every source to every receiver and :math:`\mathbf{s}` is a vector containing the slowness map. 

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

The function :func:`pestoseis.ttimerays.tracerayhorlay` provides a way to compute ray paths, traveltime and distance covered in a horizontally layered medium, provided depths of layers and their velocity. The geometrical setup is as following:

.. figure::  images/geom-rays-horlayers.png
   :align:   center
   :width: 300px
    
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
