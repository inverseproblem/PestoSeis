

.. _userguide:

User guide
=============


Installation --- loading the module
-------------------------------------

.. warning:: Installing the package globally with, for instance, ``pip`` is still untested.
	     
A simple way to load the package/module from a local folder is to start a script with adding to the path the location the folder containing the package. This way Python will be able to find the module and therefore load it.

If we currently are in the directory ``dir1`` and the package is located in the directory ``pyteachseismo`` ::

  dir1/
     myscript.py
     
     pyteachseismo/
         ttimerays/
	    ttimerays.py
	    ...
	
Then something along the along the following lines should work: ::

  import sys
  # add location of ttimerays package to path
  sys.path.append('pyteachseismo/')
  import ttimerays as TR

At this point, the functions from ``ttimerays`` should be available as, for instance::

  TR.traceallrays(gridpar, srccoo, reccoo, grdsttime)


Using ``ttimerays``
---------------------

Setup of grid parameters and models
++++++++++++++++++++++++++++++++++++

In order to compute traveltimes and rays, a 2D grid must be set up by defining its dimensions in terms of number of cells in the two ``x`` and ``y`` directions and cell size.
See :func:`ttimerays.setupgrid` for how to create the *dictionary* holding the grid parameters needed for subsequent computations.

Example::
  
  import sys
  # add location of ttimerays package to path
  sys.path.append('pyteachseismo/')
  import ttimerays as TR
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

Visualization
++++++++++++++++++

Some functions for different kinds of plots are provided (click on the function
name to get the docstring):

* plot the grid: :func:`ttimerays.plotgrid` 
* plot the traveltime array: :func:`ttimerays.plotttimemod` 
* plot the velocity model: :func:`ttimerays.plotvelmod`
* plot (previously traced) rays: :func:`ttimerays.plotrays`  


Example of plotting the grid::

  import sys
  # add location of ttimerays package to path
  sys.path.append('pyteachseismo/')
  import ttimerays as TR
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

which produces

.. figure::  images/gridpl.png
   :align:   center



Traveltimes
++++++++++++++++

Traveltime calculation given a velocity model and one or more sources and related receivers can be performed using the function :func:`ttimerays.traveltime`. By default the function returns both the traveltimes at the receivers and also the entire 2D traveltime array(s) for subsequent ray tracing.
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


Rays 
++++++++++++++++

In order to trace rays, the traveltime arrays are needed first. So, after computing traveltimes, it is possible to trace (approximatively) the rays using the function :func:`ttimerays.traceallrays`

Example::

  [...]
  ## compute traveltimes
  ttpick,ttime = TR.traveltime(velmod,gridpar,sources,receivers)
  ## now trace rays (ttime contains a set of 2D traveltime arrays)
  rays = TR.traceallrays(gridpar,sources,receivers,ttime)



Ray tomography
++++++++++++++++


Example::

  [...]

  # trace rays
  rays = TR.traceallrays(gridpar,sources,receivers,bkgttimegrd)
  # build the tomography matrix 
  tomomat,residualsvector = TR.buildtomomat(gridpar, rays, residuals)
  
  # Perform the actual inversion using a "least-squares" approach
  postm,postC_m = TR.lininv(tomomat,cov_m,cov_d,mprior,residualsvector)


API, list of functions
-----------------------------
.. automodule:: ttimerays
   :members:
   :imported-members:  
.. Need to use :imported-members: since ttimerays imports
   the functions in __init__.py, so has no own members.

