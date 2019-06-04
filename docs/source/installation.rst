
.. _installation:


************************************
Installation --- loading the module
************************************

.. warning:: Installing the package globally with, for instance, ``pip`` is still untested.
	     
A simple way to load the package/module from a local folder is to start a script with adding to the path the location the folder containing the package. This way Python will be able to find the module and therefore load it.

If we currently are in the directory ``dir1`` and the package is located in the directory ``pyteachseismo`` ::

  dir1/
     myscript.py
     
     pyteachseismo/
         pestoseis/
	     ttimerays/
	     ttimerays.py
	     ...
	     
	     seismicwaves2d/
	     elasticwaveprop2D.py
	     ...
	     
	
Then something along the along the following lines in ``myscript.py`` should work: ::
 
  import sys
  # add location of ttimerays package to path
  sys.path.append('pyteachseismo/')
  import pestoseis as PS

At this point, the functions from ``pestoseis`` should be available as, for instance::

  PS.traceallrays(gridpar, srccoo, reccoo, grdsttime)

Altenatively, if, for instance, interested in only the traveltime/rays routines: ::
  
  import pestoseis.ttimerays as TR

Otherwise, only the wave propagation routines: ::
  
  import pestoseis.seismicwaves2d as SW

