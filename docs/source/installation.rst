
.. _installation:


************************************
Installation --- loading the module
************************************

The simplest way is to use the `pip` installer:
::
   pip install git+https://gitlab.com/anzun/PestoSeis

Otherwise, if a local copy is available, a simple way to load the package/module from a local folder is to start a script with adding to the path the location the folder containing the package. This way Python will be able to find the module and therefore load it.
Then something along the along the following lines in ``myscript.py`` should work: ::
 
  import sys
  # add location of PestoSeis package to path
  sys.path.append('PestoSeis/')
  import pestoseis as PS

At this point, the functions from ``pestoseis`` should be available as, for instance::

  PS.traceallrays(gridpar, srccoo, reccoo, grdsttime)

Altenatively, if, for instance, interested in only the traveltime/rays routines: ::
  
  import pestoseis.ttimerays as TR

Otherwise, only the wave propagation routines: ::
  
  import pestoseis.seismicwaves2d as SW

