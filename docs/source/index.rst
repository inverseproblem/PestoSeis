.. PestSeis documentation master file, created by
   sphinx-quickstart on Wed May  1 10:06:22 2019.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.


.. 
    # with overline, for parts
    * with overline, for chapters
    =, for sections
    -, for subsections
    ^, for subsubsections
    ", for paragraphs



##############################
PestoSeis's documentation
##############################

The PestoSeis package is a set of functions to perform calculations of seismic traveltimes, rays and wave propagation in two dimensions (2D). PestoSeis covers the range from forward modeling wave propagation through a known model to reconstructing the structural map of an unknown medium from some observed (synthetic) data by providing algorithms to solve a 2D inverse problem by simple means.

Its main purpose is educational, it aims to be a simple tool to play around with some aspects of seismology. Due to the fact that many fields of research use the fundamentals of wave propagation, the functionalities within PestoSeis are not limited to seismology, but can also be applied for instance to problems in medical ultrasound (check out the tutorial on medical ultrasound). Due to the educational purpose of the code, no optimization for speed has been perfomed and it has *not* been tested extensively. Therefore, use of this code for research it is *not* encouraged.

The functions contained in this package aim at the following:

* Setting up a grid containing a velocity model;
* Computing traveltimes for given sources and receivers, optionally providing the entire 2D traveltime model;
* Tracing seismic rays given a 2D traveltime array and position of sources and receivers;
* Performing very simple traveltime tomography using rays;
* Computing acoustic and elastic wave propagation in 2D by finite differences;
* Performing basic seismic processing for exploration seismology. 


:Authors: 
   - Andrea Zunino 
   - Patrick Marty 
   - Ines E. Ulrich

==============
Installation
==============
     
See the :ref:`installation` for more details.

==================
User guide
==================
The User guide provides information about functionality and how to use the package with examples.

-----------------------------
Traveltimes and rays in 2D
-----------------------------

See :ref:`ttimerays_guide`.

--------------------------
Wave propagation in 2D
--------------------------

See :ref:`seismicwaves2d_guide`.

  
==================
Tutorials
==================

A set of examples in the form of Jupyter notebooks illustrating the capabilities of the code.
See :ref:`tutorials`.


==================
For developers
==================

See :ref:`developerdocs`.


#####################
Table of contents
#####################

.. toctree::
      :maxdepth: 3

      installation
      ttimerays_guide
      seismicwaves2d_guide
      tutorials
      API/index      
      developerdocs

======================
Indices and tables
======================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
