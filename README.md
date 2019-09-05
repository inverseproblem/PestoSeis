
# PestoSeis #

[![Docs](https://img.shields.io/badge/docs-blue.svg)](https://inverseproblem.github.io/PestoSeis/)

The PestoSeis package is a set of functions to perform calculations of seismic traveltimes, rays and wave propagation in two dimensions (2D).

Its main purpose is educational, it aims to be a simple tool to play around with some aspects of seismology. Because of that, no optimization of the code for speed has been perfomed and it has *not* been tested extensively. Therefore, do not use this code for research.

The functions contained in this package aim at the following:

* Setting up a grid containing a velocity model;
* Computing traveltimes for given sources and receivers, optionally providing the entire 2D traveltime model;
* Tracing seismic rays given a 2D traveltime array and position of sources and receivers;
* Performing very simple traveltime tomography using rays;
* Computing acoustic and elastic wave propagation in 2D by finite differences. 

PestoSeis' docs: <https://github.com/inverseproblem/PestoSeis.git>


## Authors: ##
   Andrea Zunino



