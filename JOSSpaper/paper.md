---
title: 'PestoSeis: A Python package for basic and educational seismology'
tags:
  - Python
  - seismology
  - geophysics
  - traveltimes
  - rays
  - acoustic waves
  - elastic waves
  - seismic processing
authors:
  - name: Andrea Zunino^[corresponding author]
    orcid: 0000-0002-3415-162X
    affiliation: "1"
  - name: Patrick Marty
    affiliation: 1
  - name: Ines Elisa Ulrich
    affiliation: 1
  - name: Andreas Fichtner
    affiliation: 1 
affiliations:
 - name: Institute of Geophysics, ETH Zurich, Switzerland
   index: 1
date: 14 November 2022
bibliography: PestoSeis.bib

---

# Summary

`PestoSeis` is a Python package which contains a collection of solvers and 
processing tools commonly used in acoustic wave physics.  With a particular
emphasis on seismological applications, `PestoSeis` contains
tools to solve two-dimensional seismic problems in terms of traveltimes, rays,
acoustic and elastic wave propagation, and related plotting functions. Moreover,
a set of functions performing basic seismic processing for exploration
geophysics is provided. The entire package is written in the Python language to
allow users to explore the code and to gain a better understanding of the numerical algorithms
used. Third-party dependencies are kept to a minimum, thus
simplifying the installation process. A set of illustrative examples using Jupyter notebooks
are provided to demonstrate how the included algorithms work.


# Statement of Need

`PestoSeis` is a numerical laboratory written in the Python language and 
contains a suite of functionalities which can be used to perform commonly used calculations
in acoustics, with a focus on seismic applications. 
More specifically, `PestoSeis` contains functions ranging from simple straight ray 
computations to a full-waveform elastic solver for solving the wave equation in
two-dimensional (2D) problems. 

One of the primary design goals of `PestoSeis` is simplicity and ease of use. It thus collects several routines within a single
package, allowing the user to explore different topics in seismology and acoustic wave physics using a
hands-on approach. Moreover, a complete set of examples demonstrating the
functionalities of the code are provided, with the combined aim of illustrating some
common applications in seismology in addition to other potential applications such as in medical ultrasound.

The entire package is written in Python to avoid the need of linking to external
libraries written in a compiled language as well as to allow the user to directly
explore the code in order to understand how the algorithms work. Furthermore,
dependencies on other third-party packages are kept to a minimum to simplify the installation
process while also enhancing the interpretability of the source code for the user.
Basic concepts of seismology, optics, various problems modeled by the eikonal equation or 
acoustic/elastic wave equation in 2D can be numerically explored with this package.
`PestoSeis` can also be used, for instance, in a course on seismology where a
set of exercises can be designed to help illustrate how certain principles can be applied in practice.

`PestoSeis` addresses the following set of 2D seismic problems:

- computation of traveltimes for given source and receiver positions provided a layered or grid-based velocity model and assuming straight rays;
- computation of traveltimes for given source and receiver positions provided a grid-based velocity model using a finite-difference method (fast marching);
- backtracing rays from the result of traveltime calculations in heterogeneous models;
- performing simple linearized ray tomography using a least squares approach;
- computing acoustic or elastic seismic wave propagation using a finite-difference method;
- performing some basic seismic processing for exploration seismology;
- plotting various input and output data related to the aforementioned applications.


# Package Content

`PestoSeis` contains three main submodules:

1. traveltimes and rays (`ttimerays`);
2. seismic wave propagation (`seismicwaves2D`);
3. seismic processing (`reflectionseismo`).

Each one of these categories considers a specific use case and provides a
variety of functions tailored at solving problems within the considered
application. In the following, an overview of the core functionalities of the
submodules are provided.

## Traveltimes and Rays in 2D

`PestoSeis` provides functions to perform the following computations in heterogeneous media:

1. traveltimes given a velocity model by solving the eikonal equation [@sethianFastMarchingLevel1996; @rawlinsonWaveFrontEvolution2004]; 
2. trace "bent" rays by following the negative gradient of traveltimes from receiver to source;
3. straight rays;
4. in the special case of a horizontally layered medium, use Snell's law to compute ray paths, traveltimes, and the distances covered by the ray.
 
The functions provided in `pestoseis.ttimerays` allow the user to compute rays
and traveltimes for a 2D velocity model which is discretized on a rectilinear
grid constructed from inputs provided by the user.  The available functions allow the user to
set up and solve simple 2D tomographic inverse problems and, thus, aims to
serve as a guide for how travel time tomography can be performed in practice. 
Ray theory is used in several fields such as seismic tomography, optics, and atmospheric and ocean acoustics. One interesting example is that of medical imaging, where traveltimes (time-of-flight in the medical literature) are used to infer the internal structure of parts of the human body [@ulrichDiffuseUltrasoundComputed2022].

Figure \ref{ttimes_and_rays} shows an illustrative example from the field of medical imaging
with computed traveltimes, bent ray paths, and straight ray paths through a
speed-of-sound model mimicking human breast tissue (e.g. @ulrichDiffuseUltrasoundComputed2022). 
The computed traveltimes can
then be used to set up a tomographic inverse problem which may, for instance, be
solved with a simple linear inversion under Gaussian assumptions (least squares
approach) [@tarantolaInverseProblemTheory2005c], as provided by `PestoSeis`.

![Visualization of computed travel times and rays (bent rays and straight rays) using the functions in `pestoseis.ttimerays`. A numerical breast phantom is employed, as often used in medical ultrasound to study the ray paths through the medium. An array of receivers surrounds the breast phantom. The shown ray paths originate from two sources at the left and right sides of the phantom. \label{ttimes_and_rays}](figs/tutorial04_results.png)


## Wave Propagation in 2D

`PestoSeis` provides means to compute wave propagation in 2D for:

1. acoustic media solving the acoustic wave equation;
2. elastic media solving the elastic wave equation.

![A snapshot of an acoustic seismic wavefield where the wavefield is superimposed on the velocity model. The location of the source is indicated by the black dot. The wavefield represents the pressure distribution at a single instance in time. The depicted overthrust model is adapted from @aminzadehSEGEAGE3D1997.\label{acouwavefield}](figs/wavefield.png)

The acoustic and elastic wave equations in 2D are solved using finite
differences on a staggered grid in both space and time in the
`pestoseis.seismicwaves2D` submodule [@bunksMultiscaleSeismicWaveform1995;@virieuxPSVWavePropagation1986]. Staggered grids generally allow
for better accuracy with a minimal increase in computational overhead. Absorbing boundary conditions are 
implemented as convolutional perfectly matched layer (C-PML) for both acoustic and elastic formulations [@pasalicConvolutionalPerfectlyMatched2010;@komatitschUnsplitConvolutionalPerfectly2007].

In order to compute wave propagation in either an acoustic or an elastic medium, the user
needs to specify the parameters of the grid used to construct the velocity
models, the source-receiver geometry, and a source time function with a given dominant
frequency. Furthermore, the user must select the desired boundary conditions as either 
reflecting, a Gauss taper (acoustic formulation only), or convolutional perfectly matched layers [@komatitschUnsplitConvolutionalPerfectly2007];
free surface boundary conditions may also be optionally set at the top of the model. 

The function returns a seismogram recorded
at the specified receiver locations as well as a set of snapshots of the wavefield (see figure
\ref{acouwavefield}). Additionally, it is possible to save a collection of 
wavefield snapshots as an `.mp4` movie file that illustrates the propagation of the seismic
waves through the medium. These functions aim to equip the user with the
ability to quickly set up and visualize small-scale 2D simulations.

## Seismic processing
 
![An example of seismic processing where reflection acoustic data from a shotgather are corrected for geometrical spreading (decay in energy with distance) and enhanced with amplitude gain (decay of amplitude with time). The plots on the top and bottom rows depict the same data but with two different plotting techniques.\label{geomspragc}](figs/geomspreagc.png)

The submodule `pestoseis.reflectionseismo` provides a series of routines for
processing the resulting data using a number of methods which are commonly used
in practice within seismology [@ozSeismicDataAnalysis2001]. Examples of processing routines implemented within `PestoSeis` include arranging the data in shot gathers, generating a
wiggle plot of the shotgathers, normal moveout (NMO) correction, correcting for
geometrical spreading, and applying automatic gain control (AGC) to a shot
gather (see figure \ref{geomspragc}). Furthermore, some functionalities which
can be used for filtering the data in the frequency-wavenumber domain are also
provided. Similar to the forward modeling tools outlined previously, the main
focus of these processing functions is to provide the user with basic 
tools which are commonly used in seismology. The users can thus experiment with
processing methods which are commonly used in both industry and academia. 

## Tutorials

Finally, a number of tutorials that provide examples on how to use the functions
within `PestoSeis` are provided in the form of Jupyter notebooks. These
tutorials showcase different numerical scenarios and can be used to get started
with `PestoSeis`.

# References

