---
title: 'PestoSeis: A Python package for educational seismology'
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
date: 4 October 2022
bibliography: PestoSeis.bib

---

# Summary

`PestoSeis` is a Python package aimed at educational seismology. It contains
tools to solve two-dimensional seismic problems in terms of traveltimes, rays,
acoustic and elastic wave propagation, and related plotting functions. Moreover,
a set of functions performing basic seismic processing for exploration
geophysics is provided. The entire package is written in the Python language to
allow users to explore the code to gain a better understand of the numerical algorithms
used while simultaneously minimizing third-party dependencies and, thus,
simplifying the installation process. A set of illustrative examples covering
all included algorithms are provided.


# Statement of Need

`PestoSeis` is a numerical laboratory written in the Python language and 
contains a suite of functionalities which can be used as learning tools for 
concepts in seismology. One of the primary design goals of `PestoSeis` is simplicity and ease of use. It
contains a set of functions related to different aspects of seismology, ranging
from simple straight ray computations to a full-waveform elastic solver for
two-dimensional (2D) problems. It thus collects several routines within a single
package allowing the user to explore different topics in seismology using a
hands-on approach. Moreover, a complete set of examples demonstrating the
functionalities of the code are provided, with the combined aim of illustrating some
common seismological applications in addition to demonstrating how to utilize 
the software.

The entire package is written in Python to avoid the need of linking to external
libraries written in a compiled language and to allow the user to directly
explore the code in order to understand how the algorithms work. Furthermore,
dependencies on other third-party packages are kept to a minimum to simplify the installation
process while also enhancing the interpretability of the source code for the user.
`PestoSeis` can thus be used, for instance, in a course on seismology where a
set of exercises can be designed to help illustrate how certain principles can be applied in practice.

`PestoSeis` addresses the following set of 2D seismic problems:

- computation of traveltimes for given source and receiver positions provided a layered or grid-based velocity model and assuming straight rays;
- computation of traveltimes for given source and receiver positions provided a grid-based velocity model using a finite-difference method (fast marching);
- backtracing rays from the result of traveltime calculations in heterogeneous models;
- performing simple linearized ray tomography using a least squares approach;
- computing acoustic or elastic seismic wave propagation using a finite-difference method;
- performing some basic seismic processing for exploration seismology;
- plotting various input and output data related to the aforementioned applications.

<!-- This package enables users to enhances their understanding of a variety of
modeling principles used within seismology using a simple yet comprehensive set
of tools. -->

# Package Content

`PestoSeis` contains three main submodules:

1. traveltimes and rays (`ttimerays`);
2. seismic wave propagation (`seismicwaves2D`);
3. seismic processing (`reflectionseismo`).

Each one of these categories considers a specific use case and provides a
variety of functions tailored at solving problems within the considered
application. In the following, an overview of the core functionalities of the
submodules are provided.

## Traveltimes and rays

`PestoSeis` provides functions to perform the following computations in heterogeneous media:

1. traveltimes given a velocity model by solving the eikonal equation [@sethianFastMarchingLevel1996; @rawlinsonWaveFrontEvolution2004]; 
2. trace "bent" rays;
3. straight rays;
4. in special case of a horizontally layered medium, use Snell's law to compute ray paths, traveltimes, and the distance covered by the ray.
 
The functions provided in `pestoseis.ttimerays` allow the user to compute rays
and traveltimes for a 2D velocity model which is discretized on a rectilinear
grid constructed from inputs provided by the user.  The available functions allow the user to
set up and solve simple 2D tomographic inverse problems and thus aims to
serve as a guide for how travel time tomography can be performed in practice.

Figure \ref{ttimes_and_rays} shows an example from the field of medical imaging
with computed traveltimes, bent ray paths, and straight ray paths through a
speed-of-sound model mimicking human breast tissue (e.g. @medicalUS). The computed traveltimes can
then be used to set up a tomographic inverse problem which may, for instance, be
solved with a simple linear inversion under Gaussian assumptions (least squares
approach).

![Visualization of computed travel times and rays (bent rays and straight rays) using the functions in `pestoseis.ttimerays`. In this example, a numerical breast phantom is considered as used in medical ultrasound to study the ray paths through the medium and subsequently compute the according travel times from an array of sources and receivers surrounding the breast phantom. \label{ttimes_and_rays}](figs/tutorial04_results.png)


## Seismic wave propagation

`PestoSeis` provides means to compute wave propagation in 2D for:

1. acoustic media solving the acoustic wave equation;
2. elastic media solving the elastic wave equation.

![A snapshot of an acoustic seismic wavefield where the wavefield is overlayed on the velocity model. The depicted overthrust model is adapted from @aminzadeh_3-d_1997.\label{acouwavefield}](figs/acouwavefield1.png)

The acoustic and elastic wave equations in 2D are solved using finite
differences on a staggered grid in both space and time in the
`pestoseis.seismicwaves2D` submodule [@virieuxPSVWavePropagation1986;
@komatitschUnsplitConvolutionalPerfectly2007]. Staggered grids generally allow
for better accuracy with a minimal increase in computational overhead. 

In order to compute wave propagation in either an acoustic or an elastic medium, the user
needs to specify the parameters of the grid used to construct the velocity
models, the source-receiver geometry, and a source time function with a given dominant
frequency. Furthermore, the user must select the desired boundary conditions as either 
reflecting, a Gauss taper (acoustic formulation only), or convolutionary perfectly matched layers [@komatitschUnsplitConvolutionalPerfectly2007];
free surface boundary conditions may also be optionally set at the top of the model. 

The function returns a seismogram recorded
at the specified receiver locations as well as a set of snapshots of the wavefield (see figure
\ref{acouwavefield}). Additionally, it is possible to save a collection of 
wavefield snapshots as a `.mp4` movie file that illustrates the propagation of the seismic
waves through the medium. These functions aim to equip the user with the
ability to quickly set up and visualize small scale 2D simulations.

## Seismic processing
 
![An example of seismic processing where the original data are corrected for geometrical spreading and amplitude gain. \label{geomspragc}](figs/geomspreagc.png)

The submodule `pestoseis.reflectionseismo` provides a series of routines for
processing the resulting data using a number of methods which are commonly used
in practice within seismology. Examples of processing routines implemented
within `PestoSeis` include arranging the data in shot gathers, generating a
wiggle plot of the shot gathers, normal moveout (NMO) correction, correcting for
geometrical spreading, and applying automatic gain control (AGC) to a shot
gather (see figure \ref{geomspragc}). Furthermore, some functionalities which
can be used for filtering the data in the frequency-wavenumber domain are also
provided. Similar to the forward modeling tools outlined previously, the main
focus of these processing functions is to provide the user with rudimentary
tools which are commonly used in seismology so that they can experiment with
processing methods which are commonly used in both industry and academia. 

## Tutorials

Finally, a number of tutorials that provide examples on how to use the functions
within `PestoSeis` are provided in the form of Jupyter notebooks. These
tutorials showcase different numerical scenarios and can be used to get started
with `PestoSeis`.

# References

