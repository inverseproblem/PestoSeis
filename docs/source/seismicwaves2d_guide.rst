.. role:: raw-math(raw)
    :format: latex html
.. _seismicwaves2d_guide:


*******************************************************
Seismic wave propagation -- using ``seismicwaves2d``
*******************************************************

=============================================
Acoustic wave equation in 2D
=============================================

We model acoustic wave propagation through a medium by relating the time and space dependent pressure wavefield :math:`p(\mathbf{x},t)` to some external force :math:`f(\mathbf{x},t)` via theconstant density acoustic wave equation

.. math:: 
        -\frac{1}{c(\mathbf{x})^2}\frac{\partial}{\partial t^2}p(\mathbf{x},t)+\nabla^2 p(\mathbf{x},t)=f(\mathbf{x},t). :label: acoustic_wave_eq

Here, the properties of the medium are parametrized in terms of velocity :math:`c(\mathbf{x},t)`. In 2D, :math:`\mathbf{x}=[x,y]^{\text{T}}` and the spatial derivatives are given by :math:`\nabla=\partial_x(\cdot)+\partial_y(\cdot)`. The acoustic wave equation :eq:`acoustic_wave_eq` explicitly shows the relationship between the velocity structure of the medium :math:`[c(\mathbf{x},t)]` and the pressure wavefield :math:`p(\mathbf{x},t)`. Equation :eq:`acoustic_wave_eq` describes the forward problem of simulating the propagation of acoustic waves in a medium, hence solving this equation means that we are interested in obtaining the pressure wavefield given a specific velocity model.

In tomography, one is interested in unraveling the properties of the medium from some measurements of the wavefield, hence going in the reverse direction of the forward equation :eq:`acoustic_wave_eq`. Here we ask the question: given some measurements of the time dependent pressure wavefield, what is the velocity structure of the medium? One simple and frequently applied method that allows us to relate our observations to the velocity of the medium is straight-ray tomography. As the name suggests, the fundamental assumption in straight-ray tomography is that waves can be modelled as rays of infinite frequency and that they propagate along straight lines from a source to a receiver. Introducing the slowness as :math:`s=c^{-1}`, the observed travel time :math:`t_i`- hence the first arrival of the wave in the seismogram - between a source and a receiver can now be related to the structure of the medium by the line integral 
.. math::
	t_i=\int_S^Rs(l)dl,
	:label:linear_SRT}
where :math:`t_i` is the travel time measured at receiver :math:`R`, :math:`S` is the source position and :math:`s(l)` the slowness along the straight path from the source to the receiver. Note that by defining the slowness as the inverse of the velocity of the medium, equation :eq:`linear_SRT` is linear. 

In practice, one needs to solve the discrete version of equation :eq:`linear_SRT`. To this end, a computational model of the domain of interest is constructed by subdividing the domain in a rectilinear grid where the slowness :math:`s_j` within each cell :math:`j` is constant. Using the straight-ray approximation a ray then propagates as a straight line from an source :math:`S` to a receiver :math:`R` through the gridded domain, covering a certain distance :math:`r_j` within a cell :math:`j`. The time-of-flight of an entire ray :math:`i` between two points is given by the line connecting them and is calculated by the sum of the cell segments over the corresponding slowness within the cell as
.. math::
	t_i=\sum_jr_js_j.
	:label:discrete_SRT
Given a discretized velocity model of the medium, we can now trace out the straight-ray paths linking a source and a receiver by solving equation :eq:`discrete_SRT`. 

=============================================
Elastic wave equation in 2D
=============================================

Bla bla bla bla bla



=============================================
API ``seismicwaves2d``, list of functions
=============================================

.. automodule:: pestoseis.seismicwaves2d
   :members:
   :imported-members:  
..   :undoc-members:
