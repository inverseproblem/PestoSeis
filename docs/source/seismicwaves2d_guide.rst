.. role:: raw-math(raw)
    :format: latex html
.. _seismicwaves2d_guide:


*******************************************************
Seismic wave propagation -- using ``seismicwaves2d``
*******************************************************

=============================================
Acoustic wave equation in 2D
=============================================

We model acoustic wave propagation through a medium by relating the time and space dependent pressure wavefield $p(\mathbf{x},t)$ to some external force $f(\mathbf{x},t)$ via the scalar acoustic wave equation for loss-less media
\begin{equation}
-\frac{1}{c(\mathbf{x})^2}\frac{\partial}{\partial t^2}p(\mathbf{x},t)+\rho(\mathbf{x})\nabla\frac{1}{\rho(\mathbf{x})}\nabla p(\mathbf{x},t)=f(\mathbf{x},t).
\label{eq:acoustic_wave_eq}
\end{equation}
Here, the properties of the medium are parametrized in terms of velocity $c(\mathbf{x},t)$ and density $\rho(\mathbf{x})$. In 2D, $\mathbf{x}=[x,y]^{\text{T}}$ and the spatial derivatives are given by $\nabla=\partial_x(\cdot)+\partial_y(\cdot)$. The acoustic wave equation \ref{eq:acoustic_wave_eq} explicitly shows the relationship between parameters of the medium $[c(\mathbf{x},t);\rho(\mathbf{x})]$ and the pressure wavefield $p(\mathbf{x},t)$. Equation \eqref{eq:acoustic_wave_eq} describes the forward problem of simulating the propagation of acoustic waves in a medium, hence solving this equation means that we are interested in obtaining the pressure wavefield given a specific velocity and density model.

In tomography, one is interested in unraveling the properties of the medium from some measurements of the wavefield, hence going in the reverse direction of the forward equation \eqref{eq:acoustic_wave_eq}. Here we ask the question: given some measurements of the time dependent pressure wavefield, what is the velocity structure of the medium? One simple and frequently applied method that allows us to relate our observations to the velocity of the medium is straight-ray tomography. As the name suggests, the fundamental assumption in straight-ray tomography is that waves can be modelled as rays of infinite frequency and that they propagate along straight lines from a source to a receiver. Introducing the slowness as $s=c^{-1}$, the observed travel time $t_i$- hence the first arrival of the wave in the seismogram - between a source and a receiver can now be related to the structure of the medium by the line integral 
\begin{equation}
	t_i=\int_S^Rs(l)dl,
	\label{eq:linear_SRT}
\end{equation}
where $t_i$ is the travel time measured at receiver $R$, $S$ is the source position and $s(l)$ the slowness along the straight path from the source to the receiver. Note that by defining the slowness as the inverse of the velocity of the medium, equation \eqref{eq:linear_SRT} is linear. 

In practice, one needs to solve the discrete version of equation \eqref{eq:linear_SRT}. To this end, a computational model of the domain of interest is constructed by subdividing the domain in a rectilinear grid where the slowness $s_j$ within each cell $j$ is constant. Using the straight-ray approximation a ray then propagates as a straight line from an source $S$ to a receiver $R$ through the gridded domain, covering a certain distance $r_j$ within a cell $j$. The time-of-flight of an entire ray $i$ between two points is given by the line connecting them and is calculated by the sum of the cell segments over the corresponding slowness within the cell as
\begin{equation}
	t_i=\sum_jr_js_j.
	\label{eq:discrete_SRT}
\end{equation}
Given a discretized velocity model of the medium, we can now trace out the straight-ray paths linking a source and a receiver by solving equation \eqref{eq:discrete_SRT}. 

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
