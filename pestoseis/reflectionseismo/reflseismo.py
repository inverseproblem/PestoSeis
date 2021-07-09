
#------------------------------------------------------------------------
#
#    PestoSeis, a numerical laboratory to learn about seismology, written
#    in the Python language.
#    Copyright (C) 2021  Andrea Zunino 
#
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
#------------------------------------------------------------------------

"""
Functions to generate and process seismic reflection data to mimic an 
  exploration seismology setup.
"""

#######################################################################
#######################################################################

import numpy as __np
import matplotlib.pyplot as __plt
from scipy.interpolate import CubicSpline as __CubicSpline
from skimage.draw import polygon
from scipy.ndimage import gaussian_filter

#######################################################################

def fwdconvolve(refle,wavelet,dt):
    """
    Convolve a reflectivity series with a wavelet to compute a seismogram.

    Args:
       refle (ndarray): reflectivity series
       wavelet (ndarray): the wavelet
       dt (float): time interval
    
    Returns:
       tarr (ndarry): time array
       convo (ndarry): seismic trace resulting from the convolution

    """
    convo = __np.convolve(refle,wavelet,mode="full")
    tarr = __np.array([i*dt for i in range(convo.size)])

    return tarr,convo

#######################################################################

def calcreflectivity(density,vel,z,dt):
    """
    Compute the reflectivity series given density and velocity as a function 
      of depth and a time interval.

    Args:
       density (ndarray): a vector of density values
       vel (ndarray): a vector of velocity values
       z (ndarray): a vector of depths
       dt (float): the time interval

    Returns:
       twt (ndarray): the two-way travel time
       refltwt (ndarray): the reflectivity series as a function of two-way travel time
    """
    assert(vel.ndim==1)
    assert(density.ndim==1)

    impdz = density*vel
    reflz = (impdz[1:]-impdz[:-1])/(impdz[1:]+impdz[:-1])

    twtfromv = _depth2time(z,vel)
    twt = __np.arange(0.0,twtfromv.max(),dt)
    npts = twt.size    

    refltwt = __np.zeros(npts)
    for i in range(reflz.size):
        if __np.abs(reflz[i])>0.0 :
            idx = __np.argmin(__np.abs(twtfromv[i]-twt))
            refltwt[idx] = reflz[i]
        
    return twt,refltwt

#######################################################################

def _depth2time(z,vel):
    """
    Convert depth to time. 
    
    Args: 
       z (ndarray): the depth array
       vel (ndarray): the velocity array

    Returns:
       twt (ndarray): the two-way travel time

    """
    assert(vel.ndim==1)
    assert(vel.size==z.size)

    npts=vel.size
    twt = __np.zeros(npts)

    for i in range(npts-1):
        deltaz=z[i+1]-z[i]
        assert(deltaz>=0.0)
        twt[i+1] = twt[i] + 2.0*deltaz/vel[i]
        
    return twt

#######################################################################

def imgshotgath(seisdata,dt,offset,amplitudeclip=1.0):
    """
    Create an image of a shotgather.

    Args:
       seisdata (ndarray): seismic data, i.e., shotgather, *rows* contain traces
       dt (float): time interval
       offset (ndarray): array of offsets 
       amplitudeclip (float,optional): clip the amplitude to a desired value

    """
    twt = __np.array([i*dt for i in range(seisdata.shape[1])])  
    vmax = amplitudeclip*__np.abs(seisdata).max()
    vmin = -vmax
    extent = [offset.min(),offset.max(),twt.max(),twt.min()]
    __plt.imshow(seisdata.T,vmin=vmin,vmax=vmax,cmap=__plt.cm.gray,
               extent=extent,aspect='auto',interpolation='bilinear')
    __plt.colorbar()
    __plt.xlabel('Offset [m]')
    __plt.ylabel('TWT [s]')
    return

#######################################################################

def wiggle(data,dt,offset=None,skiptr=1,scal=None,title=None,filltrace=True):

    """
    Create a 'wiggle' plot of a shotgather.

    Args:
       data (ndarray): seismic data, i.e., shotgather, *rows* contain traces
       dt (float): time interval
       offset (ndarray,optional): array of offsets 
       skiptr (integer,optiona): plot only every 'skiptr' traces
       title (string,optional): add a title to the plot

    """
    t = __np.array([i*dt for i in range(data.shape[1])])        
    lwidth = 1.0 # line width
    ntraces = data.shape[0]

    if not isinstance(offset,__np.ndarray) :
        ofs = __np.array([float(i) for i in range(ntraces)])
        maxval=1.0/((abs(data).max()))
    else :
        ofs = offset
        maxseis = abs(data).max()
        maxoff = (abs(__np.diff(offset).max()))
        maxval = maxoff/maxseis

    if scal!=None :
        maxval *= scal    

    if data.ndim==1:
        data=data.reshape(1,-1)

    for i in range(0,ntraces,skiptr):
        trace=ofs[i]+data[i,:]*maxval
        __plt.plot(trace,t,color='black',linewidth=lwidth)
        if filltrace :
            __plt.fill_betweenx(t,ofs[i],trace,where=(trace>=ofs[i]),color='black')
        #__plt.fill_betweenx(t,ofs[i],trace,where=(trace<ofs[i]),color='red')
    ax = __plt.gca()

    __plt.xlim(ofs[0]-data[0,:].max()*maxval,ofs[-1]+data[-1,:].max())
    __plt.ylim(t.min(),t.max())
    ax.invert_yaxis()

    if title!=None :
        __plt.title(title)
    __plt.xlabel('Offset [m]')
    __plt.ylabel('TWT [s]')
    return

#######################################################################

def geometrical_spreading(seis,twt):
    """
    Apply geometrical spreading correction to a shotgather.

    Args:
        seis (ndarray): seismic shotgather, *rows* contain traces
        twt (ndarray): two-way traveltime 
            
    Returns
        seis_gs (ndarray): seismograms corrected for geometrical spreading
    """
    ns = seis.shape[1]
    ntraces = seis.shape[0]
    
    seis_gs = seis.copy()
        
    for i in range(ntraces):
        ## twt**2 is a rough APPROXIMATION to geom. spreading...
        seis_gs[i,:] = seis_gs[i,:] * twt**2
        ## scale max amplitude to 1
        seis_gs[i,:] /= __np.abs(seis_gs[i,:]).max()
        
    return seis_gs

#######################################################################

def agc(seis, w=100, rho=0, type='uniform'):
    """
    Apply Automatic Gain Control to a shotgather.

    Args:
        seis (ndarray): seismic shotgather, *rows* contain traces
        w (float): window width
        rho (ndarray,optional): density
        agctype (string): 'linear' or 'gaussian', weight kernel type

    Returns
        seis_gs (ndarray): seismograms corrected with AGC
    """

    ns = seis.shape[1]
    ntraces = seis.shape[0];
    
    if rho>0:
        nw = int(__np.ceil(w/rho))
    else:
        nw = w
        
    # select weight kernel type
    if type.lower() == 'gaussian':
        from scipy.stats import norm
        iw=__np.linspace(-2,2,2*nw);
        weight_kernel = norm.pdf(iw)
        weight_kernel = weight_kernel/__np.sum(weight_kernel)
    else:
        weight_kernel = __np.ones(nw)
        weight_kernel = weight_kernel/__np.sum(weight_kernel)
 
    seis_agc=__np.zeros_like(seis)
    for i in range(ntraces):
        
        seis_abs = __np.abs(seis[i,:])
        seis_weight = __np.convolve(__np.abs(seis[i,:]),weight_kernel, mode='same')

        if (seis_weight==0.0).any() :
            # print("agc(): Dividing by zero in AGC correction... ")
            seis_weight = __np.where(seis_weight==0.0,1e-6,seis_weight)      
            
        ##-----
        seis_agc[i,:] = seis[i,:]/seis_weight
  
    return seis_agc
   
#######################################################################

def nmocorrection(velnmo,dt,offset,seisdat):
    """
    Common Mid Point (CMP) normal moveout correction.
  
    Args:
        velnmo (ndarray): 1D array of "NMO" velocity, representing the 
                          average velocity of the layers above
        dt (float): the time interval
        offset (ndarray): offsets for each trace
        seisdat (ndarray): input seismic data, *rows* contain traces

    Returns
        seisnmo (ndarray): seismograms NMO corrected
    """

    ntraces,nsmpl = seisdat.shape
    assert ntraces==offset.size
    seisnmo = __np.zeros_like(seisdat)
    timearr = __np.array([(i-1)*dt for i in range(nsmpl)])

    ## build the nmo array trace by trace
    for itr in range(ntraces):
        # compute the time for nmo for all time points of a single trace
        # velocity array is nx x nz so slice in depth is vel[?,:]
        tnmo = __np.sqrt(timearr**2 + offset[itr]**2/velnmo**2) 
        # new trace
        seisnmo[itr,:] =  _resampletrace(timearr,tnmo,seisdat[itr,:])
    return seisnmo

#######################################################################

def _resampletrace(torig,tnmo,seistr):
    """
    Resample a trace
    
    Args:
        torig
        tnmo
        seistr
        
    Returns
        seisnew (ndarray): resampled trace
    """
    seisnew = __np.zeros(torig.size)
    itp = __CubicSpline(torig,seistr)
    # time outside bounds?
    outbou = (tnmo>torig.max()).nonzero()[0]
    if outbou.size>0:
        idx = __np.min(outbou)
    else :
        idx = seisnew.size
    seisnew[:idx] = itp(tnmo[:idx])
    return seisnew

#######################################################################

def get_pm_mask(seisdat):
    """
    Gets a mask of the sign (+/-) for a given set of seismic traces. Useful
    for recoving the correct sign of the seismic data after FFTing the data.

    Args:
        seisdat (ndarray): input seismic data in time-space domain, *rows*
                           contain traces
    
    Returns
        pm_mask (ndarray): mask containing the sign (+1 or -1) of the trace
                           samples
    """
    pm_mask = __np.zeros_like(seisdat)
    pm_mask[__np.where(seisdat > 0)] = 1
    pm_mask[__np.where(seisdat < 0)] = -1

    return pm_mask

#######################################################################

def transform_tx_to_fk(seisdat):
    """
    Transform seismograms from time-space domain to frequency-wavenumber
    domain.

    Args:
        seisdat (ndarray): input seismic data, *rows* contain traces

    Returns
        fk_seis (ndarray): transformed seismic data
    """
    # Transform to frequency-space domain
    fx_seis = __np.fft.fftshift(__np.fft.fft(seisdat, axis=1), axes=1)

    # Transform to frequency-wavenumber domain
    fk_seis = __np.fft.fftshift(__np.fft.ifft(fx_seis, axis=0), axes=0)

    return fk_seis

#######################################################################

def transform_fk_to_tx(seisdat):
    """
    Transform seismograms from frequency-wavenumber domain to time-space
    domain.

    Args:
        seisdat (ndarray): seismic traces in frequency-wavenumber domain,
                           *rows* contain traces
    
    Returns
        tx_seis (ndarray): transformed seismic data
    """
    # Transform to frequency-space domain
    fx_seis = __np.fft.fft(__np.fft.ifftshift(seisdat, axes=0), axis=0)

    # Transform to time-space domain
    tx_seis = __np.fft.ifft(__np.fft.ifftshift(fx_seis, axes=1), axis=1)

    return __np.abs(tx_seis)

#######################################################################

import numpy as np
class _Canvas:
    def __init__(self, ax):
        """
        Constructs a canvas object for drawing on a plot

        Args:
            ax (matplotlib.axes._subplots.AxesSubplot): plot axis 
        """
        self.ax = ax
        
        # Create handle for a path of connected points
        self.path, = ax.plot([], [], 'o-', lw=2.5, color='red')
     
        self.vert = np.empty(shape=[0,2])
        self.ax.set_title('Interactive Plot\nLEFT: new point, MIDDLE: close polygon, RIGHT: delete last point')

        self.x = [] 
        self.y = []

        self.mouse_button = {1: self._add_point, 2: self._close_polygon, 3: self._delete_point}
    
    def set_location(self,event):
        if event.inaxes:
            self.x = event.xdata
            self.y = event.ydata
    
    def _add_point(self):
        self.vert = np.vstack((self.vert, np.array([self.x, self.y])))
    
    def _delete_point(self):
        if (self.vert).shape[0] > 0:
            self.vert = np.delete(self.vert, -1, 0)
    
    def _close_polygon(self):
        self.vert = np.vstack((self.vert, self.vert[0,:]))
    
    def update_path(self,event):

        # If the mouse pointer is not on the canvas, ignore buttons
        if not event.inaxes: return

        # Do whichever action correspond to the mouse button clicked
        # mouse_button: dictionary {1: self._add_point, 2: self._delete_point, 3: self._close_polygon}
        self.mouse_button[event.button]()
        print("vert:",self.vert)
        
        x = [self.vert[k,0] for k in range((self.vert).shape[0])]
        y = [self.vert[k,1] for k in range((self.vert).shape[0])]
        self.path.set_data(x,y)
        __plt.draw()

#######################################################################

def _tracepoly(fig) :
    ax = fig.gca()
    cnv = _Canvas(ax)

    # Get the information on the mouse events
    __plt.connect('button_press_event', cnv.update_path)
    __plt.connect('motion_notify_event', cnv.set_location)

    # Display the plot
    __plt.tight_layout()
    __plt.show()
        
    return cnv

#######################################################################

def plot_fk_spectrum(seisdat_fk, dt, dx, ylim=None, interactive=True):
    """
    Plots the frequency-wavenumber spectrum of the seismic traces.  Data
    must already be in frequency-wavenumber domain.

    Args:
        seisdat_fk (ndarray): seismic traces in frequency-wavenumber domain,
                              *rows* contain traces
        ylim (float, optional): upper y-limit (frequency) to plot
        interactive (bool, optional): whether to plot an interactive polygon
                                      drawing tool

    Returns
        canvas (_Canvas): canvas object containing the user-drawn polygon, only
                          returned when `interactive=True`.
    """
    # Verify that the data is complex; rudamentary check to see if the data
    # is still in time-space domain (but obviously not fool-proof)
    if not __np.any(__np.iscomplex(seisdat_fk)):
        raise Exception("Input data is not complex; verify that the input\
            data is in frequency-wavenumber domain")

    # Construct plot
    fig = __plt.figure(figsize=[10, 8])
    ax = __plt.gca()
    ax.imshow(
        __np.abs(seisdat_fk).T, 
        aspect="auto", 
        cmap="gray_r",
        origin="lower", 
    )

    freqs = np.fft.fftshift(np.fft.fftfreq(seisdat_fk.shape[1], d=dt))
    __plt.yticks(np.arange(0, seisdat_fk.shape[1], 25), np.round(freqs[::25], 0))

    ks = np.fft.fftshift(np.fft.fftfreq(seisdat_fk.shape[0], d=dx))
    __plt.xticks(np.arange(0, seisdat_fk.shape[0], 25), np.round(ks[::25], 3))


    if ylim is None:
        ylim = 0.1 * seisdat_fk.shape[1]
    __plt.ylim([seisdat_fk.shape[1]/2, seisdat_fk.shape[1]/2 + ylim])

    ax.set_title("Frequency-Wavenumber Spectrum")
    ax.set_xlabel("Wavenumber [1/m]")
    ax.set_ylabel("Frequency [Hz]")

    if interactive:
        # Set canvas for drawing polygons onto
        canvas = _tracepoly(fig)
        return canvas

#######################################################################

def mask_from_polygon(img, poly, smooth=True, smoothing_strength=2.5, invert=False):
    """
    Constructs a boolean mask from a polygon for a given image.  Regions within
    the polygon are by default assigned `0` and values outside of the polygon
    are assigned `1`.

    Args:
        img (ndarray): image to construct the mask for
        poly (ndarray, _Canvas): polygon to construct the boolean
                                    mask of
        smooth (bool, optional): flag for whether to apply gaussian smoothing
                                 to the edges of the mask
        smoothing_strength (int, optional): strength of gaussian smoothing
        invert (bool, optional): flag for inverting the boolean selection

    Returns
        mask (ndarray): boolean mask with same dimensions of `img`
    """
    if type(poly) is _Canvas:
        poly = __np.array(poly.vert)

    if poly.shape[1] != 2 or len(poly.shape) != 2:
        raise ValueError(
            f"Polygon must have shape `[n, 2]` -> Got polygon with shape {poly.shape}"
        )

    # Create the mask from the polygon
    mask = __np.ones_like(img, dtype=float)
    rr, cc = polygon(poly[:, 0], poly[:, 1])
    mask[rr, cc] = 0

    if invert:
        mask = __np.invert(mask)

    # Smooth the edges of the mask
    if smooth:
        mask = gaussian_filter(mask, smoothing_strength)

    # Flip the maps along the axes of symmetry
    inds = [int(mask.shape[0]/2), int(mask.shape[1]/2)]
    mask[-inds[0]:, :inds[1]] = __np.fliplr(mask[-inds[0]:, -inds[1]:])
    mask[:inds[0], :] = __np.flipud(mask[-inds[0]:, :])

    return mask
    

#######################################################################
