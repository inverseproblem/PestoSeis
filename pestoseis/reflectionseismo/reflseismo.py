"""Functions to generate and process reflection data
"""

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

#######################################################################
#######################################################################

import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import CubicSpline 

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
    convo = np.convolve(refle,wavelet,mode="full")
    tarr = np.array([i*dt for i in range(convo.size)])

    return tarr,convo

#######################################################################

def calcreflectivity(density,vel,z,dt):
    """
    Compute the reflectivity series.
    """
    assert(vel.ndim==1)
    assert(density.ndim==1)

    impdz = density*vel
    reflz = (impdz[1:]-impdz[:-1])/(impdz[1:]+impdz[:-1])

    twtfromv = _depth2time(z,vel)
    twt = np.arange(0.0,twtfromv.max(),dt)
    npts = twt.size    

    refltwt = np.zeros(npts)
    for i in range(reflz.size):
        if np.abs(reflz[i])>0.0 :
            idx = np.argmin(np.abs(twtfromv[i]-twt))
            refltwt[idx] = reflz[i]
        
    return twt,refltwt

#######################################################################

def _depth2time(z,vel):
    """
     Convert depth to time.
    """
    assert(vel.ndim==1)
    assert(vel.size==z.size)

    npts=vel.size
    twt = np.zeros(npts)

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
    twt = np.array([i*dt for i in range(seisdata.shape[1])])  
    vmax = amplitudeclip*np.abs(seisdata).max()
    vmin = -vmax
    extent = [offset.min(),offset.max(),twt.max(),twt.min()]
    plt.imshow(seisdata.T,vmin=vmin,vmax=vmax,cmap=plt.cm.RdGy_r,
               extent=extent,aspect='auto',interpolation='bilinear')
    plt.colorbar()
    plt.xlabel('Offset [m]')
    plt.ylabel('TWT [s]')
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
    t = np.array([i*dt for i in range(data.shape[1])])        
    lwidth = 1.0 # line width
    ntraces = data.shape[0]

    if not isinstance(offset,np.ndarray) :
        ofs = np.array([float(i) for i in range(ntraces)])
        maxval=1.0/((abs(data).max()))
    else :
        ofs = offset
        maxseis = abs(data).max()
        maxoff = (abs(np.diff(offset).max()))
        maxval = maxoff/maxseis

    if scal!=None :
        maxval *= scal    

    if data.ndim==1:
        data=data.reshape(1,-1)

    for i in range(0,ntraces,skiptr):
        trace=ofs[i]+data[i,:]*maxval
        plt.plot(trace,t,color='black',linewidth=lwidth)
        if filltrace :
            plt.fill_betweenx(t,ofs[i],trace,where=(trace>=ofs[i]),color='black')
        #plt.fill_betweenx(t,ofs[i],trace,where=(trace<ofs[i]),color='red')
    ax = plt.gca()

    plt.xlim(ofs[0]-data[0,:].max()*maxval,ofs[-1]+data[-1,:].max())
    plt.ylim(t.min(),t.max())
    ax.invert_yaxis()

    if title!=None :
        plt.title(title)
    plt.xlabel('Offset [m]')
    plt.ylabel('TWT [s]')
    return

#######################################################################

def geometrical_spreading(seis,twt):
    """
    Apply geometrical spreading correction to a shotgather.

    Args:
        seis (ndarray): seismic shotgather, *rows* contain traces
        twt (ndarray): two-way traveltime 
            
    Returns
    ----------
        seis_gs (ndarray): seismograms corrected for geometrical spreading
    """
    ns = seis.shape[1]
    ntraces = seis.shape[0]
    
    seis_gs = seis.copy()
        
    for i in range(ntraces):
        ## twt**2 is a rough APPROXIMATION to geom. spreading...
        seis_gs[i,:] = seis_gs[i,:] * twt**2
        ## scale max amplitude to 1
        seis_gs[i,:] /= np.abs(seis_gs[i,:]).max()
        
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
    ----------
        seis_gs (ndarray): seismograms corrected with AGC
    """

    ns = seis.shape[1]
    ntraces = seis.shape[0];
    
    if rho>0:
        nw = int(np.ceil(w/rho))
    else:
        nw = w
        
    # select weight kernel type
    if type.lower() == 'gaussian':
        from scipy.stats import norm
        iw=np.linspace(-2,2,2*nw);
        weight_kernel = norm.pdf(iw)
        weight_kernel = weight_kernel/np.sum(weight_kernel)
    else:
        weight_kernel = np.ones(nw)
        weight_kernel = weight_kernel/np.sum(weight_kernel)
 
    seis_agc=np.zeros_like(seis)
    for i in range(ntraces):
        
        seis_abs = np.abs(seis[i,:])
        seis_weight = np.convolve(np.abs(seis[i,:]),weight_kernel, mode='same')

        if (seis_weight==0.0).any() :
            # print("agc(): Dividing by zero in AGC correction... ")
            seis_weight = np.where(seis_weight==0.0,1e-6,seis_weight)      
            
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
    ----------
        seisnmo (ndarray): seismograms NMO corrected
    """

    ntraces,nsmpl = seisdat.shape
    assert ntraces==offset.size
    seisnmo = np.zeros_like(seisdat)
    timearr = np.array([(i-1)*dt for i in range(nsmpl)])

    ## build the nmo array trace by trace
    for itr in range(ntraces):
        # compute the time for nmo for all time points of a single trace
        # velocity array is nx x nz so slice in depth is vel[?,:]
        tnmo = np.sqrt(timearr**2 + offset[itr]**2/velnmo**2) 
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
    ----------
        seisnew (ndarray): resampled trace
    """
    seisnew = np.zeros(torig.size)
    itp = CubicSpline(torig,seistr)
    # time outside bounds?
    outbou = (tnmo>torig.max()).nonzero()[0]
    if outbou.size>0:
        idx = np.min(outbou)
    else :
        idx = seisnew.size
    seisnew[:idx] = itp(tnmo[:idx])
    return seisnew

#######################################################################



