"""Animate result of 2D finite difference simulation

"""

#------------------------------------------------------------------------
#
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
# ######################################################################

import matplotlib.pyplot as plt
from matplotlib import animation
from matplotlib import gridspec

import h5py as h5

import numpy as np
import scipy.integrate as spi
import scipy.signal as sps

#######################################################################
# ######################################################################

def animateacousticwaves(inpfile,clipamplitude=0.1,showanim=True) :
    """
     Function to 'animate' the results of a 2D acoustic finite difference simulation. It produces (and saves to a file) an .mp4 movie.

    :parameter inpfile: input HDF5 file name
    :parameter clipamplitude: amplitude clipping factor
    :parameter showanim: show plot or not
    
    """
    ##============================
    every = 1 
    kind = "acoustic"
    fl=inpfile #"outdata/acoustic_snapshots.h5"

    ##=========================================

    h=h5.File(fl,"r")
    press = h["press"][:,:,:].copy()
    srctf = h["srctf"][:].copy()
    dt = h["dt"][()]
    dh = h["dh"][()] #[()]
    nx = h["nx"][()]
    nz = h["nz"][()]
    vel = h["vel"][:,:]
    snapevery = h["snapevery"][()]
    recs = h["recpos"][:,:].copy()
    h.close()
    x = press

    N=x.shape[2]
    pmax = clipamplitude*abs(x).max()
    pmin = -pmax
    cmap = plt.cm.RdGy_r #hot_r #gray_r #jet #RdGy_r #viridis

    ######################################################

    def updatefig_acou(n):
        reddot.set_data(t[(n-1)*snapevery],srctf[(n-1)*snapevery])
        sp2.set_title("Iter {}".format(n))
        wav.set_array(x[:,:,n].T)
        return (wav,reddot,sp2)

    ######################################################

    nt = srctf.size
    t = np.arange(0.0,dt*nt,dt)

    gs = gridspec.GridSpec(1, 4)
    gs.update(left=0.05, right=0.99, wspace=0.15,hspace=0.25)

    fig1 = plt.figure(figsize=(12,5))

    sp1 = plt.subplot(gs[0, 0])
    plt.title("Source time function")
    sep1 = plt.plot(t,srctf,'k')
    reddot, = plt.plot(0,srctf[0],'or')


    sp2 = plt.subplot(gs[0, 1:])
    plt.title("Amplitude scaling factor: {}".format(clipamplitude))
    # sp2 = plt.subplot(122)
    extent = [0.0,dh*(nx-1),dh*(nz-1),0.0]
    vp = plt.imshow(vel[:,:].T,cmap=plt.cm.jet,extent=extent,
                     interpolation="nearest", alpha=0.5)
    #srcpos, = plt.plot(  ,'^b')
    plt.scatter(recs[:,0],recs[:,1],marker="v",color="k")
    wav = plt.imshow(x[:,:,1].T,vmin=pmin,vmax=pmax,cmap=cmap,extent=extent,
                     interpolation="nearest", animated=True, alpha=0.8)
    plt.colorbar()
    #plt.tight_layout()

    ###
    ani = animation.FuncAnimation(fig1, updatefig_acou, frames=range(0,N,every), interval=150, blit=False)

    # Writer = animation.writers['ffmpeg']
    # mywriter = Writer(fps=15, metadata=dict(artist='PestoSeis'), bitrate=1800)
    mywriter = animation.FFMpegWriter()
    ani.save('animation_{}.mp4'.format(kind),dpi=96,writer=mywriter)


    ##################
    if showanim:
        plt.show()
        
    return ani


""
def animateelasticwaves(inpfile,showwhatela="VxVz",clipamplitude=0.1,showanim=True) :
    """
     Function to 'animate' the results of a 2D elastic finite difference simulation. It produces (and saves to a file) an .mp4 movie.

    :parameter inpfile: input HDF5 file name
    :parameter showwhatela: what to show, either 'VxVz' or 'PS'
    :parameter clipamplitude: amplitude clipping factor
    :parameter showanim: show plot or not
    
    """
    ##============================
    every = 1 
    kind = "elastic"
    fl=inpfile 

    ##=========================================

    h=h5.File(fl,"r")
    srctf = h["srctf"][:].copy()    
    dt = h["dt"][()]
    dh = h["dh"][()]
    nx = h["nx"][()]
    nz = h["nz"][()]
    vx = h["vx"][:,:].copy()
    vz = h["vz"][:,:].copy()
    lamb = h["lambda"][:,:].copy()
    mu = h["mu"][:,:].copy()
    rho = h["rho"][:,:].copy()
    snapevery = h["snapevery"][()]
    recs = h["recpos"][:,:].copy()
    h.close()

    # N=x.shape[2]
    # pmax = clipamplitude*abs(x).max()
    # pmin = -pmax
    cmap = plt.cm.RdGy_r #hot_r #gray_r #jet #RdGy_r #viridis

    if showwhatela=='PS':
        ## get displacement from velocity
        vxdetr = sps.detrend(vx,type='linear')
        vzdetr = sps.detrend(vz,type='linear')
        ux = spi.cumtrapz(vxdetr,dx=dt)
        uz = spi.cumtrapz(vzdetr,dx=dt)
        ## calculate derivatives
        tmpux = (ux[1:,1:,:]+ux[1:,:-1,:]+ux[:-1,1:,:]+ux[:-1,:-1,:])/4.0 # staggered grid...
        tmpuz = uz[:-1,:-1,:]
        ## divengence of displacement
        duxdx = np.diff(tmpux,axis=0)[:,:-1,:]/dh
        duzdz = np.diff(tmpuz,axis=1)[:-1,:,:]/dh
        Pw = duxdx[:,:,:] + duzdz[:,:,:]
        ## curl of displacement
        duxdz = np.diff(tmpux,axis=1)[:-1,:,:]/dh
        duzdx = np.diff(tmpuz,axis=0)[:,:-1,:]/dh
        Sw = duxdz[:,:,:] - duzdx[:,:,:]

        N=Pw.shape[2]
        pmax1 = clipamplitude * abs(Pw).max() #abs(x).max()
        pmin1 = -pmax1 #x.min()
        pmax2 = clipamplitude * abs(Pw).max() #abs(x).max()
        pmin2 = -pmax2 #x.min()

        title1 = "P waves"
        title2 = "S waves"
        data1 = Pw
        data2 = Sw

    elif showwhatela=='VxVz':
        data1 = vx
        data2 = vz
        title1 = "Vx"
        title2 = "Vz"

        N=vx.shape[2]
        pmax1 = clipamplitude * abs(vx).max() #abs(x).max()
        pmin1 = -pmax1 #x.min()
        pmax2 = clipamplitude * abs(vz).max() #abs(x).max()
        pmin2 = -pmax2 #x.min()

    elif showwhatela == "Vmag":
        data1 = np.sqrt(vx**2 + vz**2)
        title1 = "V Combined"
        data2 = 0 * vx
        title2 = ""

        N=vx.shape[2]
        pmax1 = clipamplitude * abs(data1).max() #abs(x).max()
        # pmin1 = -pmax1 #x.min()
        pmin1 = 0
        pmax2 = 1
        pmin2 = -pmax2 #x.min()

    ######################################################

    def updatefig_ela(n):
        reddot.set_data(t[(n-1)*snapevery],srctf[(n-1)*snapevery])
        sp2.set_title(title1+"  Snapshot: {} ".format(n))
        wav1.set_array(data1[:,:,n].T)
        sp3.set_title(title2+"  Snapshot: {} ".format(n))
        wav2.set_array(data2[:,:,n].T)
        return (wav1,sp2,wav2,sp3)

    ######################################################

    nt = srctf.size
    t = np.arange(0.0,dt*nt,dt)

    gs = gridspec.GridSpec(2, 4)
    gs.update(left=0.05, right=0.99, wspace=0.15,hspace=0.25)

    fig1 = plt.figure(figsize=(11,8))

    ## source-time function
    sp1 = plt.subplot(gs[0, 0])
    plt.title("Source time function")
    sep1 = plt.plot(t,srctf,'k')
    reddot, = plt.plot(0,srctf[0],'or')

    ## data1 waves
    sp2 = plt.subplot(gs[0, 1:])
    extent = [0.0,dh*(nx-1),dh*(nz-1),0.0]
    prop1 = plt.imshow(rho[:,:].T,cmap=plt.cm.jet,extent=extent,
                     interpolation="nearest", alpha=0.5)
    plt.scatter(recs[:,0],recs[:,1],marker="v",color="k")
    #plt.scatter(srcs[:,0],srcs[:,1],marker="o",color="k")
    wav1 = plt.imshow(data1[:,:,1].T,vmin=pmin1,vmax=pmax1,cmap=cmap,extent=extent,
                      interpolation="nearest", animated=True, alpha=0.7)
    plt.colorbar()

    ## data2 waves
    sp3 = plt.subplot(gs[1, 1:])
    extent = [0.0,dh*(nx-1),dh*(nz-1),0.0]
    prop2 = plt.imshow(rho[:,:].T,cmap=plt.cm.jet,extent=extent,
                     interpolation="nearest", alpha=0.5)
    plt.scatter(recs[:,0],recs[:,1],marker="v",color="k")
    #plt.scatter(srcs[:,0],srcs[:,1],marker="o",color="k")
    wav2 = plt.imshow(data2[:,:,1].T,vmin=pmin2,vmax=pmax2,cmap=cmap,extent=extent,
                      interpolation="nearest", animated=True, alpha=0.7)
    plt.colorbar()
    #plt.tight_layout(

    ###
    ani = animation.FuncAnimation(fig1, updatefig_ela, frames=range(0,N,every), interval=150, blit=False)

    # Writer = animation.writers['ffmpeg']
    # mywriter = Writer(fps=15, metadata=dict(artist='PestoSeis'), bitrate=1800)
    mywriter = animation.FFMpegWriter()
    ani.save('animation_{}.mp4'.format(kind),dpi=96,writer=mywriter)


    ##################
    if showanim:
        plt.show()

    return ani
