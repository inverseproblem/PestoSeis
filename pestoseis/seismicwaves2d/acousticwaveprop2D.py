
#------------------------------------------------------------------------
#
#    Copyright (C) 2019  Andrea Zunino 
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
import sys
import h5py as h5

import matplotlib.pyplot as plt
from matplotlib import animation
from matplotlib import gridspec

#############################################################################
#######################################################################

def animateacousticwaves(inpfile,clipamplitude=0.1) :
    """
     Function to 'animate' the results of a 2D acoustic finite difference simulation. It produces (and saves to a file) an .mp4 movie.

    :parameter inpfile: input HDF5 file name
    :parameter clipamplitude: amplitude clipping factor

    
    """
    ##============================
    every = 1 
    kind = "acoustic"
    fl=inpfile #"outdata/acoustic_snapshots.h5"

    ##=========================================

    h=h5.File(fl,"r")
    press = h["press"][:,:,:].copy()
    srctf = h["srctf"][:].copy()
    dt = h["dt"].value
    dh = h["dh"].value #[()]
    nx = h["nx"].value
    nz = h["nz"].value
    vel = h["vel"][:,:]
    snapevery = h["snapevery"].value
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
    plt.show()

    return ani

#############################################################################

def _bilinear_interp(f,hgrid, pt):
    """
     Bilinear interpolation (2D).
    """
    xreq=pt[0]
    yreq=pt[1]
    xh=xreq/hgrid
    yh=yreq/hgrid
    i=int(np.floor(xh)) # index starts from 0 so no +1
    j=int(np.floor(yh)) # index starts from 0 so no +1
    xd=xh-i # index starts from 0 so no +1
    yd=yh-j # index starts from 0 so no +1
    intval=f[i,j]*(1.0-xd)*(1.0-yd)+f[i+1,j]*(1.0-yd)*xd+f[i,j+1]*(1.0-xd)*yd+f[i+1,j+1]*xd*yd
    #print xreq,yreq,xh,yh,i,j,xd,yd    
    return intval

#############################################################################

def _initGaussboundcon(nptsgau=60 ) :
    
    ## Damping region size in grid points
    ##nptsgau = 50 #21
    ### decay = 0.015 and 20-30 pts is a tipical value found in literature...
    ### decay = 0.0053 and 60 pts was found on the web...
    decay = 0.25/nptsgau  ## 0.24/nptsgau  #0.48/nptsgau
    
    xdist = np.arange(1,nptsgau+1)
    damp = np.exp( -((decay * (nptsgau - xdist))**2))

    leftdp   = damp.copy()
    rightdp  = damp[::-1].copy()
    bottomdp = damp[::-1].copy()
    topdp    = damp.copy()

    # import matplotlib.pyplot as plt
    # plt.figure()
    # plt.plot(damp)
    # plt.show()

    return nptsgau,leftdp,rightdp,bottomdp,topdp

#########################################################################

def _calc_Kab_CPML(nptspml,gridspacing,dt,Npower,d0,
                       alpha_max_pml,K_max_pml,onwhere ) :

    # L = thickness of adsorbing layer
    if onwhere=="grdpts" :
        L = nptspml*gridspacing
        # distances 
        x = np.arange(0.0,nptspml*gridspacing,gridspacing)
    elif onwhere=="halfgrdpts" :
        L = nptspml*gridspacing
        # distances 
        x = np.arange(gridspacing/2.0,nptspml*gridspacing,gridspacing)
        
    d = d0 * (x/L)**Npower    
    alpha =  alpha_max_pml * (1.0 - (x/L)) # + 0.1 * alpha_max_pml ????

    K = 1.0 + (K_max_pml - 1.0) * (x/L)**Npower
    b = np.exp( - (d / K + alpha) * dt )
    a = d * (b-1.0)/(K*(d+K*alpha))
    
    return K,a,b


###################################################################

def solveacoustic2D( inpar, ijsrc, velmod, sourcetf, srcdomfreq, recpos, saveh5=True,
                     outfileh5="acoustic_snapshots.h5"):
    """
    Solve the acoustic wave equation in 2D using finite differences on a staggered grid. 
    Wrapper function for various boundary conditions.

    Args:
        inpar (dict): dictionary containing various input parameters

                      * inpar["ntimesteps"] (int) number of time steps
                      * inpar["nx"] (int) number of grid nodes in the x direction
                      * inpar["nz"] (int) number of grid nodes in the z direction
                      * inpar["dt"] (float) time step for the simulation
                      * inpar["dh"] (float) grid spacing (same in x and z)
                      * inpar["savesnapshot"] (bool) switch to save snapshots of the entire wavefield
                      * inpar["snapevery"] (int) save snapshots every "snapevery" iterations
                      * inpar["freesurface"] (bool) True for free surface boundary condition at the top, False for PML
                      * inpar["boundcond"] (string) Type of boundary conditions "PML","GaussTap" or "ReflBou"
        ijsrc (ndarray(int,int)): integers representing the position of the source on the grid
        velmod (ndarray(nx,nz)): two-dimensional velocity model
        sourcetf (ndarray): source time function
        srcdomfreq (float): source dominant frequency
        recpos (ndarray): position of the receivers, a two-column array [x,z]
        saveh5 (bool): whether to save results to HDF5 file or not
        outfileh5 (string): name of the output HDF5 file

    Returns:
        seism (ndarray): seismograms recorded at the receivers
        psave (ndarray): set of snapshots of the wavefield (if inpar["savesnapshot"]==True)

    """

    if inpar["boundcond"]=="PML" :
        seism,psave = _solveacouwaveq2D_CPML( inpar, ijsrc, velmod, sourcetf, srcdomfreq, recpos ) 

    elif inpar["boundcond"]=="GaussTap" :
        seism,psave = _solveacouwaveq2D_GaussTaper( inpar, ijsrc, velmod, sourcetf, srcdomfreq, recpos ) 

    elif inpar["boundcond"]=="ReflBou" :
        seism,psave = _solveacouwaveq2D_ReflBound( inpar, ijsrc, velmod, sourcetf, srcdomfreq, recpos )

    else :
        raise("Wrong boundary condition type")

    ##############################
    if saveh5:
        ## save stuff
        hf = h5.File(outfileh5,"w")
        if inpar["savesnapshot"]==True :
            hf["press"] = psave
            hf["snapevery"] = inpar["snapevery"]
        hf["seism"] = seism
        hf["vel"] = velmod
        hf["srctf"] = sourcetf
        hf["dh"] = inpar["dh"]
        hf["dt"] = inpar["dt"]
        hf["nx"] = inpar["nx"]
        hf["nz"] = inpar["nz"]
        hf["recpos"] = recpos
        hf["srcij"] = ijsrc
        hf.close()
        print("Saved acoustic simulation and parameters to ",outfileh5)

    return seism,psave

###################################################################

def _solveacouwaveq2D_CPML( inpar, ijsrc, vel, sourcetf, srcdomfreq, recpos ) :

    """
    Solve the acoustic wave equation in 2D using finite differences on a staggered grid. 
    CPML boundary conditions.

    Args:
        inpar (dict): dictionary containing various input parameters:

                      * inpar["ntimesteps"] (int) number of time steps
                      * inpar["nx"] (int) number of grid nodes in the x direction
                      * inpar["nz"] (int) number of grid nodes in the z direction
                      * inpar["dt"] (float) time step for the simulation
                      * inpar["dh"] (float) grid spacing (same in x and z)
                      * inpar["savesnapshot"] (bool) switch to save snapshots of the entire wavefield
                      * inpar["snapevery"] (int) save snapshots every "snapevery" iterations
                      * inpar["freesurface"] (bool) True for free surface boundary condition at the top, False for PML
                      * inpar["boundcond"] (string) Type of boundary conditions "PML" 
        ijsrc (ndarray(int,int)): integers representing the position of the source on the grid
        vel (ndarray(nx,nz)): two-dimensional velocity model
        sourcetf (ndarray): source time function
        srcdomfreq (float): source dominant frequency
        recpos (ndarray): position of the receivers, a two-column array [x,z]

    Returns:
        seism (ndarray): seismograms recorded at the receivers
        psave (ndarray): set of snapshots of the wavefield (if inpar["savesnapshot"]==True)

    """

    assert(inpar["boundcond"]=="PML")
    print("Starting ACOUSTIC solver with CPML boundary condition.")
    
    ##
    ## The second-order staggered-grid formulation of Madariaga (1976) and Virieux (1986) is used:
    ##
    ##          p          dp/dx       p   
    ##            +---------+---------+ ---> x  
    ##            |         |         |
    ##            |         |         |
    ##            |         |         |
    ##            |         |         |
    ##            |         |         |
    ##      dp/dz +---------+         |
    ##            |                   |
    ##            |                   |
    ##            |                   |
    ##            |                   |
    ##            |                   |
    ##            +-------------------+  
    ##           p                     p
    ##            |
    ##            |
    ##            \/ z
    ##


    #############################
    dh = inpar["dh"]
    dx = dh
    dz = dh
    dt = inpar["dt"]

    ##############################
    ## Check stability criterion
    ##############################
    maxvp = vel.max()
    print(" Stability criterion, CFL number:",(maxvp*dt*np.sqrt(1/dx**2+1/dz**2)))
    assert(maxvp*dt*np.sqrt(1/dx**2 + 1/dz**2) < 1.0)
    

    ##############################
    #   Parameters
    ##############################
    ## pml = c-pml, else reflsurf
    ## freesur = free surf for top
    if inpar["freesurface"]==True :
        print(" * Free surface at the top *")
        freeboundtop=True
    else :
        print(" * Absorbing boundary at the top *")
        freeboundtop=False

    ##---------------------------
    f0 = srcdomfreq # source
    
    ## keep it simple for now...
    nx = inpar["nx"]
    nz = inpar["nz"]
    
    ##############################
    #   Parameters for PML
    ##############################
    nptspml_x = 21 ## 21 ref coeff set for 21 absorbing nodes
    nptspml_z = 21
    assert(nptspml_x<ijsrc[0]<(nx-1-nptspml_x))
    assert(ijsrc[1]<(nz-1-nptspml_z))
    if freeboundtop==False:
        assert(nptspml_z<ijsrc[1])

    print((" Size of PML layers in grid points: {} in x and {} in z".format(nptspml_x,nptspml_z)))
    
    Npower = 2.0
    assert Npower >= 1
    
    K_max_pml = 1.0
    alpha_max_pml = 2.0*np.pi*(f0/2.0) 
    
    # reflection coefficient (INRIA report section 6.1)
    # http://hal.inria.fr/docs/00/07/32/19/PDF/RR-3471.pdf
    Rcoef = 0.0001  # for 20 nodes   ## 0.001 for a PML thickness of 10 nodes
       
    # thickness of the PML layer in meters
    thickness_pml_x = nptspml_x * dx
    thickness_pml_z = nptspml_z * dz
    
    # compute d0 from INRIA report section 6.1 http://hal.inria.fr/docs/00/07/32/19/PDF/RR-3471.pdf
    d0_x = - (Npower + 1) * maxvp * np.log(Rcoef) / (2.0 * thickness_pml_x)
    d0_z = - (Npower + 1) * maxvp * np.log(Rcoef) / (2.0 * thickness_pml_z)
    
    ##############################
    #   Damping parameters
    ##############################    
    # --- damping in the x direction ---
    # assuming the number of grid points for PML is the same on 
    #    both sides    
    # damping profile at the grid points
    K_xpml,a_xpml,b_xpml = _calc_Kab_CPML(nptspml_x,dh,dt,Npower,d0_x,alpha_max_pml,K_max_pml,"grdpts")
    # damping profile at half the grid points
    K_xpml_half,a_xpml_half,b_xpml_half = _calc_Kab_CPML(nptspml_x,dh,dt,Npower,d0_x,alpha_max_pml,K_max_pml,"halfgrdpts")
    
    # --- damping in the z direction ---
    # assuming the number of grid points for PML is the same on
    # both sides    
    # damping profile at the grid points
    K_zpml,a_zpml,b_zpml = _calc_Kab_CPML(nptspml_z,dh,dt,Npower,d0_z,alpha_max_pml,K_max_pml,"grdpts")
    # damping profile at half the grid points
    K_zpml_half,a_zpml_half,b_zpml_half = _calc_Kab_CPML(nptspml_z,dh,dt,Npower,d0_z,alpha_max_pml,K_max_pml,"halfgrdpts")

    
    #######################
    #  x direction
    #######################
    K_x     = np.ones(nx)
    a_x     = np.zeros(nx)
    b_x     = np.ones(nx)
    K_x_half = np.ones(nx)
    a_x_half = np.zeros(nx)
    b_x_half = np.ones(nx)
    
    # reverse coefficients to get increasing damping away from inner model
    # left boundary
    K_x[:nptspml_x]  = K_xpml[::-1] #[end:-1:1]
    a_x[:nptspml_x]  = a_xpml[::-1]
    b_x[:nptspml_x]  = b_xpml[::-1]
    # half stuff...
    # reverse coefficients to get increasing damping away from inner model
    #  One less element on left boundary (end-1...)
    #    because of the staggered grid
    K_x_half[:nptspml_x-1] = K_xpml_half[-2::-1] #[end:-1:1]
    a_x_half[:nptspml_x-1] = a_xpml_half[-2::-1]
    b_x_half[:nptspml_x-1] = b_xpml_half[-2::-1]
    
    # right boundary 
    rightpml = nx-nptspml_x  # julia: nx-nptspml_x+1 
    K_x[rightpml:]   = K_xpml
    a_x[rightpml:]   = a_xpml
    b_x[rightpml:]   = b_xpml
    # half
    K_x_half[rightpml:]   = K_xpml_half
    a_x_half[rightpml:]   = a_xpml_half
    b_x_half[rightpml:]   = b_xpml_half
    
    #######################
    #  z direction
    #######################
    K_z     = np.ones(nz)
    a_z     = np.zeros(nz)
    b_z     = np.zeros(nz)
    # half
    K_z_half   = np.ones(nz)
    a_z_half   = np.zeros(nz)
    b_z_half   = np.zeros(nz)
    
    # bottom 
    bottompml=nz-nptspml_z  # julia: nz-nptspml_z+1
    K_z[bottompml:]  = K_zpml
    a_z[bottompml:]  = a_zpml
    b_z[bottompml:]  = b_zpml
    # half
    K_z_half[bottompml:]  = K_zpml_half
    a_z_half[bottompml:]  = a_zpml_half
    b_z_half[bottompml:]  = b_zpml_half
    
    # if PML also on top of model...
    if freeboundtop==False :
        # on grid
        K_z[:nptspml_z]  = K_zpml[::-1]
        a_z[:nptspml_z]  = a_zpml[::-1]
        b_z[:nptspml_z]  = b_zpml[::-1]    
        # half
        #  One less element on top boundary (end-1...)
        #    because of the staggered grid
        K_z_half[:nptspml_z-1]  = K_zpml_half[-2::-1]
        a_z_half[:nptspml_z-1]  = a_zpml_half[-2::-1]
        b_z_half[:nptspml_z-1]  = b_zpml_half[-2::-1]

    #####################################################

    ## Arrays to export snapshots
    if inpar["savesnapshot"]==True :
        ntsave = inpar["ntimesteps"]//inpar["snapevery"]
        psave = np.zeros((inpar["nx"],inpar["nz"],ntsave+1))
        tsave=1

    ## Arrays to return seismograms
    nrecs = recpos.shape[0]
    receiv = np.zeros((inpar["ntimesteps"],nrecs))

    # Source time function
    #lensrctf = sourcetf.size
    isrc = ijsrc[0]
    jsrc = ijsrc[1]
    
    ## Initialize arrays
    dpdx = np.zeros((nx,nz))
    dpdz = np.zeros((nx,nz))
    pcur = np.zeros((nx,nz))
    pold = np.zeros((nx,nz))
    pnew = np.zeros((nx,nz))
    
    ##--------------------------------------------
    # PML arrays
    # Arrays with size of PML areas would be sufficient and save memory,
    #   however allocating arrays with same size than model simplifies
    #   the code in the loops
    psi_x = np.zeros((nx,nz))
    psi_z = np.zeros((nx,nz))
    xi_x = np.zeros((nx,nz))
    xi_z = np.zeros((nx,nz))

    ##---------------------------------------------------
    ## to make the Numpy broadcast Hadamard product correct along x
    a_x = a_x.reshape(-1,1)
    b_x = b_x.reshape(-1,1)
    K_x = K_x.reshape(-1,1)
    a_x_half = a_x_half.reshape(-1,1)
    b_x_half = b_x_half.reshape(-1,1)
    K_x_half = K_x_half.reshape(-1,1)
    
    ################################
    fact = vel**2 * (dt**2/dh**2)

       
    ## time loop
    print((" Time step dt: {}".format(dt)))

    for t in range(inpar["ntimesteps"]) :

        if t%100==0 :
            sys.stdout.write("\r Time step {} of {}".format(t,inpar["ntimesteps"]))
            sys.stdout.flush()


        if freeboundtop==True :
            ##-------------------------------------------
            ## free surface boundary cond.
            pcur[:,0] = 0.0
            #dpdx = pcur[2:,0]-2.0*pcur[1:-1,0]+pcur[:-2,0]
            #dpdz = pcur[1:-1,1]-2.0*pcur[1:-1,0] +0.0 #pcur[1:-1,:-2]
            pnew[:,0] = 0.0 #2.0*pcur[1:-1,0] -pold[1:-1,0] + fact[1:-1,0]*(dpdx) + fact[1:-1,0]*(dpdz) 
            ##-------------------------------------------


        ##=================================================
        ## second order stencil
        ##-----------------------------------------
        ## first derivatives
        dpdx = pcur[2:,1:-1]-pcur[1:-1,1:-1] # forward 
        dpdz = pcur[1:-1,2:]-pcur[1:-1,1:-1]
        psi_x[1:-1,1:-1] = b_x_half[1:-1] / K_x_half[1:-1] * psi_x[1:-1,1:-1] + a_x_half[1:-1]*dpdx
        psi_z[1:-1,1:-1] = b_z_half[1:-1] / K_z_half[1:-1] * psi_z[1:-1,1:-1] + a_z_half[1:-1]*dpdz

        # second derivatives 
        dpdx2 = pcur[2:,1:-1]-2.0*pcur[1:-1,1:-1]+pcur[:-2,1:-1]
        dpdz2 = pcur[1:-1,2:]-2.0*pcur[1:-1,1:-1]+pcur[1:-1,:-2]

        dpsidx = psi_x[1:-1,1:-1] - psi_x[:-2,1:-1]
        dpsidz = psi_z[1:-1,1:-1] - psi_z[1:-1,:-2]
        xi_x[1:-1,1:-1] = b_x[1:-1] / K_x[1:-1] * xi_x[1:-1,1:-1] + a_x[1:-1] * (dpdx2 + dpsidx)
        xi_z[1:-1,1:-1] = b_z[1:-1] / K_z[1:-1] * xi_z[1:-1,1:-1] + a_z[1:-1] * (dpdz2 + dpsidz)
        
        damp = fact[1:-1,1:-1] * (dpsidx + dpsidz + xi_x[1:-1,1:-1] + xi_z[1:-1,1:-1])
        

        # update pressure
        pnew[1:-1,1:-1] = 2.0*pcur[1:-1,1:-1] -pold[1:-1,1:-1] + fact[1:-1,1:-1]*(dpdx2 + dpdz2) + damp

        

        # inject source
        pnew[isrc,jsrc] = pnew[isrc,jsrc] + dt**2*sourcetf[t]

        # assign the new pold and pcur
        pold[:,:] = pcur[:,:]
        pcur[:,:] = pnew[:,:]

        ##=================================================
        
        ##### receivers
        for r in range(nrecs) :
            rec_press = _bilinear_interp(pcur,dh,recpos[r,:])
            receiv[t,r] = rec_press
        
        #### save snapshots
        if (inpar["savesnapshot"]==True) and (t%inpar["snapevery"]==0) :
            psave[:,:,tsave] = pcur
            tsave=tsave+1

    ##================================
    print(" ")
    if inpar["savesnapshot"]==False :
        psave = None
        
    return receiv,psave
  
###################################################################

def _solveacouwaveq2D_Vp_density_CPML( inpar, ijsrc, Vp, density, sourcetf, srcdomfreq, recpos ) :
    """
    Solve the acoustic wave equation in 2D using finite differences on a staggered grid. 
    Velocity and density as input parameters. CPML boundary conditions.

    Args:
        inpar (dict): dictionary containing various input parameters:

                      * inpar["ntimesteps"] (int) number of time steps
                      * inpar["nx"] (int) number of grid nodes in the x direction
                      * inpar["nz"] (int) number of grid nodes in the z direction
                      * inpar["dt"] (float) time step for the simulation
                      * inpar["dh"] (float) grid spacing (same in x and z)
                      * inpar["savesnapshot"] (bool) switch to save snapshots of the entire wavefield
                      * inpar["snapevery"] (int) save snapshots every "snapevery" iterations
                      * inpar["freesurface"] (bool) True for free surface boundary condition at the top, False for PML
                      * inpar["boundcond"] (string) Type of boundary conditions "PML"
        ijsrc (ndarray(int,int)): integers representing the position of the source on the grid
        density (ndarray(nx,nz)): two-dimensional density model
        Vp (ndarray(nx,nz)): two-dimensional velocity model
        sourcetf (ndarray): source time function
        srcdomfreq (float): source dominant frequency
        recpos (ndarray): position of the receivers, a two-column array [x,z]

    Returns:
        seism (ndarray): seismograms recorded at the receivers
        psave (ndarray): set of snapshots of the wavefield (if inpar["savesnapshot"]==True)

    """
    assert(inpar["boundcond"]=="PML")
    print("\n WARNING: UNTESTED!!!!! \nStarting ACOUSTIC 'Vp_density' solver with CPML boundary condition.")
    
    ##
    ## The second-order staggered-grid formulation of Madariaga (1976) and Virieux (1986) is used:
    ##
    ##          p          dp/dx       p   
    ##            +---------+---------+ ---> x  
    ##            |         |         |
    ##            |         |         |
    ##            |         |         |
    ##            |         |         |
    ##            |         |         |
    ##      dp/dz +---------+         |
    ##            |                   |
    ##            |                   |
    ##            |                   |
    ##            |                   |
    ##            |                   |
    ##            +-------------------+  
    ##           p                     p
    ##            |
    ##            |
    ##            \/ z
    ##


    #############################
    dh = inpar["dh"]
    dx = dh
    dz = dh
    dt = inpar["dt"]
    kappa = density * Vp**2

    ##############################
    ## Check stability criterion
    ##############################
    maxvp = Vp.max()
    print(" Stability criterion, CFL number:",(maxvp*dt*np.sqrt(1/dx**2+1/dz**2)))
    assert(maxvp*dt*np.sqrt(1/dx**2 + 1/dz**2) < 1.0)
    

    ##############################
    #   Parameters
    ##############################
    ## pml = c-pml, else reflsurf
    ## freesur = free surf for top
    if inpar["freesurface"]==True :
        print(" * Free surface at the top *")
        freeboundtop=True
    else :
        print(" * Absorbing boundary at the top *")
        freeboundtop=False

    ##---------------------------
    f0 = srcdomfreq # source
    
    ## keep it simple for now...
    nx = inpar["nx"]
    nz = inpar["nz"]
    
    ##############################
    #   Parameters for PML
    ##############################
    nptspml_x = 21 ## 21 ref coeff set for 21 absorbing nodes
    nptspml_z = 21
    assert(nptspml_x<ijsrc[0]<(nx-1-nptspml_x))
    assert(ijsrc[1]<(nz-1-nptspml_z))
    if freeboundtop==False:
        assert(nptspml_z<ijsrc[1])

    print((" Size of PML layers in grid points: {} in x and {} in z".format(nptspml_x,nptspml_z)))
    
    Npower = 2.0
    assert Npower >= 1
    
    K_max_pml = 1.0
    alpha_max_pml = 2.0*np.pi*(f0/2.0) 
    
    # reflection coefficient (INRIA report section 6.1)
    # http://hal.inria.fr/docs/00/07/32/19/PDF/RR-3471.pdf
    Rcoef = 0.0001  # 0.0001 for 20 nodes   ## 0.001 for a PML thickness of 10 nodes
       
    # thickness of the PML layer in meters
    thickness_pml_x = nptspml_x * dx
    thickness_pml_z = nptspml_z * dz
    
    # compute d0 from INRIA report section 6.1 http://hal.inria.fr/docs/00/07/32/19/PDF/RR-3471.pdf
    d0_x = - (Npower + 1) * maxvp * np.log(Rcoef) / (2.0 * thickness_pml_x)
    d0_z = - (Npower + 1) * maxvp * np.log(Rcoef) / (2.0 * thickness_pml_z)
    
    ##############################
    #   Damping parameters
    ##############################    
    # --- damping in the x direction ---
    # assuming the number of grid points for PML is the same on 
    #    both sides    
    # damping profile at the grid points
    K_xpml,a_xpml,b_xpml = _calc_Kab_CPML(nptspml_x,dh,dt,Npower,d0_x,alpha_max_pml,K_max_pml,"grdpts")
    # damping profile at half the grid points
    K_xpml_half,a_xpml_half,b_xpml_half = _calc_Kab_CPML(nptspml_x,dh,dt,Npower,d0_x,alpha_max_pml,K_max_pml,"halfgrdpts")
    
    # --- damping in the z direction ---
    # assuming the number of grid points for PML is the same on
    # both sides    
    # damping profile at the grid points
    K_zpml,a_zpml,b_zpml = _calc_Kab_CPML(nptspml_z,dh,dt,Npower,d0_z,alpha_max_pml,K_max_pml,"grdpts")
    # damping profile at half the grid points
    K_zpml_half,a_zpml_half,b_zpml_half = _calc_Kab_CPML(nptspml_z,dh,dt,Npower,d0_z,alpha_max_pml,K_max_pml,"halfgrdpts")

    
    #######################
    #  x direction
    #######################
    K_x     = np.ones(nx)
    a_x     = np.zeros(nx)
    b_x     = np.ones(nx)
    K_x_half = np.ones(nx)
    a_x_half = np.zeros(nx)
    b_x_half = np.ones(nx)
    
    # reverse coefficients to get increasing damping away from inner model
    # left boundary
    K_x[:nptspml_x]  = K_xpml[::-1] #[end:-1:1]
    a_x[:nptspml_x]  = a_xpml[::-1]
    b_x[:nptspml_x]  = b_xpml[::-1]
    # half stuff...
    # reverse coefficients to get increasing damping away from inner model
    #  One less element on left boundary (end-1...)
    #    because of the staggered grid
    K_x_half[:nptspml_x-1] = K_xpml_half[-2::-1] #[end:-1:1]
    a_x_half[:nptspml_x-1] = a_xpml_half[-2::-1]
    b_x_half[:nptspml_x-1] = b_xpml_half[-2::-1]
    
    # right boundary 
    rightpml = nx-nptspml_x  # julia: nx-nptspml_x+1 
    K_x[rightpml:]   = K_xpml
    a_x[rightpml:]   = a_xpml
    b_x[rightpml:]   = b_xpml
    # half
    K_x_half[rightpml:]   = K_xpml_half
    a_x_half[rightpml:]   = a_xpml_half
    b_x_half[rightpml:]   = b_xpml_half
    
    #######################
    #  z direction
    #######################
    K_z     = np.ones(nz)
    a_z     = np.zeros(nz)
    b_z     = np.zeros(nz)
    # half
    K_z_half   = np.ones(nz)
    a_z_half   = np.zeros(nz)
    b_z_half   = np.zeros(nz)
    
    # bottom 
    bottompml=nz-nptspml_z  # julia: nz-nptspml_z+1
    K_z[bottompml:]  = K_zpml
    a_z[bottompml:]  = a_zpml
    b_z[bottompml:]  = b_zpml
    # half
    K_z_half[bottompml:]  = K_zpml_half
    a_z_half[bottompml:]  = a_zpml_half
    b_z_half[bottompml:]  = b_zpml_half
    
    # if PML also on top of model...
    if freeboundtop==False :
        # on grid
        K_z[:nptspml_z]  = K_zpml[::-1]
        a_z[:nptspml_z]  = a_zpml[::-1]
        b_z[:nptspml_z]  = b_zpml[::-1]    
        # half
        #  One less element on top boundary (end-1...)
        #    because of the staggered grid
        K_z_half[:nptspml_z-1]  = K_zpml_half[-2::-1]
        a_z_half[:nptspml_z-1]  = a_zpml_half[-2::-1]
        b_z_half[:nptspml_z-1]  = b_zpml_half[-2::-1]

    #####################################################

    ## Arrays to export snapshots
    if inpar["savesnapshot"]==True :
        ntsave = inpar["ntimesteps"]//inpar["snapevery"]
        psave = np.zeros((inpar["nx"],inpar["nz"],ntsave+1))
        tsave=1

    ## Arrays to return seismograms
    nrecs = recpos.shape[0]
    receiv = np.zeros((inpar["ntimesteps"],nrecs))

    # Source time function
    #lensrctf = sourcetf.size
    isrc = ijsrc[0]
    jsrc = ijsrc[1]
    
    ## Initialize arrays
    dpdx = np.zeros((nx,nz))
    dpdz = np.zeros((nx,nz))
    pcur = np.zeros((nx,nz))
    pold = np.zeros((nx,nz))
    pnew = np.zeros((nx,nz))
    
    ##--------------------------------------------
    # PML arrays
    # Arrays with size of PML areas would be sufficient and save memory,
    #   however allocating arrays with same size than model simplifies
    #   the code in the loops
    psi_x = np.zeros((nx,nz))
    psi_z = np.zeros((nx,nz))
    xi_x = np.zeros((nx,nz))
    xi_z = np.zeros((nx,nz))

    ##---------------------------------------------------
    ## to make the Numpy broadcast Hadamard product correct along x
    a_x = a_x.reshape(-1,1)
    b_x = b_x.reshape(-1,1)
    K_x = K_x.reshape(-1,1)
    a_x_half = a_x_half.reshape(-1,1) 
    b_x_half = b_x_half.reshape(-1,1)
    K_x_half = K_x_half.reshape(-1,1)
    
    ################################
    
    fact = dt**2/dh**2
    
    density_ihalf = (density[1:,:]+density[:-1,:])/2.0
    density_jhalf = (density[:,1:]+density[:,:-1])/2.0

    print(density_ihalf.shape,density_jhalf.shape) 

    ## time loop
    print((" Time step dt: {}".format(dt)))

    for t in range(inpar["ntimesteps"]) :

        if t%100==0 :
            sys.stdout.write("\r Time step {} of {}".format(t,inpar["ntimesteps"]))
            sys.stdout.flush()


        if freeboundtop==True :
            ##-------------------------------------------
            ## free surface boundary cond.
            pcur[:,0] = 0.0
            #dpdx = pcur[2:,0]-2.0*pcur[1:-1,0]+pcur[:-2,0]
            #dpdz = pcur[1:-1,1]-2.0*pcur[1:-1,0] +0.0 #pcur[1:-1,:-2]
            pnew[:,0] = 0.0 #2.0*pcur[1:-1,0] -pold[1:-1,0] + fact[1:-1,0]*(dpdx) + fact[1:-1,0]*(dpdz) 
            ##-------------------------------------------


        ##=================================================
        ## second order stencil
        ##  ACOUSTIC 'Vp_density'
        ##-----------------------------------------
        ## first derivatives
        dpdx = (pcur[2:,1:-1]-pcur[1:-1,1:-1]) / density_ihalf[:-1,1:-1] # forward 
        dpdz = (pcur[1:-1,2:]-pcur[1:-1,1:-1]) / density_jhalf[1:-1,:-1]
        psi_x[1:-1,1:-1] = b_x_half[1:-1] / K_x_half[1:-1] * psi_x[1:-1,1:-1] + a_x_half[1:-1]*dpdx
        psi_z[1:-1,1:-1] = b_z_half[1:-1] / K_z_half[1:-1] * psi_z[1:-1,1:-1] + a_z_half[1:-1]*dpdz

        # second derivatives 
        dpdx2 = pcur[2:,1:-1]-2.0*pcur[1:-1,1:-1]+pcur[:-2,1:-1]
        dpdz2 = pcur[1:-1,2:]-2.0*pcur[1:-1,1:-1]+pcur[1:-1,:-2]

        dpsidx = psi_x[1:-1,1:-1] - psi_x[:-2,1:-1]
        dpsidz = psi_z[1:-1,1:-1] - psi_z[1:-1,:-2]
        xi_x[1:-1,1:-1] = b_x[1:-1] / K_x[1:-1] * xi_x[1:-1,1:-1] + a_x[1:-1] * (dpdx2 + dpsidx)
        xi_z[1:-1,1:-1] = b_z[1:-1] / K_z[1:-1] * xi_z[1:-1,1:-1] + a_z[1:-1] * (dpdz2 + dpsidz)
        
        damp = fact * (dpsidx + dpsidz + xi_x[1:-1,1:-1] + xi_z[1:-1,1:-1])
        

        # update pressure
        pnew[1:-1,1:-1] = 2.0*pcur[1:-1,1:-1] -pold[1:-1,1:-1] + fact*(dpdx2+dpdz2)*kappa[1:-1,1:-1] + damp
        

        # inject source
        pnew[isrc,jsrc] = pnew[isrc,jsrc] + dt**2*sourcetf[t]

        # assign the new pold and pcur
        pold[:,:] = pcur[:,:]
        pcur[:,:] = pnew[:,:]

        ##=================================================
        
        ##### receivers
        for r in range(nrecs) :
            rec_press = _bilinear_interp(pcur,dh,recpos[r,:])
            receiv[t,r] = rec_press
        
        #### save snapshots
        if (inpar["savesnapshot"]==True) and (t%inpar["snapevery"]==0) :
            psave[:,:,tsave] = pcur
            tsave=tsave+1

    ##================================
    print(" ")
    if inpar["savesnapshot"]==False :
        psave = None
        
    return receiv,psave
  
#########################################################

def _solveacouwaveq2D_ReflBound( inpar, ijsrc, vel, sourcetf, srcdomfreq, recpos ) :
    """
    Solve the acoustic wave equation in 2D using finite differences on a staggered grid. 
    Reflective boundary conditions.

    Args:
        inpar (dict): dictionary containing various input parameters:

                      * inpar["ntimesteps"] (int) number of time steps
                      * inpar["nx"] (int) number of grid nodes in the x direction
                      * inpar["nz"] (int) number of grid nodes in the z direction
                      * inpar["dt"] (float) time step for the simulation
                      * inpar["dh"] (float) grid spacing (same in x and z)
                      * inpar["savesnapshot"] (bool) switch to save snapshots of the entire wavefield
                      * inpar["snapevery"] (int) save snapshots every "snapevery" iterations
                      * inpar["freesurface"] (bool) True for free surface boundary condition at the top, False for PML
                      * inpar["boundcond"] (string) Type of boundary conditions "ReflBou" 
        ijsrc (ndarray(int,int)): integers representing the position of the source on the grid
        vel (ndarray(nx,nz)): two-dimensional velocity model
        sourcetf (ndarray): source time function
        srcdomfreq (float): source dominant frequency
        recpos (ndarray): position of the receivers, a two-column array [x,z]

    Returns:
        seism (ndarray): seismograms recorded at the receivers
        psave (ndarray): set of snapshots of the wavefield (if inpar["savesnapshot"]==True)

    """

    assert(inpar["boundcond"]=="ReflBou")
    print("Starting ACOUSTIC solver with reflective boundaries all around.")

    #############################
    dh = inpar["dh"]
    dx = dh
    dz = dh
    dt = inpar["dt"]
    
    ##############################
    ## Check stability criterion
    ##############################
    maxvp = vel.max()
    print(" Stability criterion, CFL number:",(maxvp*dt*np.sqrt(1/dx**2+1/dz**2)))
    assert(maxvp*dt*np.sqrt(1/dx**2 + 1/dz**2) < 1.0)
    

    ##############################
    #   Parameters
    ##############################
    f0 = srcdomfreq # source
    
    ## keep it simple for now...
    nx = inpar["nx"]
    nz = inpar["nz"]
    
    #####################################################

    ## Arrays to export snapshots
    if inpar["savesnapshot"]==True :
        ntsave = inpar["ntimesteps"]//inpar["snapevery"]
        psave = np.zeros((inpar["nx"],inpar["nz"],ntsave+1))
        tsave=1

    ## Arrays to return seismograms
    nrecs = recpos.shape[0]
    receiv = np.zeros((inpar["ntimesteps"],nrecs))

    # Source time function
    lensrctf = sourcetf.size
    isrc = ijsrc[0]
    jsrc = ijsrc[1]
    
    ## Initialize arrays
    dpdx = np.zeros((nx,nz))
    dpdz = np.zeros((nx,nz))
    pcur = np.zeros((nx,nz))
    pold = np.zeros((nx,nz))
    pnew = np.zeros((nx,nz))

    ################################
    fact = vel**2 * (dt**2/dh**2)
       
    ## time loop
    print((" Time step dt: {}".format(dt)))

    for t in range(inpar["ntimesteps"]) :

        if t%100==0 :
            sys.stdout.write("\r Time step {} of {}".format(t,inpar["ntimesteps"]))
            sys.stdout.flush()

        ##=================================================
        ## second order stencil
        dpdx = pcur[2:,1:-1]-2.0*pcur[1:-1,1:-1]+pcur[:-2,1:-1]
        dpdz = pcur[1:-1,2:]-2.0*pcur[1:-1,1:-1]+pcur[1:-1,:-2]
        
        # update pressure
        pnew[1:-1,1:-1] = 2.0*pcur[1:-1,1:-1] -pold[1:-1,1:-1] + fact[1:-1,1:-1]*(dpdx) + fact[1:-1,1:-1]*(dpdz) 
        
        # inject source
        pnew[isrc,jsrc] = pnew[isrc,jsrc] + dt**2*sourcetf[t]

        # assign the new pold and pcur
        pold[:,:] = pcur[:,:]
        pcur[:,:] = pnew[:,:]

        ##=================================================
        
        ##### receivers
        for r in range(nrecs) :
            rec_press = _bilinear_interp(pcur,dh,recpos[r,:])
            receiv[t,r] = rec_press
        
        #### save snapshots
        if (inpar["savesnapshot"]==True) and (t%inpar["snapevery"]==0) :
            psave[:,:,tsave] = pcur
            tsave=tsave+1

    ##================================
    print(" ")

    if inpar["savesnapshot"]==False :
        psave = None
        
    return receiv,psave


#########################################################

def _solveacouwaveq2D_GaussTaper( inpar, ijsrc, vel, sourcetf, srcdomfreq, recpos ) :
    """
    Solve the acoustic wave equation in 2D using finite differences on a staggered grid. 
    Gaussian taper boundary conditions.

    Args:
        inpar (dict): dictionary containing various input parameters:\n

                      * inpar["ntimesteps"] (int) number of time steps\n
                      * inpar["nx"] (int) number of grid nodes in the x direction\n
                      * inpar["nz"] (int) number of grid nodes in the z direction\n
                      * inpar["dt"] (float) time step for the simulation\n
                      * inpar["dh"] (float) grid spacing (same in x and z)\n
                      * inpar["savesnapshot"] (bool) switch to save snapshots of the entire wavefield\n
                      * inpar["snapevery"] (int) save snapshots every "snapevery" iterations\n
                      * inpar["freesurface"] (bool) True for free surface boundary condition at the top, False for PML\n
                      * inpar["boundcond"] (string) Type of boundary conditions "GaussTap"\n 
        ijsrc (ndarray(int,int)): integers representing the position of the source on the grid
        vel (ndarray(nx,nz)): two-dimensional velocity model
        sourcetf (ndarray): source time function
        srcdomfreq (float): source dominant frequency
        recpos (ndarray): position of the receivers, a two-column array [x,z]

    Returns:
        seism (ndarray): seismograms recorded at the receivers
        psave (ndarray): set of snapshots of the wavefield (if inpar["savesnapshot"]==True)

    """
    assert(inpar["boundcond"]=="GaussTap")
    print("Starting ACOUSTIC solver with Gaussian taper boundary condition.")

    #############################
    dh = inpar["dh"]
    dx = dh
    dz = dh
    dt = inpar["dt"]
    
    ##############################
    ## Check stability criterion
    ##############################
    maxvp = vel.max()
    counum = maxvp * dt*  np.sqrt(1/dh**2 + 1/dh**2)
    print("Courant number: ",counum)
    assert(counum < 1.0)

    ##############################
    #   Parameters
    ##############################
    if inpar["freesurface"]==True :
        print("freesurface at the top")
        freeboundtop=True
    else :
        freeboundtop=False

    ##----------------
    f0 = srcdomfreq # source
    
    ## keep it simple for now...
    nx = inpar["nx"]
    nz = inpar["nz"]
    
    ## Gaussian taper
    nptsgau = 50
    assert(nptsgau<ijsrc[0]<(nx-1-nptsgau))
    assert(ijsrc[1]<(nz-1-nptsgau))
    if freeboundtop==False:
        assert(nptsgau<ijsrc[1])

        
    nptsgau,gtleft,gtright,gtbottom,gttop = _initGaussboundcon(nptsgau=nptsgau )

    print((" Size of GaussTaper layers in grid points: {} in both x and z".format(nptsgau)))

    #####################################################
    
    ## Arrays to export snapshots
    if inpar["savesnapshot"]==True :
        ntsave = inpar["ntimesteps"]//inpar["snapevery"]
        psave = np.zeros((inpar["nx"],inpar["nz"],ntsave+1))
        tsave=1

    ## Arrays to return seismograms
    nrecs = recpos.shape[0]
    receiv = np.zeros((inpar["ntimesteps"],nrecs))

    # Source time function
    lensrctf = sourcetf.size
    isrc = ijsrc[0]
    jsrc = ijsrc[1]
    
    ## Initialize arrays
    dpdx = np.zeros((nx,nz))
    dpdz = np.zeros((nx,nz))
    pcur = np.zeros((nx,nz))
    pold = np.zeros((nx,nz))
    pnew = np.zeros((nx,nz))

    ################################
    fact = vel**2 * (dt**2/dh**2)
       
    ## time loop
    print((" Time step dt: {}".format(dt)))

    for t in range(inpar["ntimesteps"]) :

        if t%100==0 :
            sys.stdout.write("\r Time step {} of {}".format(t,inpar["ntimesteps"]))
            sys.stdout.flush()

        ##=================================================
                         
        ## second order stencil
        dpdx2 = pcur[2:,1:-1]-2.0*pcur[1:-1,1:-1]+pcur[:-2,1:-1]
        dpdz2 = pcur[1:-1,2:]-2.0*pcur[1:-1,1:-1]+pcur[1:-1,:-2]

        # update pressure
        pnew[1:-1,1:-1] = 2.0*pcur[1:-1,1:-1] -pold[1:-1,1:-1] + fact[1:-1,1:-1]*(dpdx2 + dpdz2) 
             
        # inject source
        pnew[isrc,jsrc] = pnew[isrc,jsrc] + dt**2*sourcetf[t]

        ##---------------------------------------------
        ## Damp using Gaussian taper function
        pnew[:nptsgau,:]  *= gtleft.reshape(-1,1)
        pnew[-nptsgau:,:] *= gtright.reshape(-1,1)
        pnew[:,-nptsgau:] *= gtbottom.reshape(1,-1)

        pcur[:nptsgau,:]  *= gtleft.reshape(-1,1)
        pcur[-nptsgau:,:] *= gtright.reshape(-1,1)
        pcur[:,-nptsgau:] *= gtbottom.reshape(1,-1)

        if freeboundtop==False :
             pnew[:,:nptsgau]  *= gttop.reshape(1,-1)
             pcur[:,:nptsgau]  *= gttop.reshape(1,-1)
        ##----------------------------------------------

        # assign the new pold and pcur
        pold[:,:] = pcur[:,:] ## using [:,:] there is an explicit copy
        pcur[:,:] = pnew[:,:]

        ##=================================================
        
        ##### receivers
        for r in range(nrecs) :
            rec_press = _bilinear_interp(pcur,dh,recpos[r,:])
            receiv[t,r] = rec_press
        
        #### save snapshots
        if (inpar["savesnapshot"]==True) and (t%inpar["snapevery"]==0) :
            psave[:,:,tsave] = pcur
            tsave=tsave+1

    ##================================
    print(" ")
 
    if inpar["savesnapshot"]==False :
        psave = None
        
    return receiv,psave


#########################################################
####################################################

def testacou():
    
    # time
    nt = 3000
    dt = 0.0004 #s
    t = np.arange(0.0,nt*dt,dt) # s
    
    # space
    nx = 600
    nz = 400
    dh = 2.5 # m

    # source
    t0 = 0.03 # s
    f0 = 32.0 # 20.0 # Hz
    ijsrc = np.array([280,290])
 
    # source type
    from pestoseis.seismicwaves2d.sourcetimefuncs import gaussource, rickersource
    sourcetf = rickersource( t, t0, f0 )
    #sourcetf = gaussource( t, t0, f0 )

    ## velocity model
    velmod = 2000.0*np.ones((nx,nz))
    # velmod[:,:40] = 2000.0
    # velmod[:,50:200] = 2500.0
    #velmod[50:80,50:70] = 3800.0
    
    ###############################
    nrec = 8
    recpos = np.zeros((nrec,2))
    recpos[:,1] = 100.0
    recpos[:,0] = np.linspace(200.0,nx*dh-200.0,nrec)

    print(("Receiver positions:\n{}".format(recpos)))
    
    # pressure field in space and time
    inpar={}
    inpar["ntimesteps"] = nt
    inpar["nx"] = nx
    inpar["nz"] = nz
    inpar["dt"] = dt
    inpar["dh"] = dh
    inpar["savesnapshot"] = True
    inpar["snapevery"] = 10
    inpar["freesurface"]   = True
    inpar["boundcond"]= "PML" #"GaussTap"  ## "PML", "GaussTap","ReflBou"

    #--------------------------
    import time
    t1 = time.time()
    
    seism,psave = solveacoustic2D( inpar, ijsrc, velmod, sourcetf, f0, recpos )
        
    t2 = time.time()
    print(("Solver time: {}".format(t2-t1)))

    return

#########################################################
####################################################


if __name__  == "__main__" :

    testacou()
