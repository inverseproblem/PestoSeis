
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
  Functions to calculate elastic wave propagation in 2D
"""

import numpy as np
import sys
import h5py as h5

#############################################################################
#############################################################################

def bilinear_interp(f,hgrid, pt, gridshift_xy=np.array([0.0,0.0])):
    """
     Bilinear interpolation (2D).
    """
    xshift = gridshift_xy[0] 
    yshift = gridshift_xy[1]
    xreq=pt[0]
    yreq=pt[1]
    xh=(xreq-xshift)/hgrid
    yh=(yreq-yshift)/hgrid
    # assert(xh>=0.0)
    # assert(xh<=(hgrid*(f.shape[0]-1)))
    # assert(yh>=0.0)
    # assert(yh<=(hgrid*(f.shape[1]-1)))
    i=int(np.floor(xh)) # index starts from 0 so no +1
    j=int(np.floor(yh)) # index starts from 0 so no +1
    xd=xh-i
    yd=yh-j
    intval=f[i,j]*(1.0-xd)*(1.0-yd)+f[i+1,j]*(1.0-yd)*xd+f[i,j+1]*(1.0-xd)*yd+f[i+1,j+1]*xd*yd
    #print xreq,yreq,xh,yh,i,j,xd,yd    
    return intval

#############################################################################

def calc_Kab_CPML(nptspml,gridspacing,dt,Npower,d0,
                       alpha_max_pml,K_max_pml,onwhere ) :
    """
      Initialize/calculate the C-PML coefficients for boundary conditions.
    """
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

#############################################################################

def solveelastic2D(inpar, rockprops, ijsrc, sourcetf, srcdomfreq, recpos, saveh5=True,
                     outfileh5="elastic_snapshots.h5") :
    """
      Solve the elastic wave equation in 2D using finite differences on a staggered grid. Wrapper function for various boundary conditions.
       
    Args:
        inpar (dict): dictionary containing various input parameters:

                      * inpar["ntimesteps"] (int) number of time steps
                      * inpar["nx"] (int) number of grid nodes in the x direction
                      * inpar["nz"] (int) number of grid nodes in the z direction\n
                      * inpar["dt"] (float) time step for the simulation
                      * inpar["dh"] (float) grid spacing (same in x and z)
                      * inpar["savesnapshot"] (bool) switch to save snapshots of the entire wavefield
                      * inpar["snapevery"] (int) save snapshots every "snapevery" iterations
                      * inpar["freesurface"] (bool) True for free surface boundary condition at the top, False for PML\n
                      * inpar["boundcond"] (string) Type of boundary conditions "PML" or "ReflBou"
        rockprops (dict): rho (ndarray) density array
                          lambda (ndarray) lambda module array
                          mu (ndarray) shear module array
        ijsrc (ndarray(int,int)): integers representing the position of the source on the grid
        sourcetf (ndarray): source time function
        srcdomfreq (float): source dominant frequency
        recpos (ndarray): position of the receivers, a two-column array [x,z]
        saveh5 (bool): whether to save results to HDF5 file or not
        outfileh5 (string): name of the output HDF5 file

    Returns: 
        receiv (ndarray): seismograms, Vx and Vz components, recorded at the receivers
        (vxsave,vzsave) (ndarray,ndarray): set of Vx and Vz snapshots of the wavefield (if inpar["savesnapshot"]==True)

    """
    if inpar["boundcond"]=="PML" :
        receiv,vxzsave = _solveelawaveq2D_CPML(inpar, rockprops, ijsrc,
                                              sourcetf, srcdomfreq, recpos)
    elif inpar["boundcond"]=="ReflBou" :
        receiv,vxzsave = _solveelawaveq2D_ReflBound(inpar, rockprops, ijsrc,
                                                   sourcetf, srcdomfreq, recpos)
    else :
        raise ValueError("Error, wrong inpar['boundcond']: {}".format(inpar['boundcond']))

     ##############################
    if saveh5:
        ## save stuff
        hf = h5.File(outfileh5,"w")
        hf["seism"] = receiv
        if inpar["savesnapshot"]==True :
            hf["vx"] = vxzsave[0]
            hf["vz"] = vxzsave[1]
            hf["snapevery"] = inpar["snapevery"]
        hf["seismogrkind"] = inpar["seismogrkind"]
        hf["srctf"] = sourcetf
        hf["dh"] = inpar["dh"]
        hf["lambda"] = rockprops["lambda"]
        hf["mu"]  = rockprops["mu"]
        hf["rho"] = rockprops["rho"]
        hf["dt"]  = inpar["dt"]
        hf["nx"]  = inpar["nx"]
        hf["nz"]  = inpar["nz"]
        hf["recpos"] = recpos
        hf["srcij"] = ijsrc
        hf.close()
        print("Saved elastic simulation and parameters to ",outfileh5)

    return receiv,vxzsave

#############################################################################

def _solveelawaveq2D_CPML(inpar, rockprops, ijsrc, sourcetf, srcdomfreq, recpos) :
    """
    Solve the elastic wave equation in 2D using finite differences on a staggered grid. 
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
        rockprops (dict): rho (ndarray) density array
                          lambda (ndarray) lambda module array
                          mu (ndarray) shear module array
        ijsrc (ndarray(int,int)): integers representing the position of the source on the grid
        sourcetf (ndarray): source time function
        srcdomfreq (float): source dominant frequency
        recpos (ndarray): position of the receivers, a two-column array [x,z]

    Returns: 
        receiv (ndarray): seismograms, Vx and Vz components, recorded at the receivers
        (vxsave,vzsave) (ndarray,ndarray): set of Vx and Vz snapshots of the wavefield (if inpar["savesnapshot"]==True)

    """
    #  
    # Wave elastic staggered grid 2D solver 
    #
    #
    # Staggered grid with equal spacing, A. Levander (1988)
    # Second order time, fourth order space
    #
    # Convolutionary Perfectly Matched Layer (C-PML) boundary conditions
    #   
    #  See seismic_cpml/seismic_CPML_2D_isotropic_fourth_order.f90
    #
    #
    #   STAGGERED GRID 
    #
    #                      x
    #     +---------------------------------------------------->
    #     |
    #     |
    #     |                (i)    (i+1/2)   (i+1)  (i+3/2)
    #     |                 |        |        |       |
    #     |                 |        |        |       |
    #  z  |       (j) ---vx,rho-----Txx----vx,rho -----
    #     |              lam,mu     Tzz    lam,mu     |
    #     |                 |        |        |       |
    #     |                 |        |        |       |
    #     |  (j+1/2)  -----Txz------vz-------Txz-------
    #     |                 |        |        |       |
    #     |                 |        |        |       |
    #     |                 |        |        |       |
    #     |    (j+1) ----vx,rho-----Txx-----vx,rho-----
    #     |              lam,mu     Tzz     lam,mu
    #     v
    #      
    #   Where
    #
    #   Txx Stress_xx (Normal stress x)
    #   Tzz: Stress_zz (Normal stress z)
    #   Txz: Stress_xz (Shear stress)
    #   vx: Velocity_x (x component of velocity) 
    #   vz: Velocity_z (z component of velocity) 
    #
    #
    # Node indexing:
    # ------------------------------------------------
    # |      Code         |   Staggered grid         |
    # ------------------------------------------------
    # | rho(i,j)          | rho(i,j)                 |
    # | lam(i,j),mu(i,j)  | lam(i,j),mu(i,j)         |
    # | Txx(i,j),Tzz(i,j) | Txx(i+1/2,j),Tzz(i+1/2,j)|
    # | Txz(i,j)          | Txz(i,j+1/2)             |
    # | vx(i,j)           | vx(i,j)                  | 
    # | vz(i,j)           | vz(i+1/2,j+1/2)          |
    # ------------------------------------------------
        
    assert(inpar["boundcond"]=="PML")
    print("Starting ELASTIC solver with CPML boundary conditions.")


    ##############################  
    # Lame' parameters
    ##############################
    lamb = rockprops["lambda"]
    mu = rockprops["mu"]
    ##############################
    # Density
    ##############################
    rho = rockprops["rho"]

    #############################
    dh = inpar["dh"]
    dx = dh
    dz = dh
    dt = inpar["dt"]
    if inpar["sourcetype"] == "MomTensor" :
        MTens = inpar['momtensor']
    elif  inpar["sourcetype"] == "ExtForce" :
        ExtForce = inpar['extforce']
    print(" Source type: ",inpar["sourcetype"])

    ##############################
    ## Check stability criterion
    ##############################
    maxvp = ( np.sqrt( (lamb+2.0*mu)/rho )).max()
    Courant_number = maxvp * inpar["dt"] * np.sqrt(1.0/dx**2 + 1.0/dz**2)
    if Courant_number > 1.0 :
        print(" The stability criterion is violated. Quitting.")
        return
    
    #minvp = minimum( sqrt( (lambda+2.0*mu)./rho ))
   
    ##############################
    #   Parameters
    ##############################
    harmonicaver_mu = False
    
    ## pml = c-pml, else reflsurf
    ## freesur = free surf for top
    if inpar["freesurface"]==True :
        print(" * Free surface at the top *")
        freeboundtop=True
    else :
        print(" * Absorbing boundary at the top *")
        freeboundtop=False

    f0 = srcdomfreq # source
    
    ## keep it simple for now...
    nx = inpar["nx"]
    nz = inpar["nz"]
    
    ##############################
    #   Parameters for PML
    ##############################
    nptspml_x = 21 ## 21 ref coeff set for 21 absorbing nodes
    nptspml_z = 21 ## 21
    assert(nptspml_x<ijsrc[0]<(nx-1-nptspml_x))
    assert(ijsrc[1]<(nz-1-nptspml_z))
    if freeboundtop==False:
        assert(nptspml_z<ijsrc[1])
        
    print(" Size of PML layers in grid points: {} in x and {} in z".format(nptspml_x,nptspml_z))
    
    Npower = 2.0
    assert Npower >= 1
    K_max_pml = 1.0
    alpha_max_pml = 2.0*np.pi*(f0/2.0) 
    
    # reflection coefficient (INRIA report section 6.1)
    # http://hal.inria.fr/docs/00/07/32/19/PDF/RR-3471.pdf
    Rcoef = 0.0001
       
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
    K_xpml,a_xpml,b_xpml = calc_Kab_CPML(nptspml_x,dh,dt,Npower,d0_x,alpha_max_pml,K_max_pml,"grdpts")
    # damping profile at half the grid points
    K_xpml_half,a_xpml_half,b_xpml_half = calc_Kab_CPML(nptspml_x,dh,dt,Npower,d0_x,alpha_max_pml,K_max_pml,"halfgrdpts")
    
    # --- damping in the z direction ---
    # assuming the number of grid points for PML is the same on
    # both sides    
    # damping profile at the grid points
    K_zpml,a_zpml,b_zpml = calc_Kab_CPML(nptspml_z,dh,dt,Npower,d0_z,alpha_max_pml,K_max_pml,"grdpts")
    # damping profile at half the grid points
    K_zpml_half,a_zpml_half,b_zpml_half = calc_Kab_CPML(nptspml_z,dh,dt,Npower,d0_z,alpha_max_pml,K_max_pml,"halfgrdpts")
    

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
        vxsave = np.zeros((inpar["nx"],inpar["nz"],ntsave+1))
        vzsave = np.zeros((inpar["nx"],inpar["nz"],ntsave+1))
        tsave=1


    ## Arrays to return seismograms
    nrecs = recpos.shape[0]
    receiv = np.zeros((nrecs,inpar["ntimesteps"],2))


    # Source time function
    lensrctf = sourcetf.size
    isrc = ijsrc[0]
    jsrc = ijsrc[1]

    
    ## Initialize arrays
    vx = np.zeros((nx,nz))
    vz = np.zeros((nx,nz))
    Txx = np.zeros((nx,nz))
    Tzz = np.zeros((nx,nz))
    Txz = np.zeros((nx,nz))

    ## derivatives
    Dxb_Txx = np.zeros((nx-3,nz-3))
    Dzb_Txz = np.zeros((nx-3,nz-3))
    Dxf_Txz = np.zeros((nx-3,nz-3))
    Dzf_Tzz = np.zeros((nx-3,nz-3))
    Dxf_vx  = np.zeros((nx-3,nz-3))
    Dzb_vz  = np.zeros((nx-3,nz-3))
    Dzf_vx  = np.zeros((nx-3,nz-3))
    Dxb_vz  = np.zeros((nx-3,nz-3))


    # PML arrays
    # Arrays with size of PML areas would be sufficient and save memory,
    #   however allocating arrays with same size than model simplifies
    #   the code in the loops
    psi_DxTxx = np.zeros((nx,nz))
    psi_DzTxz = np.zeros((nx,nz))
    
    psi_DxTxz = np.zeros((nx,nz))
    psi_DzTzz = np.zeros((nx,nz))
    
    psi_DxVx = np.zeros((nx,nz))
    psi_DzVz = np.zeros((nx,nz))

    psi_DzVx = np.zeros((nx,nz))
    psi_DxVz = np.zeros((nx,nz))

    
    ##############################
    #   Derivative operators
    ##############################
    #
    #  o => point where to take the derivative
    #
    # forward operator: 1/24 * (f[i-1]-27*f[i]+27*f[i+1]-f[i+2])
    #
    #  i-1   i    i+1  i+2
    #   |    |  o  |    |
    #
    # backward operator: 1/24 * (f[i-2]-27*f[i-1]+27*f[i]-f[i+1])
    #
    #  i-2  i-1    i   i+1
    #   |    |  o  |    |
    #    
    # Weigths for taking derivarives for incresing indices
    #Dweights = 1.0/inpar["dh"] * [1/24.0, -27.0/24.0, 27.0/24.0, -1/24.0]

    ###############################################################
    # pre-interpolate properties at half distances between nodes
    ###############################################################
    # rho_ihalf_jhalf (nx-1,ny-1) ??
    rho_ihalf_jhalf = (rho[1:,1:]+rho[1:,:-1]+rho[:-1,1:]+rho[:-1,:-1])/4.0
    # mu_ihalf (nx-1,ny) ??
    # mu_jhalf (nx,ny-1) ??
    if harmonicaver_mu==True :
        mu_ihalf = 2.0/( 1.0/mu[1:,:] + 1.0/mu[:-1,:] )
        mu_jhalf = 2.0/( 1.0/mu[:,1:] + 1.0/mu[:,:-1] )
    else :
        mu_ihalf = (mu[1:,:]+mu[:-1,:])/2.0 ###?????
        mu_jhalf = (mu[:,1:]+mu[:,:-1])/2.0 ###?????
    # lamb_ihalf (nx-1,ny) ??
    lamb_ihalf = (lamb[1:,:]+lamb[:-1,:])/2.0 ###?????
    
    ##======================================##

    fact = 1.0/(24.0*inpar["dh"])
    gridshift_vz = np.array([dh/2.0,dh/2.0])
    
    ## to make the Numpy broadcast Hadamard product correct along x
    a_x = a_x.reshape(-1,1)
    b_x = b_x.reshape(-1,1)
    K_x = K_x.reshape(-1,1)
    a_x_half = a_x_half.reshape(-1,1)
    b_x_half = b_x_half.reshape(-1,1)
    K_x_half = K_x_half.reshape(-1,1)
     
    ## time loop
    dt = inpar["dt"]
    print(" Time step dt: {}".format(dt))
    for t in range(inpar["ntimesteps"]) :

        if t%100==0 :
            sys.stdout.write("\r Time step {} of {}".format(t,inpar["ntimesteps"]))
            sys.stdout.flush()
            #print("\r Time step {} of {}".format(t,inpar["ntimesteps"]))

        
        ## Inject the source 
        if t<=lensrctf :
            
            if inpar["sourcetype"]=="MomTensor":
                Txx[isrc,jsrc] = Txx[isrc,jsrc] + MTens['xx'] * sourcetf[t]* dt / dh**2
                Tzz[isrc,jsrc] = Tzz[isrc,jsrc] + MTens['zz'] * sourcetf[t]* dt / dh**2
                Txz[isrc,jsrc] = Txz[isrc,jsrc] + MTens['xz'] * sourcetf[t]* dt / dh**2

            elif inpar["sourcetype"]=="ExtForce":
                # rho_half_x_half_y = 0.25d0 * (rho(i,j) + rho(i+1,j) + rho(i+1,j+1) + rho(i,j+1))
                # vx(i,j) = vx(i,j) + force_x * DELTAT / rho(i,j)
                # vy(i,j) = vy(i,j) + force_y * DELTAT / rho_half_x_half_y
                
                vx[isrc,jsrc] = vx[isrc,jsrc] + ExtForce['x'] * dt / rho[isrc,jsrc]
                rho_half_x_half_y = 0.25 * (rho[isrc,jsrc] + rho[isrc+1,jsrc] + rho[isrc+1,jsrc+1] + rho[isrc,jsrc+1])
                vz[isrc,jsrc] = vz[isrc,jsrc] + ExtForce['z'] * dt / rho_ihalf_jhalf[isrc,jsrc]

            else :
                print("Error source type badly defined.")
                return
                
        ## space loops excluding boundaries
        
        #########################################
        # update velocities from stresses
        #########################################

        if freeboundtop==True :
            ## j=0,1 

            ### Vx
            Dxb_Txx = fact * ( Txx[:-3,:2] -27.0*Txx[1:-2,:2] +27.0*Txx[2:-1,:2] -Txx[3:,:2] )
            # Dzb_Txz = fact * ( -Txz[2:-1,1:3] +27.0*Txz[2:-1,:2] +27.0*Txz[2:-1,:2] -Txz[2:-1,1:3] )
            Dzb_Txz = np.zeros_like(Txz[2:-1,:2])
            Dzb_Txz[:, 0] = fact * ( -Txz[2:-1,1] +27.0*Txz[2:-1,0] +27.0*Txz[2:-1,0] -Txz[2:-1,1] )
            Dzb_Txz[:, 1] = fact * ( -Txz[2:-1,0] -27.0*Txz[2:-1,0] +27.0*Txz[2:-1,1] -Txz[2:-1,2] )
            # vx
            vx[2:-1,:2] = vx[2:-1,:2] + (dt/rho[2:-1,:2]) * (Dxb_Txx + Dzb_Txz)

            ###---------------------------------------------------
            # Vz
            Dxf_Txz = fact * ( Txz[:-3,:2] -27.0*Txz[1:-2,:2] +27.0*Txz[2:-1,:2] -Txz[3:,:2] )
            # Dzf_Tzz = fact * ( -Tzz[1:-2,2:4] +27.0*Tzz[1:-2,1:3] +27.0*Tzz[1:-2,1:3] -Tzz[1:-2,2:4] )
            Dzf_Tzz = np.zeros_like(Tzz[1:-2,2:4])
            Dzf_Tzz[:, 0] = fact * ( -Tzz[1:-2,2] +27.0*Tzz[1:-2,1] +27.0*Tzz[1:-2,1] -Tzz[1:-2,2] )
            Dzf_Tzz[:, 1] = fact * ( -Tzz[1:-2,1] -27.0*Tzz[1:-2,1] +27.0*Tzz[1:-2,2] -Tzz[1:-2,3] )
    
            # update velocity (rho has been interpolated in advance)
            # rho_ihalf_jhalf[1:-1,1:-1] because its size is (nx-1,ny-1) 
            vz[1:-2,:2] = vz[1:-2,:2] + (dt/rho_ihalf_jhalf[1:-1,:2]) * (Dxf_Txz + Dzf_Tzz)
            
            ## END free surface
            ##============================================

        #-----------------------------------------------------
        ### Vx
        Dxb_Txx = fact * ( Txx[:-3,2:-1] -27.0*Txx[1:-2,2:-1] +27.0*Txx[2:-1,2:-1] -Txx[3:,2:-1] )
        Dzb_Txz = fact * ( Txz[2:-1,:-3] -27.0*Txz[2:-1,1:-2] +27.0*Txz[2:-1,2:-1] -Txz[2:-1,3:] )

        # C-PML stuff 
        psi_DxTxx[2:-1,2:-1] = b_x[2:-1] * psi_DxTxx[2:-1,2:-1] + a_x[2:-1] *Dxb_Txx
        psi_DzTxz[2:-1,2:-1] = b_z[2:-1] * psi_DzTxz[2:-1,2:-1] + a_z[2:-1] *Dzb_Txz
        
        Dxb_Txx = Dxb_Txx / K_x[2:-1] + psi_DxTxx[2:-1,2:-1]
        Dzb_Txz = Dzb_Txz / K_z[2:-1] + psi_DzTxz[2:-1,2:-1]
                
        # update velocity
        vx[2:-1,2:-1] = vx[2:-1,2:-1] + (dt/rho[2:-1,2:-1]) * (Dxb_Txx + Dzb_Txz)
    
        
        #-----------------------------------------------------
        # Vz
        Dxf_Txz = fact * ( Txz[:-3,1:-2] -27.0*Txz[1:-2,1:-2] +27.0*Txz[2:-1,1:-2] -Txz[3:,1:-2] )
        Dzf_Tzz = fact * ( Tzz[1:-2,:-3] -27.0*Tzz[1:-2,1:-2] +27.0*Tzz[1:-2,2:-1] -Tzz[1:-2,3:] )
    
        # C-PML stuff 
        psi_DxTxz[1:-2,1:-2] = b_x_half[1:-2] * psi_DxTxz[1:-2,1:-2] + a_x_half[1:-2]*Dxf_Txz
        psi_DzTzz[1:-2,1:-2] = b_z_half[1:-2] * psi_DzTzz[1:-2,1:-2] + a_z_half[1:-2]*Dzf_Tzz
        
        Dxf_Txz = Dxf_Txz / K_x_half[1:-2] + psi_DxTxz[1:-2,1:-2]
        Dzf_Tzz = Dzf_Tzz / K_z_half[1:-2] + psi_DzTzz[1:-2,1:-2]
        
        # update velocity (rho has been interpolated in advance)
        # rho_ihalf_jhalf[1:-1,1:-1] because its size is (nx-1,ny-1) 
        vz[1:-2,1:-2] = vz[1:-2,1:-2] + (dt/rho_ihalf_jhalf[1:-1,1:-1]) * (Dxf_Txz + Dzf_Tzz)
  
        
        #########################################
        # update stresses from velocities 
        #########################################
        
        if freeboundtop==True :
            ##--------------------------------------------
            ## j=0

            # Txx,Tzz
            Dxf_vx = fact * (vx[:-3,0] -27.0*vx[1:-2,0] +27.0*vx[2:-1,0] -vx[3:,0])
            Dzb_vz = -(1.0-2.0*mu_ihalf[1:-1,0]/lamb_ihalf[1:-1,0])*Dxf_vx
            # Txx
            # Txx[1:-2,0] = Txx[1:-2,0] + (lamb_ihalf[1:-1,0]+2.0*mu_ihalf[1:-1,0]) * dt * Dxf_vx + lamb_ihalf[1:-1,0] * dt * Dzb_vz
            Txx[1:-2,0] = Txx[1:-2,0] + (lamb_ihalf[1:-1,0] - lamb_ihalf[1:-1,0] / (lamb_ihalf[1:-1,0] + 2 * mu_ihalf[1:-1,0]) + 2 * mu_ihalf[1:-1,0]) * dt * Dxf_vx
            # Txx[1:-2,0] = Txx[1:-2,0] + dt * (4 * mu_ihalf[1:-1,0] - (2 * mu_ihalf[1:-1,0]) ** 2 / (lamb_ihalf[1:-1,0] + 2 * mu_ihalf[1:-1,0])) * Dxf_vx
            # Tzz
            Tzz[1:-2,0] = 0.0

            ##----------------------------------------------
            ## j=1
            # Txx,Tzz
            Dxf_vx = fact * (vx[:-3,1] -27.0*vx[1:-2,1] +27.0*vx[2:-1,1] -vx[3:,1])
            Dzb_vz = fact * ( 0.0 -27.0*vz[1:-2,0] +27.0*vz[1:-2,1] -vz[1:-2,2] )
            
            # Txx
            Txx[1:-2,1] = Txx[1:-2,1] + (lamb_ihalf[1:-1,1]+2.0*mu_ihalf[1:-1,1]) * dt * Dxf_vx + lamb_ihalf[1:-1,1] * dt * Dzb_vz
            # Tzz
            Tzz[1:-2,1] = Tzz[1:-2,1] + (lamb_ihalf[1:-1,1]+2.0*mu_ihalf[1:-1,1]) * dt* Dzb_vz + lamb_ihalf[1:-1,1] * dt * Dxf_vx


            ##===============================================
            ## j=0
            # Txz
            Dzf_vx = fact * (0.0 -27.0*vx[2:-1,0] +27.0*vx[2:-1,1] -vx[2:-1,2])
            Dxb_vz = fact * (vz[:-3,0] -27.0*vz[1:-2,0] +27.0*vz[2:-1,0] -vz[3:,0]) 
            # Txz
            Txz[2:-1,0] = Txz[2:-1,0] + mu_jhalf[2:-1,0] * dt * (Dzf_vx + Dxb_vz)
            # Txz[2:-1,0] = 0.0
            ##--------------------------------------------

            ## END free surface
            ##============================================
            
      
        #-----------------------------------------------------
        # Txx,Tzz
        Dxf_vx = fact * (vx[:-3,2:-1] -27.0*vx[1:-2,2:-1] +27.0*vx[2:-1,2:-1] -vx[3:,2:-1])
        Dzb_vz = fact * (vz[1:-2,:-3] -27.0*vz[1:-2,1:-2] +27.0*vz[1:-2,2:-1] -vz[1:-2,3:])
                
        # C-PML stuff 
        psi_DxVx[1:-2,2:-1] = b_x_half[1:-2] * psi_DxVx[1:-2,2:-1] + a_x_half[1:-2]*Dxf_vx
        psi_DzVz[1:-2,2:-1] = b_z[2:-1] * psi_DzVz[1:-2,2:-1] + a_z[2:-1]*Dzb_vz

        Dxf_vx = Dxf_vx / K_x_half[1:-2] + psi_DxVx[1:-2,2:-1]
        Dzb_vz = Dzb_vz / K_z[2:-1] + psi_DzVz[1:-2,2:-1]
        
        # Txx
        # ?_ihalf[1:-1,2:-1], etc. because its size is (nx-1,ny-1)
        #print Txx[1:-2,2:-1].shape, lamb_ihalf[1:-1,2:-1].shape
        Txx[1:-2,2:-1] = Txx[1:-2,2:-1] + (lamb_ihalf[1:-1,2:-1]+2.0*mu_ihalf[1:-1,2:-1]) * dt * Dxf_vx + lamb_ihalf[1:-1,2:-1] * dt * Dzb_vz
        
        ## derivatives are the same than for Txx
        # Tzz
        Tzz[1:-2,2:-1] = Tzz[1:-2,2:-1] + (lamb_ihalf[1:-1,2:-1]+2.0*mu_ihalf[1:-1,2:-1]) * dt* Dzb_vz + lamb_ihalf[1:-1,2:-1] * dt * Dxf_vx

        #-----------------------------------------------------
        # Txz
        Dzf_vx = fact * (vx[2:-1,:-3] -27.0*vx[2:-1,1:-2] +27.0*vx[2:-1,2:-1] -vx[2:-1,3:])
        Dxb_vz = fact * (vz[:-3,1:-2] -27.0*vz[1:-2,1:-2] +27.0*vz[2:-1,1:-2] -vz[3:,1:-2]) 
                      
        # C-PML stuff
        psi_DzVx[2:-1,1:-2] = b_z_half[1:-2] * psi_DzVx[2:-1,1:-2] + a_z_half[1:-2]*Dzf_vx
        psi_DxVz[2:-1,1:-2] = b_x[2:-1] * psi_DxVz[2:-1,1:-2] + a_x[2:-1]*Dxb_vz

        Dzf_vx = Dzf_vx / K_z[1:-2] + psi_DzVx[2:-1,1:-2]
        Dxb_vz = Dxb_vz / K_x[2:-1] + psi_DxVz[2:-1,1:-2]
        
        # Txz
        Txz[2:-1,1:-2] = Txz[2:-1,1:-2] + mu_jhalf[2:-1,1:-1] * dt * (Dzf_vx + Dxb_vz)
                         

        ##------------------------------------------------
        ##### receivers
        for r in range(nrecs) :
            rec_vx = bilinear_interp(vx,dh,recpos[r,:])
            rec_vz = bilinear_interp(vz,dh,recpos[r,:],gridshift_xy=gridshift_vz)
            receiv[r,t,:] = np.array([rec_vx, rec_vz])
        
        #### save snapshots
        if (inpar["savesnapshot"]==True) and (t%inpar["snapevery"]==0) :
            vxsave[:,:,tsave] = vx
            vzsave[:,:,tsave] = vz
            tsave=tsave+1
    
    ##### End time loop ################
    print(" ")
 
    if inpar["savesnapshot"]==False :
        psave = None
        
    if (inpar["savesnapshot"]==True):
        return receiv,(vxsave,vzsave)
    else:
        return receiv, ()


#############################################################################

def _solveelawaveq2D_ReflBound(inpar, rockprops, ijsrc, sourcetf, srcdomfreq, recpos) :
    """
    Solve the elastic wave equation in 2D using finite differences on a staggered grid. 
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
        rockprops (dict): rho (ndarray) density array
                          lambda (ndarray) lambda module array
                          mu (ndarray) shear module array
        ijsrc (ndarray(int,int)): integers representing the position of the source on the grid
        sourcetf (ndarray): source time function
        srcdomfreq (float): source dominant frequency
        recpos (ndarray): position of the receivers, a two-column array [x,z]

    Returns: 
        receiv (ndarray): seismograms, Vx and Vz components, recorded at the receivers
        (vxsave,vzsave) (ndarray,ndarray): set of Vx and Vz snapshots of the wavefield (if inpar["savesnapshot"]==True)

    """
    #  
    # Wave elastic staggered grid 2D solver 
    #
    #
    # Staggered grid with equal spacing, A. Levander (1988)
    # Second order time, fourth order space
    #
    # 
    #
    #   STAGGERED GRID 
    #
    #                      x
    #     +---------------------------------------------------->
    #     |
    #     |
    #     |                (i)    (i+1/2)   (i+1)  (i+3/2)
    #     |                 |        |        |       |
    #     |                 |        |        |       |
    #  z  |       (j) ---vx,rho-----Txx----vx,rho -----
    #     |              lam,mu     Tzz    lam,mu     |
    #     |                 |        |        |       |
    #     |                 |        |        |       |
    #     |  (j+1/2)  -----Txz------vz-------Txz-------
    #     |                 |        |        |       |
    #     |                 |        |        |       |
    #     |                 |        |        |       |
    #     |    (j+1) ----vx,rho-----Txx-----vx,rho-----
    #     |              lam,mu     Tzz     lam,mu
    #     v
    #      
    #   Where
    #
    #   Txx Stress_xx (Normal stress x)
    #   Tzz: Stress_zz (Normal stress z)
    #   Txz: Stress_xz (Shear stress)
    #   vx: Velocity_x (x component of velocity) 
    #   vz: Velocity_z (z component of velocity) 
    #
    #
    # Node indexing:
    # ------------------------------------------------
    # |      Code         |   Staggered grid         |
    # ------------------------------------------------
    # | rho(i,j)          | rho(i,j)                 |
    # | lam(i,j),mu(i,j)  | lam(i,j),mu(i,j)         |
    # | Txx(i,j),Tzz(i,j) | Txx(i+1/2,j),Tzz(i+1/2,j)|
    # | Txz(i,j)          | Txz(i,j+1/2)             |
    # | vx(i,j)           | vx(i,j)                  | 
    # | vz(i,j)           | vz(i+1/2,j+1/2)          |
    # ------------------------------------------------
        
    assert(inpar["boundcond"]=="ReflBou")
    print("Starting ELASTIC solver with reflective boundary condition all around.")

    ##############################  
    # Lame' parameters
    ##############################
    lamb = rockprops["lambda"]
    mu = rockprops["mu"]
    ##############################
    # Density
    ##############################
    rho = rockprops["rho"]

    #############################
    dh = inpar["dh"]
    dx = dh
    dz = dh
    dt = inpar["dt"]
    if inpar["sourcetype"] == "MomTensor" :
        MTens = inpar['momtensor']
    elif  inpar["sourcetype"] == "ExtForce" :
        ExtForce = inpar['extforce']
    print(" Source type: ",inpar["sourcetype"])
    
    ##############################
    ## Check stability criterion
    ##############################
    maxvp = ( np.sqrt( (lamb+2.0*mu)/rho )).max()
    Courant_number = maxvp * inpar["dt"] * np.sqrt(1.0/dx**2 + 1.0/dz**2)
    if Courant_number > 1.0 :
        print(" The stability criterion is violated. Quitting.")
        return
    
    #minvp = minimum( sqrt( (lambda+2.0*mu)./rho ))
    

    ##############################
    #   Parameters
    ##############################

    harmonicaver_mu = False
    f0 = srcdomfreq # source
    
    ## keep it simple for now...
    nx = inpar["nx"]
    nz = inpar["nz"]
        
    #####################################################

    ## Arrays to export snapshots
    if inpar["savesnapshot"]==True :
        ntsave = inpar["ntimesteps"]//inpar["snapevery"]
        vxsave = np.zeros((inpar["nx"],inpar["nz"],ntsave+1))
        vzsave = np.zeros((inpar["nx"],inpar["nz"],ntsave+1))
        tsave=1

    ## Arrays to return seismograms
    nrecs = recpos.shape[0]
    receiv = np.zeros((nrecs,inpar["ntimesteps"],2))

    # Source time function
    lensrctf = sourcetf.size
    isrc = ijsrc[0]
    jsrc = ijsrc[1]
    
    ## Initialize arrays
    vx = np.zeros((nx,nz))
    vz = np.zeros((nx,nz))
    Txx = np.zeros((nx,nz))
    Tzz = np.zeros((nx,nz))
    Txz = np.zeros((nx,nz))

    ## derivatives
    Dxb_Txx = np.zeros((nx-3,nz-3))
    Dzb_Txz = np.zeros((nx-3,nz-3))
    Dxf_Txz = np.zeros((nx-3,nz-3))
    Dzf_Tzz = np.zeros((nx-3,nz-3))
    Dxf_vx  = np.zeros((nx-3,nz-3))
    Dzb_vz  = np.zeros((nx-3,nz-3))
    Dzf_vx  = np.zeros((nx-3,nz-3))
    Dxb_vz  = np.zeros((nx-3,nz-3))
    
    ##############################
    #   Derivative operators
    ##############################
    #
    #  o => point where to take the derivative
    #
    # forward operator: 1/24 * (f[i-1]-27*f[i]+27*f[i+1]-f[i+2])
    #
    #  i-1   i    i+1  i+2
    #   |    |  o  |    |
    #
    # backward operator: 1/24 * (f[i-2]-27*f[i-1]+27*f[i]-f[i+1])
    #
    #  i-2  i-1    i   i+1
    #   |    |  o  |    |
    #    

    ###############################################################
    # pre-interpolate properties at half distances between nodes
    ###############################################################
    # rho_ihalf_jhalf (nx-1,ny-1) ??
    rho_ihalf_jhalf = (rho[1:,1:]+rho[1:,:-1]+rho[:-1,1:]+rho[:-1,:-1])/4.0
    # mu_ihalf (nx-1,ny) ??
    # mu_jhalf (nx,ny-1) ??
    if harmonicaver_mu==True :
        mu_ihalf = 1.0/( 1.0/mu[1:,:] + 1.0/mu[:-1,:] )
        mu_jhalf = 1.0/( 1.0/mu[:,1:] + 1.0/mu[:,:-1] )
    else :
        mu_ihalf = (mu[1:,:]+mu[:-1,:])/2.0 ###?????
        mu_jhalf = (mu[:,1:]+mu[:,:-1])/2.0 ###?????
    # lamb_ihalf (nx-1,ny) ??
    lamb_ihalf = (lamb[1:,:]+lamb[:-1,:])/2.0 ###?????
    
    ##======================================##
    ##======================================##
    
    fact = 1.0/(24.0*inpar["dh"])
    gridshift_vz = np.array([dh/2.0,dh/2.0])

    ## time loop
    dt = inpar["dt"]
    print(" Time step dt: {}".format(dt))
    for t in range(inpar["ntimesteps"]) :

        if t%100==0 :
            sys.stdout.write("\r Time step {} of {}".format(t,inpar["ntimesteps"]))
            sys.stdout.flush()
            #print("\r Time step {} of {}".format(t,inpar["ntimesteps"]))

        
        ## Inject the source 
        if t<=lensrctf :
            
            if inpar["sourcetype"]=="MomTensor":
                Txx[isrc,jsrc] = Txx[isrc,jsrc] + MTens['xx'] * sourcetf[t]* dt / dh**2
                Tzz[isrc,jsrc] = Tzz[isrc,jsrc] + MTens['zz'] * sourcetf[t]* dt / dh**2
                Txz[isrc,jsrc] = Txz[isrc,jsrc] + MTens['xz'] * sourcetf[t]* dt / dh**2

            elif inpar["sourcetype"]=="ExtForce":
                # ExtForce['x'] = # dt/rho *
                # Extforce['z'] =  
                raise ValueError("ExtForce source not yet implemented! .Exiting.")

            else :
                print("Error source type badly defined.")
                return
                
        ## space loops excluding boundaries
        
        #########################################
        # update velocities from stresses
        #########################################

        #-----------------------------------------------------
        ### Vx
        Dxb_Txx = fact * ( Txx[:-3,2:-1] -27.0*Txx[1:-2,2:-1] +27.0*Txx[2:-1,2:-1] -Txx[3:,2:-1] )
        Dzb_Txz = fact * ( Txz[2:-1,:-3] -27.0*Txz[2:-1,1:-2] +27.0*Txz[2:-1,2:-1] -Txz[2:-1,3:] )

        # update velocity
        vx[2:-1,2:-1] = vx[2:-1,2:-1] + (dt/rho[2:-1,2:-1]) * (Dxb_Txx + Dzb_Txz)
    
        
        #-----------------------------------------------------
        # Vz
        Dxf_Txz = fact * ( Txz[:-3,1:-2] -27.0*Txz[1:-2,1:-2] +27.0*Txz[2:-1,1:-2] -Txz[3:,1:-2] )
        Dzf_Tzz = fact * ( Tzz[1:-2,:-3] -27.0*Tzz[1:-2,1:-2] +27.0*Tzz[1:-2,2:-1] -Tzz[1:-2,3:] )
    
        # update velocity (rho has been interpolated in advance)
        # rho_ihalf_jhalf[1:-1,1:-1] because its size is (nx-1,ny-1) 
        vz[1:-2,1:-2] = vz[1:-2,1:-2] + (dt/rho_ihalf_jhalf[1:-1,1:-1]) * (Dxf_Txz + Dzf_Tzz)
  
        
        #########################################
        # update stresses from velocities 
        #########################################
              
        #-----------------------------------------------------
        # Txx,Tzz
        Dxf_vx = fact * (vx[:-3,2:-1] -27.0*vx[1:-2,2:-1] +27.0*vx[2:-1,2:-1] -vx[3:,2:-1])
        Dzb_vz = fact * (vz[1:-2,:-3] -27.0*vz[1:-2,1:-2] +27.0*vz[1:-2,2:-1] -vz[1:-2,3:])
                
        # Txx
        # ?_ihalf[1:-1,2:-1], etc. because its size is (nx-1,ny-1)
        #print Txx[1:-2,2:-1].shape, lamb_ihalf[1:-1,2:-1].shape
        Txx[1:-2,2:-1] = Txx[1:-2,2:-1] + (lamb_ihalf[1:-1,2:-1]+2.0*mu_ihalf[1:-1,2:-1]) * dt * Dxf_vx + lamb_ihalf[1:-1,2:-1] * dt * Dzb_vz
        
        ## derivatives are the same than for Txx
        # Tzz
        Tzz[1:-2,2:-1] = Tzz[1:-2,2:-1] + (lamb_ihalf[1:-1,2:-1]+2.0*mu_ihalf[1:-1,2:-1]) * dt* Dzb_vz + lamb_ihalf[1:-1,2:-1] * dt * Dxf_vx

        #-----------------------------------------------------
        # Txz
        Dzf_vx = fact * (vx[2:-1,:-3] -27.0*vx[2:-1,1:-2] +27.0*vx[2:-1,2:-1] -vx[2:-1,3:])
        Dxb_vz = fact * (vz[:-3,1:-2] -27.0*vz[1:-2,1:-2] +27.0*vz[2:-1,1:-2] -vz[3:,1:-2]) 
                      
        # Txz
        Txz[2:-1,1:-2] = Txz[2:-1,1:-2] + mu_jhalf[2:-1,1:-1] * dt * (Dzf_vx + Dxb_vz)
                         

        ##=============================================================================
        
        ##### receivers
        for r in range(nrecs) :
            rec_vx = bilinear_interp(vx,dh,recpos[r,:])
            rec_vz = bilinear_interp(vz,dh,recpos[r,:],gridshift_xy=gridshift_vz)
            receiv[r,t,:] = np.array([rec_vx, rec_vz])
        
        #### save snapshots
        if (inpar["savesnapshot"]==True) and (t%inpar["snapevery"]==0) :
            vxsave[:,:,tsave] = vx
            vzsave[:,:,tsave] = vz
            tsave=tsave+1
    
    ##### End time loop ################
    ##--------------------------
    print(" ")

    if inpar["savesnapshot"]==False :
        psave = None
        
    return receiv,(vxsave,vzsave)
  

#############################################################################

def _testela() :
    
    # time
    nt = 4000
    dt = 0.0002 #s
    t = np.arange(0.0,nt*dt,dt) # s
    
    # space
    nx = 250
    nz = 250
    dh = 5.0 # m

    # source
    t0 = 0.03 # s
    f0 = 32.0 # 20.0 # Hz
    ijsrc = np.array([109,109])

    ###############################3
    from sourcetimefuncs import gaussource, rickersource
    sourcetf = rickersource( t, t0, f0 )
    #sourcetf = gaussource1D( t, t0, f0 )
    MTens = dict(xx=1.0, xz=0.0, zz=1.0)
    ExtForce= dict(x=1.0, z=3.5)
    
    # srctype = "ExtForce"
    # ExtForce= dict(x= , z= )


    ## lambda,mu,rho
    rho = 2700.0*np.ones((nx,nz))
    rho[:,nx//2:] = 2950.0
    rho[nx//2:,nz//2:nz//2+20] = 2350.0
    
    P_wav = 3900.0*np.ones((nx,nz))
    P_wav[:,nz//2:] = 4500.0
    P_wav[nx//2:,nz//2:nz//2+20] = 4100.0

    S_wav = 2500.0*np.ones((nx,nz))
    S_wav[:,nz//2:] = 2700.0
    S_wav[nx//2:,nz//2:nz//2+20] = 2100.0

    lamb = rho*(P_wav**2-2.0*S_wav**2) # P-wave modulus
    mu = rho*S_wav**2 # S-wave modulus (shear modulus)
    

    ###############################
    nrec = 4
    recpos = np.zeros((nrec,2))
    recpos[:,1] = 600.0
    recpos[:,0] = np.linspace(200.0,nx*dh-200.0,nrec)
    #print("Receiver positions: \n{}".format(recpos))

    
    # pressure field in space and time
    inpar={}
    inpar["ntimesteps"] = nt
    inpar["nx"] = nx
    inpar["nz"] = nz
    inpar["dt"] = dt
    inpar["dh"] = dh
    inpar["sourcetype"] = "MomTensor" # "MomTensor" "ExtForce"
    inpar["momtensor"] = MTens
    inpar["extforce"]  = ExtForce 
    inpar["savesnapshot"] = True
    inpar["snapevery"] = 50
    inpar["seismogrkind"] = "velocity"    
    inpar["freesurface"]  = True
    inpar["boundcond"] = "PML"  #"PML" "ReflBou"

    #--------------------------
    rockprops = {}
    rockprops["lambda"] = lamb
    rockprops["mu"] = mu
    rockprops["rho"] = rho

    import time
    t1 = time.time()

    seism,vxzsave = solveelastic2D(inpar, rockprops,ijsrc, sourcetf,f0, recpos)
        
    t2 = time.time()
    print("Solver time: {}".format(t2-t1))
        
    return

#############################################################################

if __name__  == "__main__" :

    _testela()
