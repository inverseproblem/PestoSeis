
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

# -*- coding: utf-8 -*

#
# Andrea Zunino
#  zunino[at]nbi[dot]dk
# 
#

import numpy as NP
import sys

#############################################################################
#############################################################################

def bilinear_interp(f,hgrid, pt):
    xreq=pt[0]
    yreq=pt[1]
    xh=xreq/hgrid
    yh=yreq/hgrid
    i=int(NP.floor(xh)) # index starts from 0 so no +1
    j=int(NP.floor(yh)) # index starts from 0 so no +1
    xd=xh-i # index starts from 0 so no +1
    yd=yh-j # index starts from 0 so no +1
    intval=f[i,j]*(1.0-xd)*(1.0-yd)+f[i+1,j]*(1.0-yd)*xd+f[i,j+1]*(1.0-xd)*yd+f[i+1,j+1]*xd*yd
    #print xreq,yreq,xh,yh,i,j,xd,yd    
    return intval

#############################################################################

def initGaussboundcon(nptsgau=50 ) :

    ## Damping region size in grid points
    ##nptsgau = 50 #21
    decay = 0.24/nptsgau  #2.3/nptsgau
    
    xdist = NP.arange(1,nptsgau+1)
    damp = NP.exp( -((decay * (nptsgau - xdist))**2))
    
    leftdp   = damp.copy()
    rightdp  = damp[::-1].copy()
    bottomdp = damp[::-1].copy()
    topdp    = damp.copy()

    return nptsgau,leftdp,rightdp,bottomdp,topdp

#########################################################################

def calc_Kab_CPML(nptspml,gridspacing,dt,Npower,d0,
                       alpha_max_pml,K_max_pml,onwhere ) :

    # L = thickness of adsorbing layer
    if onwhere=="grdpts" :
        L = nptspml*gridspacing
        # distances 
        x = NP.arange(0.0,nptspml*gridspacing,gridspacing)
    elif onwhere=="halfgrdpts" :
        L = nptspml*gridspacing
        # distances 
        x = NP.arange(gridspacing/2.0,nptspml*gridspacing,gridspacing)
        
    d = d0 * (x/L)**Npower    
    alpha =  alpha_max_pml * (1.0 - (x/L)) # + 0.1 * alpha_max_pml ????

    K = 1.0 + (K_max_pml - 1.0) * (x/L)**Npower
    b = NP.exp( - (d / K + alpha) * dt )
    a = d * (b-1.0)/(K*(d+K*alpha))
    
    return K,a,b


###################################################################

def solveacoustic2D( inpar, ijsrc, vel, sourcetf, srcdomfreq, recpos ):

    if inpar["boundcond"]=="PML" :
        seism,psave       = solveacouwaveq2D_CPML( inpar, ijsrc, velmod, sourcetf, f0, recpos ) 

    elif inpar["boundcond"]=="GaussTap" :
        seism,psave = solveacouwaveq2D_GaussTaper( inpar, ijsrc, velmod, sourcetf, f0, recpos ) 

    elif inpar["boundcond"]=="ReflBou" :
        seism,psave  = solvewacouaveq2D_ReflBound( inpar, ijsrc, velmod, sourcetf, f0, recpos )

    return seism,psave

###################################################################

def solveacouwaveq2D_CPML( inpar, ijsrc, vel, sourcetf, srcdomfreq, recpos ) :

    assert(inpar["boundcond"]=="PML")

    #############################
    dh = inpar["dh"]
    dx = dh
    dz = dh
    dt = inpar["dt"]
    
    ##############################
    ## Check stability criterion
    ##############################
    maxvp = vel.max()
    print(" Stability criterion ???????????",(maxvp*dt*NP.sqrt(1/dx**2+1/dz**2)))
    assert(maxvp*dt*NP.sqrt(1/dx**2 + 1/dz**2) < 1.0)
    

    ##############################
    #   Parameters
    ##############################
    ## pml = c-pml, else reflsurf
    ## freesur = free surf for top
    if inpar["freesurface"]==True :
        print("freesurface at the top")
        freeboundtop=True
    else :
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

    print((" Size of PML layers in grid points: {} in x and {} in z".format(nptspml_x,nptspml_z)))
    
    Npower = 2.0
    assert Npower >= 1
    
    K_max_pml = 1.0
    alpha_max_pml = 2.0*NP.pi*(f0/2.0) 
    
    # reflection coefficient (INRIA report section 6.1)
    # http://hal.inria.fr/docs/00/07/32/19/PDF/RR-3471.pdf
    Rcoef = 0.0001  # for 20 nodes   ## 0.001 for a PML thickness of 10 nodes
       
    # thickness of the PML layer in meters
    thickness_pml_x = nptspml_x * dx
    thickness_pml_z = nptspml_z * dz
    
    # compute d0 from INRIA report section 6.1 http://hal.inria.fr/docs/00/07/32/19/PDF/RR-3471.pdf
    d0_x = - (Npower + 1) * maxvp * NP.log(Rcoef) / (2.0 * thickness_pml_x)
    d0_z = - (Npower + 1) * maxvp * NP.log(Rcoef) / (2.0 * thickness_pml_z)
    
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
    K_x     = NP.ones(nx)
    a_x     = NP.zeros(nx)
    b_x     = NP.ones(nx)
    
    # reverse coefficients to get increasing damping away from inner model
    # left boundary
    K_x[:nptspml_x]     = K_xpml[::-1] #[end:-1:1]
    a_x[:nptspml_x]     = a_xpml[::-1]
    b_x[:nptspml_x]     = b_xpml[::-1]

    # right boundary
    rightpml = nx-nptspml_x  # julia: nx-nptspml_x+1 
    K_x[rightpml:]     = K_xpml
    a_x[rightpml:]     = a_xpml
    b_x[rightpml:]     = b_xpml

    #######################
    #  z direction
    #######################
    K_z     = NP.ones(nz)
    a_z     = NP.zeros(nz)
    b_z     = NP.zeros(nz)

    # bottom 
    bottompml=nz-nptspml_z  # julia: nz-nptspml_z+1
    K_z[bottompml:]     = K_zpml
    a_z[bottompml:]     = a_zpml
    b_z[bottompml:]     = b_zpml

    # if PML also on top of model...
    if freeboundtop==False :
        # on grid
        K_z[:nptspml_z]     = K_zpml[::-1]
        a_z[:nptspml_z]     = a_zpml[::-1]
        b_z[:nptspml_z]     = b_zpml[::-1]

    #####################################################

    ## Arrays to export snapshots
    if inpar["savesnapshot"]==True :
        ntsave = inpar["ntimesteps"]//inpar["snapevery"]
        psave = NP.zeros((inpar["nx"],inpar["nz"],ntsave+1))
        tsave=1

    ## Arrays to return seismograms
    nrecs = recpos.shape[0]
    receiv = NP.zeros((inpar["ntimesteps"],nrecs))

    # Source time function
    #lensrctf = sourcetf.size
    isrc = ijsrc[0]
    jsrc = ijsrc[1]
    
    ## Initialize arrays
    dpdx = NP.zeros((nx,nz))
    dpdz = NP.zeros((nx,nz))
    pcur = NP.zeros((nx,nz))
    pold = NP.zeros((nx,nz))
    pnew = NP.zeros((nx,nz))
    
    ##--------------------------------------------
    # PML arrays
    # Arrays with size of PML areas would be sufficient and save memory,
    #   however allocating arrays with same size than model simplifies
    #   the code in the loops
    psi_dpdx = NP.zeros((nx,nz))
    psi_dpdz = NP.zeros((nx,nz))

    ##---------------------------------------------------
    ## to make the Numpy broadcast Hadamard product correct along x
    a_x = a_x.reshape(-1,1)
    b_x = b_x.reshape(-1,1)
    K_x = K_x.reshape(-1,1)

    ################################
    fact = vel**2 * (dt**2/dh**2)

       
    ## time loop
    print((" Time step: {}".format(dt)))

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
        dpdx = pcur[2:,1:-1]-2.0*pcur[1:-1,1:-1]+pcur[:-2,1:-1]
        dpdz = pcur[1:-1,2:]-2.0*pcur[1:-1,1:-1]+pcur[1:-1,:-2]
        
        # C-PML stuff 
        psi_dpdx[1:-1,1:-1] = b_x[1:-1] * psi_dpdx[1:-1,1:-1] + a_x[1:-1]*dpdx
        psi_dpdz[1:-1,1:-1] = b_z[1:-1] * psi_dpdz[1:-1,1:-1] + a_z[1:-1]*dpdz
        
        dpdx = dpdx / K_x[1:-1] + psi_dpdx[1:-1,1:-1]
        dpdz = dpdz / K_z[1:-1] + psi_dpdz[1:-1,1:-1]
        
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
            rec_press = bilinear_interp(pcur,dh,recpos[r,:])
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

def solveacouwaveq2D_ReflBound( inpar, ijsrc, vel, sourcetf, srcdomfreq, recpos ) :

    assert(inpar["boundcond"]=="ReflBou")

    #############################
    dh = inpar["dh"]
    dx = dh
    dz = dh
    dt = inpar["dt"]
    
    ##############################
    ## Check stability criterion
    ##############################
    maxvp = vel.max()
    assert (maxvp * dt*  NP.sqrt(1/dh**2 + 1/dh**2) < 1.0)

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
        psave = NP.zeros((inpar["nx"],inpar["nz"],ntsave+1))
        tsave=1

    ## Arrays to return seismograms
    nrecs = recpos.shape[0]
    receiv = NP.zeros((inpar["ntimesteps"],nrecs))

    # Source time function
    lensrctf = sourcetf.size
    isrc = ijsrc[0]
    jsrc = ijsrc[1]
    
    ## Initialize arrays
    dpdx = NP.zeros((nx,nz))
    dpdz = NP.zeros((nx,nz))
    pcur = NP.zeros((nx,nz))
    pold = NP.zeros((nx,nz))
    pnew = NP.zeros((nx,nz))

    ################################
    fact = vel**2 * (dt**2/dh**2)
       
    ## time loop
    print((" Time step: {}".format(dt)))

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
            rec_press = bilinear_interp(pcur,dh,recpos[r,:])
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

def solveacouwaveq2D_GaussTaper( inpar, ijsrc, vel, sourcetf, srcdomfreq, recpos ) :

    #############################
    dh = inpar["dh"]
    dx = dh
    dz = dh
    dt = inpar["dt"]
    
    ##############################
    ## Check stability criterion
    ##############################
    maxvp = vel.max()
    counum = maxvp * dt*  NP.sqrt(1/dh**2 + 1/dh**2)
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
    nptsgau = 40
    nptsgau,gtleft,gtright,gtbottom,gttop = initGaussboundcon(nptsgau=nptsgau )

    #####################################################
    
    ## Arrays to export snapshots
    if inpar["savesnapshot"]==True :
        ntsave = inpar["ntimesteps"]//inpar["snapevery"]
        psave = NP.zeros((inpar["nx"],inpar["nz"],ntsave+1))
        tsave=1

    ## Arrays to return seismograms
    nrecs = recpos.shape[0]
    receiv = NP.zeros((inpar["ntimesteps"],nrecs))

    # Source time function
    lensrctf = sourcetf.size
    isrc = ijsrc[0]
    jsrc = ijsrc[1]
    
    ## Initialize arrays
    dpdx = NP.zeros((nx,nz))
    dpdz = NP.zeros((nx,nz))
    pcur = NP.zeros((nx,nz))
    pold = NP.zeros((nx,nz))
    pnew = NP.zeros((nx,nz))

    ################################
    fact = vel**2 * (dt**2/dh**2)
       
    ## time loop
    print((" Time step: {}".format(dt)))

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
        
        ##---------------------------------------------
        ## Damp using Gaussian taper function
        # print("shapes: ",pnew[:nptsgau,:].shape,gtleft.reshape(-1,1).shape)
        pnew[:nptsgau,:]  *= gtleft.reshape(-1,1)
        pnew[-nptsgau:,:] *= gtright.reshape(-1,1)
        pnew[:,-nptsgau:] *= gtbottom.reshape(1,-1)

        pcur[:nptsgau,:]  *= gtleft.reshape(-1,1)
        pcur[-nptsgau:,:] *= gtright.reshape(-1,1)
        pcur[:,-nptsgau:] *= gtbottom.reshape(1,-1)

        if freeboundtop==False :
             pnew[:,:nptsgau]  *= gttop.reshape(-1,1)
             pcur[:,:nptsgau]  *= gttop.reshape(-1,1)
        ##----------------------------------------------
             
        # inject source
        pnew[isrc,jsrc] = pnew[isrc,jsrc] + dt**2*sourcetf[t]

        # assign the new pold and pcur
        pold[:,:] = pcur[:,:] ## using [:,:] there is an explicit copy
        pcur[:,:] = pnew[:,:]

        ##=================================================
        
        ##### receivers
        for r in range(nrecs) :
            rec_press = bilinear_interp(pcur,dh,recpos[r,:])
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
    t = NP.arange(0.0,nt*dt,dt) # s
    
    # space
    nx = 600
    nz = 400
    dh = 2.5 # m

    # source
    t0 = 0.03 # s
    f0 = 32.0 # 20.0 # Hz
    ijsrc = NP.array([280,290])
 
    # source type
    from sourcetimefuncs import gaussource, rickersource
    sourcetf = rickersource1D( t, t0, f0 )
    #sourcetf = gaussource1D( t, t0, f0 )

    ## velocity model
    velmod = 2000.0*NP.ones((nx,nz))
    # velmod[:,:40] = 2000.0
    # velmod[:,50:200] = 2500.0
    #velmod[50:80,50:70] = 3800.0
    
    ###############################
    nrec = 8
    recpos = NP.zeros((nrec,2))
    recpos[:,1] = 100.0
    recpos[:,0] = NP.linspace(200.0,nx*dh-200.0,nrec)

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
    inpar["boundcond"]= "GaussTap"  ## "PML", "GaussTap","ReflBou"

    #--------------------------
    import time
    t1 = time.time()
    
    seism,psave = solveacoustic2D( inpar, ijsrc, vel, sourcetf, srcdomfreq, recpos )
        
    t2 = time.time()
    print(("Solver time: {}".format(t2-t1)))

    
    ##############################
    ## save stuff
    import h5py as H5
    hf = H5.File("acu_imgs_python.h5","w")
    if inpar["savesnapshot"]==True :
        hf["press"] = psave
    hf["seism"] = seism
    hf["vel"] = velmod
    hf["srctf"] = sourcetf
    hf["dh"] = dh
    hf["dt"] = dt
    hf["nx"] = nx
    hf["nz"] = nz
    hf["recpos"] = recpos
    hf["srcij"] = ijsrc
    hf["snapevery"] = inpar["snapevery"]
    hf.close()

    return

#########################################################
####################################################


if __name__  == "__main__" :

    testacou()
