
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

# -*- coding: utf-8 -*-

import numpy as NP
import sys

#############################################################################
#############################################################################

from collections import namedtuple
SrcInterp = namedtuple("SrcInterp", ['x','z','imin','imax','jmin','jmax','wind'])
RecInterp = namedtuple("RecInterp", ['x','z','imin','imax','jmin','jmax','wind'])
MomentTensor = namedtuple("MomentTensor", ['Mxx','Mzz','Mxz'])

#############################################################################

# # kaiser window parameterized by alpha
# def kaiserdef(M,alpha) :
#     ## numpy.kaiser(M, alpha) ???
#     from numpy.dual import i0
#     if M == 1:
#         return np.array([1.])
#     n = arange(0, M)
#     alpha = (M-1)/2.0
#     return i0(alpha * sqrt(1-((n-alpha)/alpha)**2.0))/i0(float(alpha))

#############################################################################

def kaiser(x, x0, b, r) :
    # Kaiser window function
    #  r is window half with
    # Rule of thumb for finite diff.:
    #   b=4.14  b=6.31
    #   r = 4.0*dx
    w=NP.zeros(x.size)
    for i in range(x.size):
        if -r<=(x[i]-x0)<=r :
        ## Bessel function first kind: numpy.i0 
            den = 1.0/NP.i0(b)
            w[i] = den*NP.i0(b*(NP.sqrt(1 -((x[i]-x0)/r)**2)))
        else :
            w[i] = 0.0
    return w

#############################################################################

def setupreceivinterp(xstart,zstart,dx,dz,nx,nz,recpos) :

    # compute coefficients
    #  ix and iz contain the indices i_min,i_max and j_min,j_max of the 2D window

    nrecs = recpos.shape[0]
    recint = NP.zeros(nrecs,dtype="object")
    for r in range(nrecs):
        ## compute window
        xzwind,ix,iz = coeffsinc2D(xstart,zstart,dx,dz,recpos[r,0],recpos[r,1])
        ## check valid position 
        if ix[0]<1 :
            raise ValueError("setuprecinter(): Error receiver too close to boundaries! Quitting.")
        elif ix[1]>nx :
            raise ValueError("setuprecinter(): Error receiver too close to boundaries! Quitting.")
        elif iz[0]<1 :
            print("setuprecinter(): Receiver close to top boundary, using mirroring technique.")
            print("setuprecinter():    See Hicks 2002, Geophysics.")
            print("setuprecinter():  To be tested!!!")
            ##-------------------------
            nn = -iz[0]  # +2??
            if xzwind.shape[1]%2==0 :
                # tmp = xzwind[:,:nn]
                # tmp = tmp[:,::-1]
                # smallwin = xzwind[:,nn:]
                # mm = min(smallwin.shape[1],tmp.shape[1])
                # smallwin[:,:mm] = smallwin[:,:mm] - tmp[:,:mm]
                # xzwind = smallwin #xzwind[:,nn:]
                # iz[0] = 0
                
                xzwind = xzwind[:,nn:]
                iz[0] = 0 # 1 Julia
                
            else :
                # tmp = xzwind[:,:nn+1]
                # tmp = tmp[:,::-1]
                # smallwin = xzwind[:,nn:]
                # mm = min(smallwin.shape[1],tmp.shape[1])
                # smallwin[:,:mm] = smallwin[:,:mm] - tmp[:,:mm]
                # xzwind = smallwin #xzwind[:,nn:]
                # iz[0] = 0
               
                xzwind = xzwind[:,nn:]
                iz[0] = 0 # 1 Julia

        elif iz[1]>nz :
            raise ValueError("setuprecinter(): Error receiver too close to boundaries! Quitting.")
        #-------------
        recint[r] = RecInterp(x=recpos[r,0],z=recpos[r,1],imin=ix[0],imax=ix[1],
                              jmin=iz[0],jmax=iz[1],wind=xzwind)
        
        
        # import matplotlib.pyplot as plt
        # plt.figure(figsize=(10,4))
        # plt.subplot(121)
        # plt.title('REC {}, xzwind  xstart: {}, zstart: {}'.format(r,xstart,zstart))
        # plt.imshow(xzwind.T)
        # plt.xlabel("x")
        # plt.ylabel("z")
        # plt.colorbar()
        # plt.subplot(122)
        # plt.plot(xzwind[4,:])
        # plt.show()

    return recint

#############################################################################

def setupsourceinterp(xstart,zstart,dx,dz,nx,nz,srcpos,deriv_x=False, deriv_z=False ) :

    # compute coefficients
    #  ix and iz contain the indices i_min,i_max and j_min,j_max of the 2D window

    nsrc = srcpos.shape[0]
    if nsrc>1 :
        raise ValueError("Only one source allowed for now...")

 
    srcint =  NP.zeros(nsrc,dtype="object")
    for r in range(nsrc) :
        ## compute sinc window for each moment tensor component
        ## xstart+dx/2 staggered grid
        xzwind,ix,iz = coeffsinc2D(xstart,zstart,dx,dz,srcpos[r,0],srcpos[r,1],
                                   deriv_x=deriv_x, deriv_y=deriv_z )
        #print("src ix iz:",ix,iz)          
        ## check valid position 
        if ix[0]<1 :
            raise ValueError("setupsrcinter(): Error source too close to boundaries! Quitting.")
        elif ix[1]>nx :
            raise ValueError("setupsrcinter(): Error source too close to boundaries! Quitting.")
        elif iz[0]<1 :
            print("setupsrcinter(): r: {} Source close to top boundary, using mirroring technique.".format(r))
            print("setupsrcinter():    See Hicks 2002, Geophysics.")
            print("setupsrcinter():  To be tested!!!")
            ##-------------------------
            nn = -iz[0]  # +2??
           
            if xzwind.shape[1]%2==0 :
                # print("Z EVEN  nn: {}, ix: {}, iz: {} ".format(nn,ix,iz))
                # tmp = xzwind[:,:nn]
                # tmp = tmp[:,::-1]
                # smallwin = xzwind[:,nn:]
                # mm = min(smallwin.shape[1],tmp.shape[1])
                # smallwin[:,:mm] = smallwin[:,:mm] - tmp[:,:mm]
                # xzwind = smallwin #xzwind[:,nn:]
                # iz[0] = 0 # 1 Julia

                xzwind = xzwind[:,nn:]
                iz[0] = 0 # 1 Julia
                
            else :
                # print("Z ODD  nn: {}, ix: {}, iz: {} ".format(nn,ix,iz))
                # tmp = xzwind[:,:nn+1]
                # tmp = tmp[:,::-1]
                # smallwin = xzwind[:,nn:]
                # mm = min(smallwin.shape[1],tmp.shape[1])
                # smallwin[:,:mm] = smallwin[:,:mm] - tmp[:,:mm]
                # xzwind = smallwin #xzwind[:,nn:]
                # iz[0] = 0

                xzwind = xzwind[:,nn:]
                iz[0] = 0 # 1 Julia

        elif iz[1]>nz :
            raise ValueError("setupsrcinter(): Error source too close to boundaries! Quitting.")
        #-------------
        srcint[r] = SrcInterp(x=srcpos[r,0],z=srcpos[r,1],imin=ix[0],imax=ix[1],
                              jmin=iz[0],jmax=iz[1],wind=xzwind)

        # import matplotlib.pyplot as plt
        # plt.figure(figsize=(10,4))
        # plt.subplot(121)
        # plt.title('SRC {}, xzwind  xstart: {}, zstart: {}'.format(r,xstart,zstart))
        # plt.imshow(xzwind.T)
        # plt.xlabel("x")
        # plt.ylabel("z")
        # plt.colorbar()
        # plt.subplot(122)
        # plt.plot(xzwind[4,:])
        # plt.show()

    return srcint

#############################################################################

def coeffsinc(xstart,dx,xcenter, deriv=False,npts=4,beta=4.14) :
    ## Coefficients for sinc interpolation
    ##  in 1D
    ##  xstart is the x coordinate of first node in the regular grid    
    ##  
    ## beta = 4.14 Hicks 2002, Geophysics
    ###
    ### Julia:  sinc(x) =
    ###    \sin(\pi x) / (\pi x) if x \neq 0, and 1 if x = 0
    ###
    rx = npts*dx
    ## Assuming x from grid starts xstart
    xh = (xcenter-xstart)/dx
    ix = int(NP.floor(xh+1))
    if (xcenter-xstart)%dx == 0.0 :
        ixsta = ix-npts
        ixend = ix+npts
    else :
        ixsta = ix-npts+1
        ixend = ix+npts

    x = NP.asarray( [xstart+dx*(i-1) for i in range(ixsta,ixend+1)] )
    
    if deriv==False :
        # interpolating sinc(x) = sin(pi*x)/(pi*x) see Julia and Numpy definition
        intrpsinc = NP.sinc((x-xcenter)/dx)
        
    elif deriv==True :
        # derivative 
        intrpsinc = NP.zeros(x.size)
        for i in range(x.size) :
            ## Hicks 2002 Geophysics
            tmp = kmax*(x[i]-xsource)
            numer = tmp*NP.cos(tmp)-NP.sin(tmp)
            denom = tmp**2
            if denom==0.0 :
                intrpsinc[i]=0.0
            else :
                intrpsinc[i] = kmax * numer / denom

    # apply Kaiser windowing
    kaix = kaiser(x,xcenter,beta,rx)
    windx = kaix * intrpsinc
    # return also indices of window
    #print("windx: {}".format(windx))
    return windx,NP.array([ixsta-1,ixend-1]) #-1 Python

#############################################################################

def coeffsinc2D(xstart,ystart,dx,dy,xcenter,ycenter,deriv_x=False,deriv_y=False, 
                npts=4,beta=4.14) :
    ## Calculate the 2D array of coefficients
    windx,ixstaend = coeffsinc(xstart,dx,xcenter, deriv=deriv_x, npts=4, beta=4.14)
    windy,iystaend = coeffsinc(ystart,dy,ycenter, deriv=deriv_y, npts=4, beta=4.14)
    # tensor product of x and y
    xywind = NP.outer(windx, windy) # transposed, outer prod
    # return also indices of window
    return xywind,ixstaend,iystaend

#############################################################################

def gaussource1D( t, t0, f0 ) :

    # boh = f0 .* (t-t0)
    # source = -8.0*boh.*exp( -boh.^2/(4.0*f0)^2 )
    a = (NP.pi*f0)**2
    source = - 8.0*a*(t-t0)*NP.exp( -a*(t-t0)**2 )
    
    return source

#############################################################################

def derivRickersource1D( t, t0, f0 ):
    # first derivative of a Ricker wavelet
    ## SymPY
    ## import sympy as SY
    ## f0,t,t0 = SY.symbols('f0 t t0')
    ## w = (1.0-2.0*(SY.pi*f0*(t-t0))**2) * SY.exp(-(SY.pi*f0*(t-t0))**2)
    ## dw = SY.simplify(SY.diff(w,t))

    a = NP.pi**2 * f0**2
    w = a*(-4.0*t+4.0*t0+2.0*(t-t0)*(2.0*a*(t-t0)**2-1.0))*NP.exp(-a*(t-t0)**2)
    return w


def rickersource1D( t, t0, f0 ):
    ## Ricker
    b = (NP.pi*f0*(t-t0))**2
    w = (1.0-2.0*b)*NP.exp(-b)
    return w

#############################################################################

# def bilinear_interp(f,hgrid, pt):
#     xreq=pt[0]
#     yreq=pt[1]
#     xh=xreq/hgrid
#     yh=yreq/hgrid
#     i=int(NP.floor(xh)) # index starts from 0 so no +1
#     j=int(NP.floor(yh)) # index starts from 0 so no +1
#     xd=xh-i
#     yd=yh-j
#     intval=f[i,j]*(1.0-xd)*(1.0-yd)+f[i+1,j]*(1.0-yd)*xd+f[i,j+1]*(1.0-xd)*yd+f[i+1,j+1]*xd*yd
#     #print xreq,yreq,xh,yh,i,j,xd,yd    
#     return intval

#############################################################################

def calc_dKa_CPML(nptspml,gridspacing,dt,Npower,d0,
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
    K = 1.0 + (K_max_pml - 1.0) * (x/L)**Npower
    #from Stephen Gedney's unpublished class notes for class EE699, lecture 8, slide 8-2
    alpha =  alpha_max_pml * (1.0 - (x/L)) + 0.1 * alpha_max_pml   
    b = NP.exp( - (d / K + alpha) * dt )
    a = d * (b-1.0)/(K*(d+K*alpha))
    
    return d,K,alpha,b,a


#############################################################################

def solveelawaveq2D_CPML(inpar, rockprops, inpsrcpos, sourcetf, momtens, srcdomfreq, inprecpos) :
    #  
    # Wave elastic staggered grid 2D solver 
    #
    # Andrea Zunino, 4/5-2017
    #
    # Staggered grid with equal spacing, A. Levander (1988)
    # Second order time, fourth order space
    #
    # Convolutionary Perfectly Matched Layer (C-PML) boundary conditions
    #
    #   References:
    #       Levander A. (1988), Fourth-order finite-difference P-SV seismograms, Geophysics.
    #       Komatitsch D. and Martin R. (2007), An unsplit convolutional
    #            perfectly matched layer improved at grazing incidence for the
    #            seismic wave equation, Geophysics.
    #       Robertsson, J.O. (1996) Numerical Free-Surface Condition for
    #            Elastic/Viscoelastic Finite-Difference Modeling in the Presence
    #            of Topography, Geophysics.
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
    print(" Grid dimensions: {} x {}".format(dx*inpar['nx'],dz*inpar['nz']))

    maxx = dx*inpar['nx']
    maxz = dz*inpar['nz']

    if (inprecpos[:,0]<0.0).any() :
        raise ValueError("inprecpos[:,0]<0.0")
    elif (inprecpos[:,0]>maxx).any() :
        raise ValueError("inprecpos[:,0]>maxx")
    if (inprecpos[:,1]<0.0).any() :
        raise ValueError("inprecpos[:,1]<0.0")
    elif (inprecpos[:,1]>maxz).any() :
        raise ValueError("inprecpos[:,1]>maxz")

    if (inpsrcpos[:,0]<0.0).any() :
        raise ValueError("inpsrcpos[:,0]<0.0")
    elif (inpsrcpos[:,0]>maxx).any() :
        raise ValueError("inpsrcpos[:,0]>maxx")
    if (inpsrcpos[:,1]<0.0).any() :
        raise ValueError("inpsrcpos[:,1]<0.0")
    elif (inpsrcpos[:,1]>maxz).any() :
        raise ValueError("inpsrcpos[:,1]>maxz")

    ##############################
    ## Check stability criterion
    ##############################
    maxvp = ( NP.sqrt( (lamb+2.0*mu)/rho )).max()
    Courant_number = maxvp * inpar["dt"] * NP.sqrt(1.0/dx**2 + 1.0/dz**2)
    print(" Courant number: {}".format(Courant_number))
    if Courant_number > 1.0 :
        print(" The stability criterion is violated. \nCourant_number={} \n Quitting.\n".format(Courant_number))
        return

        # inpar["dt"] = 0.66 / (maxvp * NP.sqrt(1.0/dx**2 + 1.0/dz**2))
        # print(" Setting dt to: {}".format(inpar["dt"]))
    
    #minvp = minimum( sqrt( (lambd+2.0*mu)./rho ))

    ###############################
    ## Check dispersion criterium
    ###############################
    fmaxsource = srcdomfreq # very approx.... # Nyquist: 1/(2.0*dt)???
    arbfact = 1.5 # * max freq to be safe... boh?
    ngridptsperwavlen = 8 # standard value for this setup
    idxmupos = NP.where(mu>0.0)
    if idxmupos[0].size > 0 :
        minvs = ( NP.sqrt(mu[idxmupos]/rho[idxmupos]) ).min()
    else :
        minvs =  1e30
    minvp = ( NP.sqrt((lamb+2.0*mu)/rho) ).min()
    minvel = min(minvs,minvp)
    maxalloweddh = minvel/(ngridptsperwavlen*arbfact*fmaxsource)
    print(" Max. allowed dh: {}, minvel {}".format(maxalloweddh,minvel))
    if dh>=maxalloweddh :
        raise ValueError(" The dispersion criterion is violated.\nmaxallowed dh = {}\n Quitting.\n".format(maxalloweddh))


    ##############################
    #   Parameters
    ##############################

    harmonicaver_mu = False
    # will be modified, so copy...
    srcpos = inpsrcpos.copy()
    recpos = inprecpos.copy()

    f0 = srcdomfreq # source
     
    ## pml = c-pml, else reflsurf
    ## freesur = free surf for top
    topbound    = inpar["topbound"]
    leftbound   = inpar["leftbound"]
    rightbound  = inpar["rightbound"]
    bottombound = inpar["bottombound"]
    assert (topbound in ["freesurf","pml","reflsurf"])
    if topbound=="freesurf" :
        freeboundtop=True
    else :
        freeboundtop=False
    assert (leftbound in ["pml","reflsurf"])
    assert (rightbound in ["pml","reflsurf"])
    assert (bottombound in ["pml","reflsurf"])

    print("\n Boundary conditions: ")
    print("   +------- {:8s} -------|--> x ".format(topbound))
    print("   |                        |")
    print("  {:8s}                 {:8s} ".format(leftbound,rightbound))
    print("   |                        | ")
    print("   |------- {:8s} -------|".format(bottombound))
    print("   | ")
    print("   v z \n")
    
    print(" Grid spacing dx,dz: {}".format(dx,dz))

    ## keep it simple for now...
    nx = inpar["nx"]
    nz = inpar["nz"]
    print(" Size of grid: nx= {} and nz={} ".format(nx,nz))

    ##############################
    #   Parameters for PML
    ##############################
    factorPMLlayers = 1.2
    print(" factorPMLlayers: {}".format(factorPMLlayers))
    npml_x = int(NP.ceil(factorPMLlayers*(maxvp/f0)/inpar["dh"]))
    npml_z = int(NP.ceil(factorPMLlayers*(maxvp/f0)/inpar["dh"]))
    print(" Size of PML layers in grid points: {} in x and {} in z".format(npml_x,npml_z))
    
    


    #######################################################

    def setupPMLbounds(inprkprp):
        # Add patches of PML layers to the sides of
        #  the grid depending on boundary conditions
        # 
        rkprp = inprkprp.copy()
        if topbound=="pml" :
            # add padding on top
            toppatch = NP.tile(inprkprp[:,0],(npml_z,1)).T
            rkprp = NP.hstack((toppatch,rkprp))
            
        if rightbound=="pml" :
            # add padding on the right
            smallright = NP.tile(inprkprp[-1,:],(npml_x,1))
            if topbound=="pml" :
                uprightcor = inprkprp[-1,0]*NP.ones((npml_x,npml_z))
                rightpatch = NP.hstack((uprightcor,smallright))
            else :
                rightpatch = smallright
            rkprp = NP.vstack((rkprp,rightpatch))

        if bottombound=="pml" :
            # add padding on the bottom
            smallbottom = NP.tile(inprkprp[:,-1],(npml_z,1)).T
            if rightbound=="pml" :
                lowrightcor = inprkprp[-1,-1]*NP.ones((npml_x,npml_z))
                bottompatch = NP.vstack((smallbottom,lowrightcor))
            else :
                bottompatch = smallbottom
            rkprp = NP.hstack((rkprp,bottompatch))

        if leftbound=="pml" : 
            # add padding on the left
            smallleft = NP.tile(inprkprp[-1,:],(npml_x,1))
            if bottombound=="pml" and topbound=="pml" :
                lowleftcor = inprkprp[0,-1]*NP.ones((npml_x,npml_z))
                upleftcor = inprkprp[0,0]*NP.ones((npml_x,npml_z))
                leftpatch = NP.hstack((upleftcor,smallleft,lowleftcor))
            elif bottombound=="pml" :
                lowleftcor = inprkprp[0,-1]*NP.ones((npml_x,npml_z))
                leftpatch = NP.hstack((smallleft,lowleftcor))
            elif topbound=="pml" :
                upleftcor = inprkprp[0,0]*NP.ones((npml_x,npml_z))
                leftpatch = NP.hstack((upleftcor,smallleft))
            else :
                leftpatch = smallleft
            rkprp = NP.vstack((rkprp,leftpatch))
        #----
        return rkprp
        
    ##-----------------------------------------------------

    # xstart = 0.0
    # zstart = 0.0
    if topbound=="pml" :
        ## shift z position because of PML layers
        srcpos[:,1] = srcpos[:,1] + npml_z*dz
        recpos[:,1] = recpos[:,1] + npml_z*dz
        # zstart = 0.0 #-npml_z*dz

    if leftbound=="pml" : 
        ## shift x position because of PML layers
        srcpos[:,0] = srcpos[:,0] + npml_x*dx
        recpos[:,0] = recpos[:,0] + npml_x*dx
        # xstart = 0.0 #-npml_x*dx
                
    # Increase grid size for PML
    lamb = setupPMLbounds(lamb)
    mu = setupPMLbounds(mu)
    rho = setupPMLbounds(rho)

    # nxpml = nx + 2*npml_x
    # nzpml = nz + 2*npml_z
    nx = lamb.shape[0]
    nz = lamb.shape[1]
    print(" New size of grid: nx= {} and nz= {} ".format(nx,nz))
    


    
    ##-----------------------------------------------------
    import matplotlib.pyplot as plt

    plt.figure()
    plt.subplot(221)
    plt.title('lamb')
    plt.imshow(lamb.T)
    plt.colorbar()
    plt.subplot(222)
    plt.title('mu')
    plt.imshow(mu.T)
    plt.colorbar()
    plt.subplot(223)
    plt.title('rho')
    plt.imshow(rho.T)
    plt.colorbar()
    plt.show(block=False)
    



    #--------------------------------
    Npower = 2.0
    assert Npower >= 1
    
    K_max_pml = 1.0
    alpha_max_pml = 2.0*NP.pi*(f0/2.0) 
    
    # reflection coefficient (INRIA report section 6.1)
    # http://hal.inria.fr/docs/00/07/32/19/PDF/RR-3471.pdf
    Rcoef = 0.001
       
    # thickness of the PML layer in meters
    thickness_pml_x = npml_x * dx
    thickness_pml_z = npml_z * dz
    
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
    d_xpml,K_xpml,alpha_xpml,b_xpml,a_xpml = calc_dKa_CPML(npml_x,dh,dt,Npower,d0_x,alpha_max_pml,K_max_pml,"grdpts")
    # damping profile at half the grid points
    d_xpml_half,K_xpml_half,alpha_xpml_half,b_xpml_half,a_xpml_half = calc_dKa_CPML(npml_x,dh,dt,Npower,d0_x,alpha_max_pml,K_max_pml,"halfgrdpts")
    
    # --- damping in the z direction ---
    # assuming the number of grid points for PML is the same on
    # both sides    
    # damping profile at the grid points
    d_zpml,K_zpml,alpha_zpml,b_zpml,a_zpml = calc_dKa_CPML(npml_z,dh,dt,Npower,d0_z,alpha_max_pml,K_max_pml,"grdpts")
    # damping profile at half the grid points
    d_zpml_half,K_zpml_half,alpha_zpml_half,b_zpml_half,a_zpml_half = calc_dKa_CPML(npml_z,dh,dt,Npower,d0_z,alpha_max_pml,K_max_pml,"halfgrdpts")
    
    #######################
    #  x direction
    #######################
    d_x     = NP.zeros(nx)
    K_x     = NP.ones(nx)
    alpha_x = NP.zeros(nx)
    b_x     = NP.zeros(nx)
    a_x     = NP.zeros(nx)
    d_x_half     = NP.zeros(nx)
    K_x_half     = NP.ones(nx)
    alpha_x_half = NP.zeros(nx)
    b_x_half     = NP.zeros(nx)
    a_x_half     = NP.zeros(nx)

    # reverse coefficients to get increasing damping away from inner model
    if leftbound=="pml" :
        # left boundary
        d_x[:npml_x]     = d_xpml[::-1]
        K_x[:npml_x]     = K_xpml[::-1] #[end:-1:1]
        alpha_x[:npml_x] = alpha_xpml[::-1]
        b_x[:npml_x]     = b_xpml[::-1]
        a_x[:npml_x]     = a_xpml[::-1]
        # half stuff...
        # reverse coefficients to get increasing damping away from inner model
        #  One less element on left boundary (end-1...)
        #    because of the staggered grid 
        d_x_half[:npml_x-1]     = d_xpml_half[-2::-1] # julia: [end-1:-1:1]
        K_x_half[:npml_x-1]     = K_xpml_half[-2::-1]
        alpha_x_half[:npml_x-1] = alpha_xpml_half[-2::-1]
        b_x_half[:npml_x-1]     = b_xpml_half[-2::-1]
        a_x_half[:npml_x-1]     = a_xpml_half[-2::-1]

    if rightbound=="pml" :
        # right boundary
        rightpml = nx-npml_x  # julia: nx-npml_x+1 
        d_x[rightpml:]     = d_xpml
        K_x[rightpml:]     = K_xpml
        alpha_x[rightpml:] = alpha_xpml
        b_x[rightpml:]     = b_xpml
        a_x[rightpml:]     = a_xpml 
        # half stuff
        # right boundary
        d_x_half[rightpml:]     = d_xpml_half
        K_x_half[rightpml:]     = K_xpml_half
        alpha_x_half[rightpml:] = alpha_xpml_half
        b_x_half[rightpml:]     = b_xpml_half
        a_x_half[rightpml:]     = a_xpml_half 


    #######################
    #  z direction
    #######################
    d_z     = NP.zeros(nz)
    K_z     = NP.ones(nz)
    alpha_z = NP.zeros(nz)
    b_z     = NP.zeros(nz)
    a_z     = NP.zeros(nz)
    # half stuff...
    d_z_half     = NP.zeros(nz)
    K_z_half     = NP.ones(nz)
    alpha_z_half = NP.zeros(nz)
    b_z_half     = NP.zeros(nz)
    a_z_half     = NP.zeros(nz)

    if bottombound=="pml" :
        # bottom 
        bottompml=nz-npml_z  # julia: nz-npml_z+1
        d_z[bottompml:]     = d_zpml
        K_z[bottompml:]     = K_zpml
        alpha_z[bottompml:] = alpha_zpml
        b_z[bottompml:]     = b_zpml
        a_z[bottompml:]     = a_zpml
        # reverse coefficients to get increasing damping away from inner model  
        # bottom
        d_z_half[bottompml:]     = d_zpml_half
        K_z_half[bottompml:]     = K_zpml_half
        alpha_z_half[bottompml:] = alpha_zpml_half
        b_z_half[bottompml:]     = b_zpml_half
        a_z_half[bottompml:]     = a_zpml_half

    # if PML also on top of model...
    if topbound=="pml" :
        # on grid
        d_z[:npml_z]     = d_zpml[::-1]
        K_z[:npml_z]     = K_zpml[::-1]
        alpha_z[:npml_z] = alpha_zpml[::-1]
        b_z[:npml_z]     = b_zpml[::-1]
        a_z[:npml_z]     = a_zpml[::-1]
        # half
        #  One less element on top boundary (end-1...)
        #    because of the staggered grid
        d_z_half[:npml_z-1]     = d_zpml_half[-2::-1]
        K_z_half[:npml_z-1]     = K_zpml_half[-2::-1]
        alpha_z_half[:npml_z-1] = alpha_zpml_half[-2::-1]
        b_z_half[:npml_z-1]     = b_zpml_half[-2::-1]
        a_z_half[:npml_z-1]     = a_zpml_half[-2::-1]
        

    #####################################################

    ## Arrays to export snapshots
    if inpar["savesnapshot"]==True :
        ntsave = inpar["ntimesteps"]//inpar["snapevery"]
        ## using original inpar['nx'], etc, to exclude
        ##   PML padding layers
        vxsave = NP.zeros((inpar['nx'],inpar['nz'],ntsave+1))
        vzsave = NP.zeros((inpar['nx'],inpar['nz'],ntsave+1))
        tsave=1


    ## Arrays to return seismograms
    nrecs = recpos.shape[0]
    receiv = NP.zeros((inpar["ntimesteps"],nrecs,2))


    ########################################
    ##  Setup interpolation for receivers  #
    ########################################
    ## separate for vx and vz because of the staggered grid,
    ##   so  xstart,ystart is different,
    ## vx is on integer grid and vz on half,half grid
    reit_vx = setupreceivinterp(0.0,   0.0,   dx,dz,nx,nz,recpos)
    reit_vz = setupreceivinterp(dx/2.0,dz/2.0,dx,dz,nx,nz,recpos)


    ########################################
    ##  Setup interpolation for sources    #
    ########################################
    ## separate for vx and vz because of the staggered grid,
    ##   so  xstart,ystart is different,
    ## vx is on integer grid and vz on half,half grid
    srcit_Mxx = setupsourceinterp(dx/2, 0.0, dx,dz,nx,nz,srcpos,False, False)
    srcit_Mzz = setupsourceinterp(dx/2, 0.0, dx,dz,nx,nz,srcpos,False, False)
    srcit_Mxz = setupsourceinterp(0.0, dz/2, dx,dz,nx,nz,srcpos,False, False)
    nsources = srcpos.shape[0]


    # Source time function
    lensrctf = sourcetf.size
 
    
    ## Initialize arrays
    vx = NP.zeros((nx,nz))
    vz = NP.zeros((nx,nz))
    Txx = NP.zeros((nx,nz))
    Tzz = NP.zeros((nx,nz))
    Txz = NP.zeros((nx,nz))

    ## derivatives
    Dxb_Txx = NP.zeros((nx-3,nz-3))
    Dzb_Txz = NP.zeros((nx-3,nz-3))
    Dxf_Txz = NP.zeros((nx-3,nz-3))
    Dzf_Tzz = NP.zeros((nx-3,nz-3))
    Dxf_vx  = NP.zeros((nx-3,nz-3))
    Dzb_vz  = NP.zeros((nx-3,nz-3))
    Dzf_vx  = NP.zeros((nx-3,nz-3))
    Dxb_vz  = NP.zeros((nx-3,nz-3))


    # PML arrays
    # Arrays with size of PML areas would be sufficient and save memory,
    #   however allocating arrays with same size than model simplifies
    #   the code in the loops
    psi_DxTxx = NP.zeros((nx,nz))
    psi_DzTxz = NP.zeros((nx,nz))
    
    psi_DxTxz = NP.zeros((nx,nz))
    psi_DzTzz = NP.zeros((nx,nz))
    
    psi_DxVx = NP.zeros((nx,nz))
    psi_DzVz = NP.zeros((nx,nz))

    psi_DzVx = NP.zeros((nx,nz))
    psi_DxVz = NP.zeros((nx,nz))

    
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

    fact = 1.0/(24.0*inpar["dh"])
 
    ###############################################################
    # pre-interpolate properties at half distances between nodes
    ###############################################################

    if inpar['interprockprop']==True :
        print(" Interpolating rock properties at half grid spacing.")
        # rho_ihalf_jhalf (nx-1,ny-1) ??
        rho_ihalf_jhalf = rho[1:,1:]   #(rho[1:,1:]+rho[1:,:-1]+rho[:-1,1:]+rho[:-1,:-1])/4.0
        # mu_ihalf (nx-1,ny) ??
        # mu_jhalf (nx,ny-1) ??
        if harmonicaver_mu==True :
            mu_ihalf = 1.0/( 1.0/mu[1:,:] + 1.0/mu[:-1,:] )
            mu_jhalf = 1.0/( 1.0/mu[:,1:] + 1.0/mu[:,:-1] )
        else :
            mu_ihalf = (mu[1:,:]+mu[:-1,:])/2.0 ###?????
            mu_jhalf = (mu[:,1:]+mu[:,:-1])/2.0 ###?????
        # lamb_ihalf (nx-1,ny) ??
        lamb_ihalf = lamb[1:,:]  #(lamb[1:,:]+lamb[:-1,:])/2.0 ###?????

    else :
        rho_ihalf_jhalf = rho[:-1,:-1]
        mu_ihalf        = mu[:-1,:]
        mu_jhalf        = mu[:,:-1]
        lamb_ihalf      = lamb[:-1,:]

    ##======================================##
   
    ## to make the Numpy broadcast Hadamard product correct along x
    a_x = a_x.reshape(-1,1)
    b_x = b_x.reshape(-1,1)
    K_x = K_x.reshape(-1,1)
    a_x_half = a_x_half.reshape(-1,1)
    b_x_half = b_x_half.reshape(-1,1)
    K_x_half = K_x_half.reshape(-1,1)
     
    ##======================================##
    ## time loop
    timearray = NP.array([])
    dt = inpar["dt"]
    print(" Time step: {}".format(dt))
    for t in range(inpar["ntimesteps"]) :

        timearray = NP.append(timearray,t*dt)
        
        if t%10==0 :
            sys.stdout.write("\r Time step {} of {}, currently at {:.6f} ".format(t,inpar["ntimesteps"],(t*dt)))
            sys.stdout.flush()
            #print("\r Time step {} of {}".format(t,inpar["ntimesteps"]))

        
        ## Inject the source 
        if t<=lensrctf :
            for s in range(nsources) :
                
                if inpar["sourcetype"]=="MomTensor":
                    ## interpolate source at nearby nodes using sinc method
                    Txx[srcit_Mxx[s].imin:srcit_Mxx[s].imax+1,srcit_Mxx[s].jmin:srcit_Mxx[s].jmax+1] = Txx[srcit_Mxx[s].imin:srcit_Mxx[s].imax+1,srcit_Mxx[s].jmin:srcit_Mxx[s].jmax+1] + (srcit_Mxx[s].wind * momtens.Mxx) * sourcetf[t] * dt
                
                    Tzz[srcit_Mzz[s].imin:srcit_Mzz[s].imax+1,srcit_Mzz[s].jmin:srcit_Mzz[s].jmax+1] = Tzz[srcit_Mzz[s].imin:srcit_Mzz[s].imax+1,srcit_Mzz[s].jmin:srcit_Mzz[s].jmax+1] + (srcit_Mzz[s].wind * momtens.Mzz) * sourcetf[t] * dt
                
                    Txz[srcit_Mxz[s].imin:srcit_Mxz[s].imax+1,srcit_Mxz[s].jmin:srcit_Mxz[s].jmax+1] = Txz[srcit_Mxz[s].imin:srcit_Mxz[s].imax+1,srcit_Mxz[s].jmin:srcit_Mxz[s].jmax+1] + (srcit_Mxz[s].wind * momtens.Mxz) * sourcetf[t] * dt 

                    
                elif inpar["sourcetype"]=="ExtForce":
                    # ExtForce['x'] = # dt/rho *
                    # Extforce['z'] =  
                    raise ValueError("ExtForce source not yet implemented! .Exiting.")

                else :
                    print("Error, source type badly defined.")
                    return
                
       
      
        #########################################
        # update velocities from stresses
        #########################################

        if freeboundtop==True :
            ## from Robertsson 1996
            ## j=0,1

            ### Vx
            Dxb_Txx = fact * ( Txx[:-3,:2] -27.0*Txx[1:-2,:2] +27.0*Txx[2:-1,:2] -Txx[3:,:2] )
            Dzb_Txz = fact * ( -Txz[2:-1,1:3] +27.0*Txz[2:-1,:2] +27.0*Txz[2:-1,:2] -Txz[2:-1,1:3] )

            # vx
            vx[2:-1,:2] = vx[2:-1,:2] + (dt/rho[2:-1,:2]) * (Dxb_Txx + Dzb_Txz)

            ###---------------------------------------------------
            # Vz
            Dxf_Txz = fact * ( Txz[:-3,:2] -27.0*Txz[1:-2,:2] +27.0*Txz[2:-1,:2] -Txz[3:,:2] )
            Dzf_Tzz = fact * ( -Tzz[1:-2,2:4] +27.0*Tzz[1:-2,1:3] +27.0*Tzz[1:-2,1:3] -Tzz[1:-2,2:4] )
    
            # update velocity (rho has been interpolated in advance)
            # rho_ihalf_jhalf[1:-1,1:-1] because its size is (nx-1,ny-1) 
            vz[1:-2,:2] = vz[1:-2,:2] + (dt/rho_ihalf_jhalf[1:-1,:2]) * (Dxf_Txz + Dzf_Tzz)
            
            ## END free surface
            ##============================================

        ## space loops excluding boundaries

        #-----------------------------------------------------
        ### Vx
        Dxb_Txx = fact * ( Txx[:-3,2:-1] -27.0*Txx[1:-2,2:-1] +27.0*Txx[2:-1,2:-1] -Txx[3:,2:-1] )
        Dzb_Txz = fact * ( Txz[2:-1,:-3] -27.0*Txz[2:-1,1:-2] +27.0*Txz[2:-1,2:-1] -Txz[2:-1,3:] )

        # C-PML stuff 
        psi_DxTxx[2:-1,2:-1] = b_x[2:-1] * psi_DxTxx[2:-1,2:-1] + a_x[2:-1] * Dxb_Txx
        psi_DzTxz[2:-1,2:-1] = b_z[2:-1] * psi_DzTxz[2:-1,2:-1] + a_z[2:-1] * Dzb_Txz
        
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
            ## from Robertsson 1996
            ## j=0

            # Txx,Tzz
            Dxf_vx = fact * (vx[:-3,0] -27.0*vx[1:-2,0] +27.0*vx[2:-1,0] -vx[3:,0])
            #Dzb_vz = -(1.0-2.0*mu_ihalf[1:-1,0]/lamb_ihalf[1:-1,0])*Dxf_vx
            ###  <<<<<<  NEW v     OLD /\
            Dzb_vz = -(1.0-2.0*mu_ihalf[1:-1,0]/(lamb_ihalf[1:-1,0]+2.0*mu_ihalf[1:-1,0]))*Dxf_vx
                
            # Txx
            Txx[1:-2,0] = Txx[1:-2,0] + (lamb_ihalf[1:-1,0]+2.0*mu_ihalf[1:-1,0]) * dt * Dxf_vx + lamb_ihalf[1:-1,0] * dt * Dzb_vz
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
                         

        #######################################################################################
        ##### receivers
        for r in range(nrecs) :

            #rec_vx = bilinear_interp(vx,dh,recpos[r,:])
            # ## vz at half grid points, so recpos - dh/2.0 in both x an z
            #rec_vz = bilinear_interp(vz,dh,recpos[r,:]- dh/2.0 ) # - dh/2.0)
            ##--------------------------------------------------
            ## interpolate seismograms using sinc method
            rec_vx = NP.sum( reit_vx[r].wind * vx[reit_vx[r].imin:reit_vx[r].imax+1,
                                             reit_vx[r].jmin:reit_vx[r].jmax+1] )
            rec_vz = NP.sum( reit_vz[r].wind * vz[reit_vz[r].imin:reit_vz[r].imax+1,
                                             reit_vz[r].jmin:reit_vz[r].jmax+1] )
            receiv[t,r,:] = NP.r_[rec_vx, rec_vz]

            #---------------------------------------------------
            if inpar["seismogrkind"] == "displacement" :
                # integrate to get displacement
                receiv[t,r,0] = NP.trapz(receiv[:,r,0],dx=dt)
                receiv[t,r,1] = NP.trapz(receiv[:,r,1],dx=dt)

                
        #### save snapshots
        if (inpar["savesnapshot"]==True) and (t%inpar["snapevery"]==0) :
            if leftbound=="pml" :
                ixsmin = npml_x
                ixsmax = npml_x+inpar['nx']
            else :
                ixsmin = 0
                ixsmax = inpar['nx']
            if topbound=="pml" :
                izsmin = npml_z
                izsmax = npml_z+inpar['nz']
            else :
                izsmin = 0
                izsmax = inpar['nz']
            ##------------------------------------
            vxsave[:,:,tsave] = vx[ixsmin:ixsmax,izsmin:izsmax]
            vzsave[:,:,tsave] = vz[ixsmin:ixsmax,izsmin:izsmax]
            tsave=tsave+1
    
    ##### End time loop ################

    ##--------------------------
    print(" ")
    if inpar["savesnapshot"]==True :
        return timearray,receiv,vxsave,vzsave
    else :
        return timearray,receiv

#####################################################################
