#    TTimeRays, a program to learn about seismic rays and traveltimes.
#    Copyright (C) 2019  Andrea Zunino
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
#     along with this program. If not, see <https://www.gnu.org/licenses/>.


"""
.. module:: traveltime_rays.py
    :synopsis: A set of functions to compute seismic traveltimes, trace rays and performing simple ray tomography.
 
.. moduleauthor:: Andrea Zunino
"""
import numpy as __NP
import matplotlib.pyplot as __PL
#import fdtime2d as TT
import sys as __sys
from scipy import linalg as __LA
from scipy import interpolate as __SPINT

from ttimerays.fmm2D import forwtt as _forwtt

### important!
tolerance = 1e-6

###########################################################################

def rollmod(mod,nx,ny) :
    """
    Reshape a flattened (a vector) model to 2D (an array).

    :param mod: input flattened model
    :type mod: numpy.ndarray (2D)
    :param nx,ny: sizes of the input 2D model
    :type nx,ny: int,int

    :returns: the reshaped model as a 2D array
    :rtype: numpy.ndarray (2D)
    """
    modr = mod.copy().reshape(nx,ny,order='F')
    return modr

#############################################################

def unrollmod(mod) :
    """
    Flatten a 2D model to a vector, using column-major Fortran order.
   
    :param mod: input 2D model (probably velocity)
    
    :returns: the flattened model as a vector (1D)

    """
    modu = mod.copy().flatten('F')
    return mody

#############################################################

def setupgrid(nx,ny,dh,xinit,yinit) :
    """    
    Setup grid parameters.

    :param nx,ny: grid dimensions in x and y
    :param dx,dy: grid spacing (cell size) in x and y
    :param xinit,yinit: x and y axes origin

    """
    gridpar = {}
    gridpar['dh']  = float(dh)
    gridpar['nx']  = int(nx)
    gridpar['ny']  = int(ny)
    gridpar['xinit']  = float(xinit)
    gridpar['yinit']  = float(yinit)
    gridpar['xttmin'] = gridpar['xinit'] - gridpar['dh']/2.0 #  float(xttmin)
    gridpar['yttmin'] = gridpar['yinit'] - gridpar['dh']/2.0 #  float(yttmin)
    gridpar['xttmax'] = gridpar['nx']*gridpar['dh']+gridpar['xttmin']
    gridpar['yttmax'] = gridpar['ny']*gridpar['dh']+gridpar['yttmin']
    gridpar['xvelmin'] = gridpar['xinit']
    gridpar['yvelmin'] = gridpar['yinit']
    gridpar['xvelmax'] = gridpar['nx']*gridpar['dh']+gridpar['xttmin']-gridpar['dh']/2.0
    gridpar['yvelmax'] = gridpar['ny']*gridpar['dh']+gridpar['yttmin']-gridpar['dh']/2.0
    return gridpar

###########################################################################

def lininv(G,cov_m,cov_d,mprior,dobs) :
    """
    Linear inversion under Gaussian assumptions.

    :param G: forward model matrix (d=Gm)
    :param cov_m,cov_d: covariances for model parameters and observed data, respectively
    :param mprior: prior model
    :param dobs: obsrved data

    :returns: the posterior mean model and the posterior covariance matrix

    """
    assert dobs.ndim==1
    assert mprior.ndim==1
    assert G.ndim==2
    assert cov_m.ndim==2
    assert cov_d.ndim==2
    if not type(dobs[0])==__NP.float64 :
        print("lininv(): Error: Observed data are not in a simple 1D array")
        print("lininv(): Maybe traveltime picks have not been unrolled properly? ")
        exit()
    nobs,nmodpar = G.shape[0],G.shape[1]
    postm = __NP.zeros(nmodpar)
    postC_m = __NP.zeros((nobs,nmodpar))

    print('Computing posterior mean model and covariance...')
        
    ## Computing posterior mean model
    T1 = __NP.dot(cov_m,G.T)
    A = __NP.dot(__NP.dot(G,cov_m),G.T) + cov_d  # matrix to be inverted
    #print 'A.shape:',A.shape
    b = dobs - __NP.dot(G,mprior)
    ## A^-1 b = y
    y =  __LA.solve(A,b)  # solve a lin system instead of inverting
    ##y = __LA.lstsq(A,b)[0] 
    postm = mprior + __NP.dot(T1,y)

    ## Computing posterior covariance
    ## A^-1 x2 = y2  -> x2 = A y2
    x2 = __NP.dot(G,cov_m)
    y2 = __LA.solve(A,x2)
    ##y2 = __LA.lstsq(A,x2)[0]
    postC_m = cov_m - __NP.dot(T1,x2)

    return postm,postC_m

###########################################################################

def __bilinear_interp(f,hgrid, pt):
    """ 
     Bilinear interpolation.
    """
    xreq=pt[0]
    yreq=pt[1]
    xh=xreq/hgrid
    yh=yreq/hgrid
    i=int(__NP.floor(xh))
    j=int(__NP.floor(yh))
    xd=xh-i
    yd=yh-j
    intval=f[i,j]*(1.0-xd)*(1.0-yd)+f[i+1,j]*(1.0-yd)*xd+f[i,j+1]*(1.0-xd)*yd+f[i+1,j+1]*xd*yd
    return intval

###########################################################################

def traveltime(velmod,gridpar,srcs,recs) :
    """
      Calculate traveltime for all sources and receivers.
      
      :param velmod: input velocity model
      :param gridpar: grid parameters dictionary (as defined by setupgrid())    
      :param srcs: coordinates of sources
      :param recs: coordinates of receivers
      
      :returns: traveltimes at the receivers and traveltime arrays
      :rtype: ndarray,ndarray

    """
    #assert gridpar['dx'] == gridpar['dy']
    hgrid = gridpar['dh']
    xinit = gridpar['xinit']
    yinit = gridpar['yinit']
    nsrc = srcs.shape[0]
    ttpicks = __NP.zeros(nsrc,dtype='object')
    ttime = __NP.zeros(nsrc,dtype='object')
    for i in range(nsrc) :
        line = "\rCalculating traveltime for source {} of {}".format(i+1,nsrc)
        __sys.stdout.write(line)
        __sys.stdout.flush()

        ttpicks[i],ttime[i] = _forwtt(velmod,hgrid,xinit,yinit,srcs[i,:],recs,ttarrout=True)
        
        ##ttpicks[i],ttime[i]=traveltimesinglesrc(velmod,gridpar,srcs[i,:],recs)
    print(' ')
    return ttpicks,ttime

###########################################################################

# def traveltimesinglesrc(velmod,gridpar,coordsrc,coordrec) :
#     """
#     Calculate travel time given a velocity model and souce position
#     """
#     assert gridpar['dx'] == gridpar['dy']
#     hgrid = gridpar['dx']
        
#     # eps=0.001
  
#     # ## grid in FDTIMES starts from 0.0...
#     # xsrc=coordsrc[0]-gridpar['xttmin']
#     # ysrc=coordsrc[1]-gridpar['yttmin']
#     # ## staggered grid in FDTIMES, so last col and row of vel are ignored
#     # ##  but need to be passed anyways
#     # vel2 = __NP.zeros((velmod.shape[0]+1,velmod.shape[1]+1))
#     # vel2[:-1,:-1] = velmod
#     # ## call Fortran subroutine
#     # # vel : input rank-2 array('f') with bounds (nx,ny)
#     # # xs : input float
#     # # ys : input float
#     # # h : input float
#     # # eps : input float
#     # tbuf,msg,errstatus = TT.time2dmod(vel2,xsrc,ysrc,hgrid,eps)
#     # nrec=coordrec.shape[0]
#     # ttpicks=__NP.zeros(nrec)
#     # for i in range(nrec):
#     #     ## interp in coord starting from 0.0, so subtract xyini
#     #     coorec2 = coordrec[i,:]-__NP.array([gridpar['xttmin'],gridpar['yttmin']])
#     #     ttpicks[i] = __bilinear_interp(tbuf,hgrid,coorec2)

#     ttpicks,ttarr = forwtt(vel,grdh,xinit,yinit,coordsrc,coordrec,ttarrout=False)



#     ##----------------
#     return ttpicks,tbuf


##########################################################

def __seg_intersect(a1,a2, b1,b2) :
    """
    Determine the intersection of two segments
    """
    ### http://www-cs.ccny.cuny.edu/~wolberg/capstone/intersection/Intersection%20point%20of%20two%20lines.html

    # x1 + ua (x2 - x1) = x3 + ub (x4 - x3)
    # y1 + ua (y2 - y1) = y3 + ub (y4 - y3)
    
    # If the denominator for the equations for ua and ub is 0 then the two lines are parallel.

    # If the denominator and numerator for the equations for ua and ub are 0 then the two lines are coincident.

    # The equations apply to lines, if the intersection of line segments is required then it is only necessary to test if ua and ub lie between 0 and 1. Whichever one lies within that range then the corresponding line segment contains the intersection point. If both lie within the range of 0 to 1 then the intersection point is within both line segments. 

    a1x = a1[0]
    a1y = a1[1]
    a2x = a2[0]
    a2y = a2[1]
    
    b1x = b1[0]
    b1y = b1[1]
    b2x = b2[0]
    b2y = b2[1]
    
    denom = (b2y-b1y)*(a2x-a1x) - (b2x-b1x)*(a2y-a1y)

    numer_a =  ( (b2x-b1x)*(a1y-b1y) - (b2y-b1y)*(a1x-b1x) )
    numer_b =  ( (a2x-a1x)*(a1y-b1y) - (a2y-a1y)*(a1x-b1x) )

    #print 'denom,numer_a,numer_b',denom,numer_a,numer_b
    status = 0
    if (abs(denom)<=tolerance) :
        if abs(numer_a)<=tolerance and abs(numer_b)<=tolerance :
            status = -1 #'coincident'ua,ub
        else :
            status = -2 #'parallel'
        interspoint = __NP.array([__NP.nan,__NP.nan])
        ua,ub = __NP.nan,__NP.nan
        return status,interspoint
      
    ua =  numer_a / denom
    ub =  numer_b / denom
    
    x = a1x + ua*(a2x-a1x)
    y = a1y + ua*(a2y-a1y)

    if (0.0<=ua<=1.0) and (0.0<=ub<=1.0) :
        status = 0 #'intersection'
    elif 0.0<=ua<=1.0 :
        status = 1 #'ua in limits'
    elif 0.0<=ub<=1.0 :
        status = 2 # 'ub in limits'
    else :
        status = 3 #'no segment intersection'
  
    # print '\nstatus ', status
    # print 'ua,ub',ua,ub
    # print 'x,y',x,y,a1,a2, b1,b2
        
    interspoint = __NP.array([x,y])
  
    return status,interspoint

#############################################################

def __pointisonLINE(a,b, pointc):
    """
    Check if point is on LINE
    """
    ## a, b points of line
    ## c test point
    crossproduct = (pointc[1] - a[1]) * (b[0] - a[0]) - (pointc[0] - a[0]) * (b[1] - a[1])
    if abs(crossproduct) > tolerance :
        return False   # (or != 0 if using integers)
    return True

#############################################################

def __segintersRECT(rect,seg,rectangle=True) :
    """
    Find intersection of segments and a polygon.

    """
    if rectangle==True :
        assert rect.shape[1]==2 
        assert rect.shape[0]==5
    npts = rect.shape[0]-1
    itsp = __NP.zeros((npts,2))
    itsp[:,:] = __NP.nan
    sta = __NP.zeros(npts,dtype='int') #,dtype='str')
    for i in range(npts) :
        #print '\nside number ',i,npts
        ## check that it be meaningfully ordered...
        if rectangle==True :
            assert ( rect[i,0]==rect[i+1,0] or rect[i,1]==rect[i+1,1] )
        ## skip the intersection due to the point laying on an edge
        if ( not __pointisonLINE(rect[i,:],rect[i+1,:], seg[0,:]) ) :
            sta[i],itsp[i,:] = __seg_intersect(rect[i,:],rect[i+1,:],seg[0,:],seg[1,:])
        else :
            # point is on the edge
            sta[i] = 10
    
    return sta,itsp

##########################################################

def __findclosestnode(xini,yini,dx,dy,pt) :
    # xini ?
    # yini ?
    x = pt[0]
    y = pt[1]
    ix = __NP.floor((x-xini)/dx)
    iy = __NP.floor((y-yini)/dy)
    rx = x-xini-ix*dx
    ry = y-yini-iy*dy
    if rx>=(dx/2.0) :
        ix = ix+1        
    if ry>=(dy/2.0) :
        iy = iy+1        
    xcp = ix*dx + xini
    ycp = iy*dy + yini
    #print '\n$$$$$',x,y,ix,iy,rx,ry,xcp,ycp
    if abs(ycp-pt[1])<=tolerance and abs(xcp-pt[0])<=tolerance :
        side = 'corner' 
    elif abs(ycp-pt[1])<=tolerance :
        side = 'edge_right' if xcp<pt[0] else 'edge_left'
    elif abs(xcp-pt[0])<=tolerance :
        side = 'edge_above' if ycp<pt[1] else 'edge_below'
    else :
        side1 = 'right' if xcp<pt[0] else 'left'
        side = side1+'_above' if ycp<pt[1] else side1+'_below'        
    return __NP.array([xcp,ycp]),side,__NP.array([ix,iy])

#############################################################

def __findcelledge(xini,yini,dx,dy,nx,ny, pt, endptgrad, source=False ) :

    """
     Find the corners defining the cell where ray will propagate next

    """
    
    cgridpt,pos,idxsttime = __findclosestnode(xini,yini,dx,dy,pt) 

    #print '>>>> position: ',pos, pt,cgridpt
    
    ## set lower left corner based on position with
    ##  respect to the closest grid point
    ## on the edge ;)
    if pos=='corner' :
        if endptgrad[0]>=pt[0] and endptgrad[1]>=pt[1] :
            lowerleft_corner = cgridpt
        elif endptgrad[0]<pt[0] and endptgrad[1]>=pt[1] :
            lowerleft_corner = cgridpt+__NP.array([-dx,0.0])
        elif endptgrad[0]>=pt[0] and endptgrad[1]<pt[1] :
            lowerleft_corner = cgridpt+__NP.array([0.0,-dy])
        elif endptgrad[0]<pt[0] and endptgrad[1]<pt[1] :
            lowerleft_corner = cgridpt+__NP.array([-dx,-dy])
    elif pos=='edge_above' :
        lowerleft_corner = cgridpt if endptgrad[0]>pt[0] else cgridpt+__NP.array([-dx,0.0])
    elif pos=='edge_below' :
        lowerleft_corner = cgridpt+__NP.array([0.0,-dy]) if endptgrad[0]>pt[0] else cgridpt+__NP.array([-dx,-dy])
    elif pos=='edge_left' :
        lowerleft_corner = cgridpt+__NP.array([-dx,0.0]) if endptgrad[1]>pt[1] else cgridpt+__NP.array([-dx,-dy])
    elif pos=='edge_right' :
        lowerleft_corner = cgridpt if endptgrad[1]>pt[1] else cgridpt+__NP.array([0.0,-dy])
    ## now we're inside a cell
    elif pos=='right_above' :
        lowerleft_corner = cgridpt
    elif pos=='right_below' :
        lowerleft_corner = cgridpt+__NP.array([0.0,-dy])
    elif pos=='left_above' :
        lowerleft_corner = cgridpt+__NP.array([-dx,0.0])
    elif pos=='left_below' :
        lowerleft_corner = cgridpt+__NP.array([-dx,-dy])
    else :
        print('side error')
        exit()

    rectangle = __NP.zeros((5,2))
    
    if source==False :
        
        ## build the rectangle starting from lower left corner anti-clockwise
        addx = [0.0,dx,dx,0.0,0.0]
        addy = [0.0,0.0,dy,dy,0.0]
  
    elif source==True :
        
        if pos=='corner' :
            lowerleft_corner = cgridpt+__NP.array([-dx,-dy])
            addx = [0.0,2.*dx,2.*dx,0.0,0.0]
            addy = [0.0,0.0,2.*dy,2.*dy,0.0]

        elif pos=='edge_above' :
            lowerleft_corner = cgridpt+__NP.array([-dx,0.0])
            addx = [0.0,2.*dx,2.*dx,0.0,0.0]
            addy = [0.0,0.0,dy,dy,0.0]

        elif pos=='edge_below' :
            lowerleft_corner = cgridpt+__NP.array([-dx,-dy])
            addx = [0.0,2.*dx,2.*dx,0.0,0.0]
            addy = [0.0,0.0,dy,dy,0.0]

        elif pos=='edge_left' :
            lowerleft_corner = cgridpt+__NP.array([-dx,-dy])
            addx = [0.0,dx,dx,0.0,0.0]
            addy = [0.0,0.0,2.*dy,2.*dy,0.0]

        elif pos=='edge_right' :
            lowerleft_corner = cgridpt+__NP.array([0.0,-dy])
            addx = [0.0,dx,dx,0.0,0.0]
            addy = [0.0,0.0,2.*dy,2.*dy,0.0]
        ## now we're inside a cell
        elif pos=='right_above' :
            lowerleft_corner = cgridpt
            addx = [0.0,dx,dx,0.0,0.0]
            addy = [0.0,0.0,dy,dy,0.0]
        elif pos=='right_below' :
            lowerleft_corner = cgridpt+__NP.array([0.0,-dy])
            addx = [0.0,dx,dx,0.0,0.0]
            addy = [0.0,0.0,dy,dy,0.0]
        elif pos=='left_above' :
            lowerleft_corner = cgridpt+__NP.array([-dx,0.0])
            addx = [0.0,dx,dx,0.0,0.0]
            addy = [0.0,0.0,dy,dy,0.0]
        elif pos=='left_below' :
            lowerleft_corner = cgridpt+__NP.array([-dx,-dy])
            addx = [0.0,dx,dx,0.0,0.0]
            addy = [0.0,0.0,dy,dy,0.0]
        else :
            print('side error')
            exit()

    #--------
    for i in range(5) :
        rectangle[i,:] = __NP.array([lowerleft_corner[0]+addx[i],
                                   lowerleft_corner[1]+addy[i]])
        
          
    #print 'cell edge: ',pt,rectangle,cgridpt,pos
    return rectangle,pos


##########################################################

def __ijinvelgrid(segment,gridpar) :
    """
     Get the indices of the velocity cell related to a given segment.
    """
    xini=gridpar['xvelmin']
    yini=gridpar['yvelmin']
    dx=gridpar['dh']
    dy=gridpar['dh']
    ## calculate midpoint of segment
    midpoint = (segment[0,:]+segment[1,:])/2.0
    ## indices corresponding to the nearest velocity model cell
    ivelcell = __NP.rint((midpoint[0]-xini)/dx)
    jvelcell = __NP.rint((midpoint[1]-yini)/dy)
    ijvel = __NP.array([ivelcell,jvelcell],dtype='int')
    return ijvel

##########################################################

def __globgrad(gridpar,ttime) :
    """
      Calculate interpolant functions for the gradient of the traveltime array.
    """
    #assert gridpar['dx']==gridpar['dy']
    h = gridpar['dh']
    globxgrad,globygrad = __NP.gradient(ttime,h) #,edge_order=2)
    # x, y : array_like
    # Arrays defining the data point coordinates.
    # If the points lie on a regular grid, x can
    #  specify the column coordinates and y the
    #  row coordinates, for example:
    x = __NP.linspace(gridpar['xttmin'],
                    gridpar['xttmax'],
                    gridpar['nx']+1)

    y = __NP.linspace(gridpar['yttmin'],
                    gridpar['yttmax'],
                    gridpar['ny']+1)
    ## create interpolant
    fxgrad = __SPINT.interp2d(x,y,globxgrad.T,kind='linear')
    fygrad = __SPINT.interp2d(x,y,globygrad.T,kind='linear')
    return fxgrad,fygrad

###########################################################

def __gradatpt(fxgrad,fygrad,pt) :
    """
     Compute gradient at point pt.
    """
    xgr = fxgrad(pt[0],pt[1])[0]
    ygr = fygrad(pt[0],pt[1])[0]
    #print '>> grad >>> ',xgr,ygr,'pt',pt
    return __NP.array([xgr,ygr])

##########################################################

def __nextptgrad(fxgrad,fygrad,steplen, pt ) :
    """
    Calculate target point using negative gradient of traveltime increment
    """
    grad = __gradatpt(fxgrad,fygrad, pt )
    #print 'grad at point',grad
    ## modulus of grad
    modgrad = __NP.sqrt(grad[0]**2+grad[1]**2)
    if modgrad==0.0 :
        print(modgrad==0.0)
        exit()
    #print 'mod of grad',modgrad
    ## direction OPPOSITE to grad
    #print pt,steplen,grad,modgrad
    endptgrad = pt - grad*steplen/modgrad
    return endptgrad
        
##########################################################

def __pointisonSEGMENT(seg, pt) :
    """
    Is point on a segment?
    """
    eps = 1e-6
    xminseg = seg[:,0].min()
    xmaxseg = seg[:,0].max()
    yminseg = seg[:,1].min()
    ymaxseg = seg[:,1].max()
    inx = (xminseg-eps <= pt[0] and pt[0] <= xmaxseg+eps)
    iny = (yminseg-eps <= pt[1] and pt[1] <= ymaxseg+eps)
    if inx and iny : return True
    return False

##########################################################

def __isonsrccelledge(srccell,curpt) :
    """
    Is point on a cell edge?
    """
    for i in range(4) :
        isonedge = __pointisonSEGMENT(__NP.vstack((srccell[i,:],srccell[i+1,:])), curpt)
        if isonedge :
            return True
    return False

##########################################################

def __ptinbounds(gridpar,pt) :
    if pt[0]>gridpar['xttmin'] and pt[0]<gridpar['xttmax'] and pt[1]>gridpar['yttmin'] and pt[1]<gridpar['yttmax'] :
        return True
    print('\n',pt,gridpar['xttmin'],gridpar['xttmax'],gridpar['yttmin'],gridpar['yttmax'])
    return False

##########################################################

def _traceray(gridpar,recpos,coordsrc, ttime) :
    """
    Back-trace a single ray using negative gradient of traveltime direction

    :param gridpar: grid parameters dictionary (as defined by setupgrid())    
    :param recpos: position of the receiver
    :param srcpos: position of the source
    :param ttime: traveltime array (2D)

    :returns: the traced ray 

    """
    xini = gridpar['xttmin'] 
    yini = gridpar['yttmin']
    dx = gridpar['dh'] 
    dy = gridpar['dh'] 
    nx = int(gridpar['nx']) 
    ny = int(gridpar['ny'])
    
    ## calculate gradient of traveltime  everywhere and
    ##   setup 2 interpolants
    fxgrad,fygrad = __globgrad(gridpar,ttime) 

    steplen = 1.5*__NP.sqrt(dx**2+dy**2)
    nearsource = tolerance #__NP.sqrt(dx/2.0**2+dy/2.0**2) #tolerance

    ## find cell including source
    endptgrad = __nextptgrad(fxgrad,fygrad,steplen, coordsrc )
    srccell,possrc = __findcelledge(xini,yini,dx,dy,nx,ny, coordsrc, endptgrad,source=True )

    # print 'srccell',srccell
    # __PL.plot(srccell[:,0],srccell[:,1],'o-')
    # __PL.plot(coordsrc[0],coordsrc[1],'o')
    # __PL.show()
    
    ## from receiver to first edge???
    curpt = recpos
    endptgrad = __nextptgrad(fxgrad,fygrad,steplen, curpt )
       
    rect,pos = __findcelledge(xini,yini,dx,dy,nx,ny, curpt, endptgrad )
    segment = __NP.vstack((curpt,endptgrad))
    status,itspt = __segintersRECT(rect,segment,rectangle=True)

    
    idx = __NP.where(status==0)[0][0]
    curpt = itspt[idx].copy()
    assert __ptinbounds(gridpar,curpt)
    
    ray = {}
    ray['xy'] = __NP.vstack((recpos,curpt))
    ijvel= __ijinvelgrid(ray['xy'][-2:,:],gridpar) 
    ray['ij'] = __NP.asarray(ijvel,dtype='int') #__NP.vstack((ijvel))
    seglen = __NP.linalg.norm(ray['xy'][-1,:]-ray['xy'][-2,:])
    ray['segment_len'] = __NP.asarray(seglen)
    
    cou=0
    notonsource = True
    while notonsource :
        cou+=1
        endptgrad = __nextptgrad(fxgrad,fygrad,steplen, curpt )

        rect,pos = __findcelledge(xini,yini,dx,dy,nx,ny, curpt, endptgrad ) 
        segment = __NP.vstack((curpt,endptgrad))

        status,itspt = __segintersRECT(rect,segment,rectangle=True)
        idx = __NP.where(status==0)[0][0]

        curpt = itspt[idx].copy()
        assert __ptinbounds(gridpar,curpt)
        ray['xy'] = __NP.vstack((ray['xy'],curpt))
        ijvel= __ijinvelgrid(ray['xy'][-2:,:],gridpar) 
        ray['ij'] = __NP.vstack((ray['ij'],ijvel))
        seglen = __NP.linalg.norm(ray['xy'][-1,:]-ray['xy'][-2,:])
        ray['segment_len'] = __NP.vstack((ray['segment_len'],seglen))
    
        
        # # print cou,curpt,coordsrc
            
        curptonsrccell = __isonsrccelledge(srccell,curpt)
        ## dist in case of ray 'ringing' nearby source
        disttosrc = __NP.linalg.norm(curpt-coordsrc)

        #print('\r cou {} disttosrc {}'.format(cou,disttosrc))

        # if cou>100 :
        #     print(cou,ray['xy'])
        #     __PL.title(cou)
        #     __PL.plot(ray['xy'][:,0],ray['xy'][:,1],'-')
        #     __PL.plot(segment[:,0],segment[:,1],'-')
        #     __PL.plot(segment[0,0],segment[0,1],'o-')
        #     __PL.plot(segment[1,0],segment[1,1],'^-')
        #     __PL.plot(rect[:,0],rect[:,1],'-')
        #     #__PL.plot(coordsrc[0],coordsrc[1],'o')
        #     __PL.gca().set_aspect('equal', 'datalim')
        #     __PL.show()
        # if cou>16 :
        #     exit()
        
        if curptonsrccell or disttosrc<=__NP.sqrt(2.0)*dx :# nearsource :
            ray['xy'] = __NP.vstack((ray['xy'],coordsrc))
            ijvel= __ijinvelgrid(ray['xy'][-2:,:],gridpar) 
            ray['ij'] = __NP.vstack((ray['ij'],ijvel))
            seglen = __NP.linalg.norm(ray['xy'][-1,:]-ray['xy'][-2,:])
            ray['segment_len'] = __NP.vstack((ray['segment_len'],seglen))
            notonsource = False
    #--

    return ray

##########################################################

def traceallrays(gridpar,srccoo,reccoo,grdsttime) :
    """
    Trace multiple rays, for all sources and receivers.

    :param gridpar: grid parameters dictionary (as defined by setupgrid())    
    :param srccoo: position of the source, a 2D array with two columns, representing x and y
    :param reccoo: position of the receiver, a 2D array with two columns, representing x and y
    :param grdsttime: array of traveltime arrays [list of 2D arrays, one per source]

    :returns: the traced rays 

    """
    assert reccoo.ndim==2
    assert srccoo.ndim==2
    assert grdsttime.shape[0]==srccoo.shape[0]
    nsrc = srccoo.shape[0]
    nrec = reccoo.shape[0]
    rays = __NP.zeros((nrec,nsrc),dtype='object')
    for j in range(nsrc) :
        ttime = grdsttime[j]
        singlesrc = srccoo[j,:]
        line = '\rtracing rays for source {} of {}    '.format(j+1,nsrc)
        __sys.stdout.write(line)
        __sys.stdout.flush()
        for i in range(nrec) :
            # print('tracing ray for receiver',i+1,'of',nrec)
            singlerec = reccoo[i,:]
            rays[i,j] = _traceray(gridpar,singlerec,singlesrc,ttime)
    print(' ')
    return rays

##########################################################
            
def _trace_straight_ray(gridpar,srcpos,recpos) :
    """
    Trace a straight ray, from receiver to source.

    """

    xini = gridpar['xttmin'] 
    yini = gridpar['yttmin']
    dx = gridpar['dh'] 
    dy = gridpar['dh'] 
    nx = int(gridpar['nx']) 
    ny = int(gridpar['ny'])
    
 
    steplen = 1.5*__NP.sqrt(dx**2+dy**2)
    nearsource = tolerance #__NP.sqrt(dx/2.0**2+dy/2.0**2) #tolerance

    ## find cell including source
    # endptgrad = __nextptgrad(fxgrad,fygrad,steplen, coordsrc )
    srccell,possrctxt = __findcelledge(xini,yini,dx,dy,nx,ny,srcpos,recpos,source=True )


    curpt = recpos.copy()
    endpt = srcpos.copy()
       
    rect,pos = __findcelledge(xini,yini,dx,dy,nx,ny, curpt, endpt )
    segment = __NP.vstack((curpt,endpt))
    status,itspt = __segintersRECT(rect,segment,rectangle=True)

    
    idx = __NP.where(status==0)[0][0]
    curpt = itspt[idx].copy()
    assert __ptinbounds(gridpar,curpt)
    
    ray = {}
    ray['xy'] = __NP.vstack((recpos,curpt))
    ijvel= __ijinvelgrid(ray['xy'][-2:,:],gridpar) 
    ray['ij'] = __NP.asarray(ijvel,dtype='int') #__NP.vstack((ijvel))
    seglen = __NP.linalg.norm(ray['xy'][-1,:]-ray['xy'][-2,:])
    ray['segment_len'] = __NP.asarray(seglen)
    
    cou=0
    notonsource = True
    while notonsource :
        cou+=1
        ##endpt = __nextptgrad(fxgrad,fygrad,steplen, curpt )

        rect,pos = __findcelledge(xini,yini,dx,dy,nx,ny, curpt, endpt ) 
        segment = __NP.vstack((curpt,endpt))

        status,itspt = __segintersRECT(rect,segment,rectangle=True)
        idx = __NP.where(status==0)[0][0]

        curpt = itspt[idx].copy()

        assert __ptinbounds(gridpar,curpt)
        ray['xy'] = __NP.vstack((ray['xy'],curpt))
        ijvel= __ijinvelgrid(ray['xy'][-2:,:],gridpar) 
        ray['ij'] = __NP.vstack((ray['ij'],ijvel))
        seglen = __NP.linalg.norm(ray['xy'][-1,:]-ray['xy'][-2,:])
        ray['segment_len'] = __NP.vstack((ray['segment_len'],seglen))

        curptonsrccell = __isonsrccelledge(srccell,curpt)
        ## dist in case of ray 'ringing' nearby source
        #disttorec = __NP.linalg.norm(curpt-recpos)

        # nearsource :
        if curptonsrccell : #or disttorec<=(0.5*__NP.sqrt(2.0)*dx) :
            ray['xy'] = __NP.vstack((ray['xy'],srcpos))
            ijvel= __ijinvelgrid(ray['xy'][-2:,:],gridpar) 
            ray['ij'] = __NP.vstack((ray['ij'],ijvel))
            seglen = __NP.linalg.norm(ray['xy'][-1,:]-ray['xy'][-2,:])
            ray['segment_len'] = __NP.vstack((ray['segment_len'],seglen))
            notonsource = False

    #--
    return ray

##########################################################

def traceall_straight_rays(gridpar,srccoo,reccoo) :
    """
    Trace multiple straight rays.

    :param gridpar: grid parameters dictionary (as defined by setupgrid())    
    :param srccoo: position of the source, a 2D array with two columns, representing x and y
    :param reccoo: position of the receiver, a 2D array with two columns, representing x and y  
           
    :returns: the coordinates of the ray paths
    :rtype: ndarray

    """
    assert reccoo.ndim==2
    assert srccoo.ndim==2
    nsrc = srccoo.shape[0]
    nrec = reccoo.shape[0]
    rays = __NP.zeros((nrec,nsrc),dtype='object')
    for j in range(nsrc) :
        singlesrc = srccoo[j,:]
        line = '\rtracing rays for source {} of {}    '.format(j+1,nsrc)
        __sys.stdout.write(line)
        __sys.stdout.flush()
        for i in range(nrec) :
            # print('tracing ray for receiver',i+1,'of',nrec)
            singlerec = reccoo[i,:]
            rays[i,j] = _trace_straight_ray(gridpar,singlesrc,singlerec)
            
    print(' ')
    return rays

##########################################################

def buildtomomat(gridpar,rays, ttpick) :
    """
    Build forward matrix from rays for traveltime tomography.

    :param gridpar: grid parameters dictionary (as defined by setupgrid())    
    :param rays: seismic rays, as outputted by traceallrays()
    :param ttpick: traveltime picks at the receivers

    :returns: the 'tomography' matrix and the vector of traveltime picks 
                  (flattened for performing tomography)
    :rtype: ndarray

    """
    print("Building forward matrix")
    nrec = rays.shape[0]
    nsrc = rays.shape[1]
    tomomat = __NP.zeros((nrec*nsrc,gridpar['nx']*gridpar['ny']))
    ttpickvector = __NP.zeros(nrec*nsrc)
    obscou = -1
    for isrc in range(nsrc) :
        for irec in range(nrec) :
            obscou += 1
            nseg = rays[irec,isrc]['ij'].shape[0]
            ttpickvector[obscou] = ttpick[isrc][irec]
            for s in range(nseg) :
                i,j = rays[irec,isrc]['ij'][s,:]
                idunr = __NP.ravel_multi_index((i,j),
                                             (gridpar['nx'],gridpar['ny']),
                                             order='F')
                tomomat[obscou,idunr] = rays[irec,isrc]['segment_len'][s]
            
    return tomomat,ttpickvector

###################################################################

def plotgrid(gridpar) :
    """
    Plot staggered grid edges: traveltime at nodes and velocity in cells

    :param gridpar: grid parameters dictionary (as defined by setupgrid())    

    """
    for j in range(gridpar['ny']+1):
        __PL.hlines(y = gridpar['yttmin'] + j*gridpar['dh'], xmin = gridpar['xttmin'],
                  xmax = gridpar['xttmin'] + gridpar['nx']*gridpar['dh'],
                  color='black',linewidth=0.6)
    for j in range(gridpar['nx']+1):
        __PL.vlines(x =  gridpar['xttmin'] + j*gridpar['dh'], ymin = gridpar['yttmin'],
                  ymax = gridpar['yttmin'] + gridpar['ny']*gridpar['dh'],
                  color='black',linewidth=0.6)
    __PL.xlim([gridpar['xttmin'],gridpar['xttmax']])
    __PL.ylim([gridpar['yttmax'],gridpar['yttmin']]) 
    return

####################################################################

def plotrays(src,rec,rays) :
    """
    Plot rays as polylines.

    :param src: coordinates of the sources
    :param rec: coordinates of the receivers
    :param rays: seismic rays, as outputted by traceallrays()

    """
    for j in range(rays.shape[1]) :
        __PL.plot(src[j,0],src[j,1],'ok',zorder=100)
        for i in range(rays.shape[0]) :
            __PL.plot(rays[i,j]['xy'][:,0],
                      rays[i,j]['xy'][:,1],'-',color='black',linewidth=0.6)
        for r in range(rec.shape[0]) :
            __PL.plot(rec[r,0],rec[r,1],'vk',zorder=100)
    return

####################################################################

def plotvelmod(gridpar,velmod,vmin=None,vmax=None) :
    """
    Plot velocity model as an image.

    :param gridpar: grid parameters dictionary (as defined by setupgrid())    
    :param velmod: velocity model
    :param vmin,vmax: optional values to clip the colorbar min and max values

    """
    if vmin==None and vmax==None :
        vmin=velmod.min()
        vmax=velmod.max()

    extent_vel = [gridpar['xvelmin']-gridpar['dh']/2.0,gridpar['xvelmax']+gridpar['dh']/2.0,
                  gridpar['yvelmax']+gridpar['dh']/2.0,gridpar['yvelmin']-gridpar['dh']/2.0 ]
    __PL.imshow(velmod.T,interpolation='nearest',extent=extent_vel,
              origin='upper',cmap=__PL.cm.rainbow,aspect='auto',vmin=vmin,vmax=vmax)
    cb=__PL.colorbar()
    cb.set_label('velocity')
    #__PL.gca().set_aspect('equal')
    return

####################################################################

def plotttimemod(gridpar,ttime) :
    """
    Plot traveltime array as an image.

    :param gridpar: grid parameters dictionary (as defined by setupgrid()) 
    :param ttime: traveltime array
   
    """
    extent_ttime = [gridpar['xttmin']-gridpar['dh']/2.0,gridpar['xttmax']+gridpar['dh']/2.0,
                    gridpar['yttmax']+gridpar['dh']/2.0,gridpar['yttmin']-gridpar['dh']/2.0 ]
    __PL.imshow(ttime.T,interpolation='nearest',extent=extent_ttime,
              origin='upper',cmap=__PL.cm.rainbow,aspect='auto')
    __PL.colorbar()
    #__PL.gca().set_aspect('equal')
    return


#####################################################################
#####################################################################


##########################################################
##########################################################


# if __name__ == '__main__' :
#     test()
  
