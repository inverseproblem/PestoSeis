
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


##########################################################################

from .binheap import BinHeapMin 

import numpy as __NP

##########################################################################

class Grid2D():
    """
    Class containing the info of a 2D grid for traveltime calculations.
    
    Args:
      nx (int): number of points/cells in the x direction
      ny (int): number of points/cells in the y direction
      h (float): grid spacing (same for x and y)
      xinit (float): grid origin for the x direction
      yinit (float): grid origin for the y direction

    """
    def __init__(self, nx, ny, h,xinit,yinit):
        self.nx = nx
        self.ny = ny
        self.hgrid = h
        self.ntx = nx+1
        self.nty = ny+1
        self.xinit = xinit
        self.yinit = yinit
        
# from collections import namedtuple
# MyStruct = namedtuple("MyStruct", "field1 field2 field3")
# m = MyStruct(field1="foo", field2="bar", field3="baz")

##########################################################################

def ttimefmm(vel,grdh,xinit,yinit,coordsrc,coordrec,ttarrout=False) :
    """
    Traveltimes calculation given a velocity model, position of sources and receivers.
      
    Args:
       vel (ndarray): input velocity model
       grdh (float): grid spacing
       xinit (float): x axis origin coordinate
       yinit (float): y axis origin coordinate
       coordsrc (ndarray): coordinates of the sources
       coordred (ndarray): coordinates of the receivers
       ttarrout (bool): optional, returns the full array of traveltimes
      
    Returns:
       ttpicks (ndarray), ttarr (ndarray,optional): the traveltimes at the 
                            receivers and, optionally, the traveltime array(s)
      
    """

    nx,ny = vel.shape
    grd = Grid2D(nx,ny,grdh,xinit,yinit)
    ttarr = _ttFMM(vel,coordsrc,grd)
    
    nrec = coordrec.shape[0]
    ttpicks = __NP.zeros(nrec)
    for i in range(nrec) :
        ttpicks[i] = __bilinear_interp(ttarr, grd.hgrid,grd.xinit,
                                     grd.yinit,coordrec[i,0],
                                     coordrec[i,1])
    if ttarrout :
        return ttpicks,ttarr
    else :
        return ttpicks

##########################################################################

def __bilinear_interp(f,hgrid,xinit,yinit,xreq,yreq) :
    """
     Bilinear interpolation (2D).
    """
    nx,ny = f.shape
    ## rearrange such that the coordinates of corners are (0,0), (0,1), (1,0), and (1,1)
    xh=(xreq-xinit)/hgrid
    yh=(yreq-yinit)/hgrid
    i=int(__NP.floor(xh)) # indices starts from 0
    j=int(__NP.floor(yh)) # indices starts from 0
  
    ## rearrange such that the coordinates of corners are (0,0), (0,1), (1,0), and (1,1)
    xh=(xreq-xinit)/hgrid
    yh=(yreq-yinit)/hgrid
    i=int(__NP.floor(xh)) # indices starts from 0
    j=int(__NP.floor(yh)) # indices starts from 0
    ## if at the edges of domain choose previous square...
    if i==nx :
        i=i-1
    if j==ny :
        j=j-1
    if i==nx :
        i=i-1
    if j==ny :
        j=j-1

    xd=xh-(i) # indices starts from 0
    yd=yh-(j) # indices starts from 0
    intval = f[i,j]*(1.0-xd)*(1.0-yd)+f[i+1,j]*(1.0-yd)*xd + f[i,j+1]*(1.0-xd)*yd+f[i+1,j+1]*xd*yd
    return intval

##########################################################################

def __fmm_findclosestnode(x,y,xinit,yinit,h) :
    """
     Find closest grid node to given coordinates.
    """
    xzero = x-xinit
    yzero = y-yinit
    ix = int(__NP.floor(xzero/h)) 
    iy = int(__NP.floor(yzero/h)) 
    rx = xzero-(ix*h)
    ry = yzero-(iy*h)
    middle = h/2.0
    if rx>=middle :
        ix = ix+1
    if ry>=middle :
        iy = iy+1  
    ## return Int(ix+1),Int(iy+1) # julia
    return int(ix),int(iy) # python

##########################################################################

def _ttFMM(vel,src,grd) :
    """   
    Fast marching method to compute traveltimes in 2D.
    Uses a staggered grid, time array is bigger than velocity array, 
    following Podvin & Lecompte scheme.
     
    Args
      vel: input velocity model
      src: source position
      grd: grid parameters
     
    Returns:
      ttime: array of traveltimes
     
    """

    epsilon = 1e-6
      
    ## ttime
    ntx,nty=grd.ntx,grd.nty # ttime grid, size(vel)+1  ## STAGGERED GRID!!!
    nvx = grd.nx # velocity grid size
    nvy = grd.ny # velocity grid size
    inittt = 1e30
    ttime = __NP.zeros((ntx,nty))
    ttime[:,:] = inittt
    
    ## source location, etc.      
    mindistsrc = 1e-5
    onsrc = __NP.zeros((ntx,nty),dtype=bool)
    onsrc[:,:] = False
    xsrc,ysrc=src[0],src[1]

    ## grd.xinit-hgr because TIME array on STAGGERED grid
    hgr = grd.hgrid/2.0
    ix,iy = __fmm_findclosestnode(xsrc,ysrc,grd.xinit-hgr,grd.yinit-hgr,grd.hgrid)
    
    rx = src[0]-((ix)*grd.hgrid+grd.xinit-hgr) # NO (i-1)*grd.., just i*grd.. python
    ry = src[1]-((iy)*grd.hgrid+grd.yinit-hgr) # NO (i-1)*grd.., just i*grd.. python
    halfg = 0.0 #hgrid/2.0
    
    dist = __NP.sqrt(rx**2+ry**2)
    #@show dist,src,rx,ry
    if dist<=mindistsrc :
        onsrc[ix,iy] = True
        ttime[ix,iy] = 0.0 
    else :
        if (rx>=halfg) and (ry>=halfg) :
            onsrc[ix:ix+2,iy:iy+2] = True
        elif (rx<halfg) and (ry>=halfg) :
            onsrc[ix-1:ix+1,iy:iy+2] = True
        elif (rx<halfg) and (ry<halfg) :
            onsrc[ix-1:ix+1,iy-1:iy+1] = True
        elif (rx>=halfg) and (ry<halfg) :
            onsrc[ix:ix+2,iy-1:iy+1] = True
            
        ## set ttime around source ONLY FOUR points!!!
        ######=======>>>>>>>
        isrc,jsrc = __NP.where(onsrc) ##ind2sub(size(onsrc),find(onsrc))
        ######=======>>>>>>>
        for (j,i) in zip(jsrc,isrc) :
            #for i in isrc
            ## grd.xinit-hgr because TIME array on STAGGERED grid
            xp = (i)*grd.hgrid+grd.xinit-hgr # NO (i-1)*grd.., just i*grd.. python
            yp = (j)*grd.hgrid+grd.yinit-hgr
            ii = i-1   #int(__NP.floor((xsrc-grd.xinit))/grd.hgrid) #+1
            jj = j-1    #int(__NP.floor((ysrc-grd.yinit))/grd.hgrid) #+1
            #### vel[isrc[1,1],jsrc[1,1]] STAGGERED GRID!!!
            ttime[i,j] = __NP.sqrt((xsrc-xp)**2+(ysrc-yp)**2) / vel[ii-1,jj-1]#[isrc[1,1],jsrc[1,1]]
            #print("source:",ii,jj,xp,yp)

    # import matplotlib.pyplot as plt
    # plt.figure()
    # extent=
    # plt.imshow(onsrc.T,extent=extent)
    # plt.plot(xsrc/grd.hgrid,ysrc/grd.hgrid,'or')
    # plt.show()
    # exit()
    
    ###################################################### 
    ## indices of points clockwise
    ## A
    cooa = __NP.array([[-1, 0],
                     [0,  1],
                     [0,  1],
                     [1,  0],
                     [1,  0],
                     [0, -1],
                     [0, -1],
                     [-1, 0] ])
    ## B
    coob = __NP.array([[-1,  1],
                     [-1,  1],
                     [1,  1],
                     [1,   1],
                     [1, -1],
                     [1,  -1],
                     [-1, -1],
                     [-1, -1] ]) 
    ## velocity in ABCD
    coovin = __NP.array([[-1,  0],
                       [-1,  0],
                       [0,   0],
                       [0,   0],
                       [0,  -1],
                       [0,  -1],
                       [-1, -1],
                       [-1, -1] ])
    # velocity in AFGC see paper
    coovadj = __NP.array([[-1, -1],
                        [0,   0],
                        [-1,  0],
                        [0,  -1],
                        [0,   0],
                        [-1, -1],
                        [0,  -1],
                        [-1,  0] ])
    
    ##================================
    slowness = 1.0/vel
    neigh = __NP.array([[1, 0],
                      [0, 1],
                      [-1, 0],
                      [0, -1]])
                     
    #-------------------------------
    ## init FMM 
                     
    status = __NP.zeros((ntx,nty),dtype=int) #Array{Int64}(ntx,nty)
    status[:,:] = 0   ## set all to far
    status[onsrc] = 2 ## set to accepted on src
    naccinit=__NP.count_nonzero(status==2)  ## count(status.==2)
    
    ## get all i,j accepted 
    iss,jss = __NP.where(status==2) ## findn(status.==2) #ind2sub((ntx,nty),find(status.==2))
    naccinit = iss.size
    iss = __NP.zeros(naccinit,dtype=int) ##Array{Int64,1}(naccinit)
    jss = __NP.zeros(naccinit,dtype=int) ##Array{Int64,1}(naccinit)
    l=0 #1
    for j in range(grd.nty):
        for i in range(grd.ntx) :
            if status[i,j]==2 :
                iss[l] = i
                jss[l] = j
                l+=1

    ## Init the max binary heap with void arrays but max size
    Nmax = ntx*nty
    ## ???? bheap = build_minheap(Array{Float64}(0),Nmax,Array{Int64}(0))
    ## initialize an empty min bin heap
    bheap = BinHeapMin(__NP.array([]),Nmax,__NP.array([]))

    
    ## pre-allocate
    tmptt = 0.0 
    ttlocmin = __NP.zeros(8)

    ## construct initial narrow band
    for l in range(naccinit) :        
        for ne in range(4) : ## four potential neighbors

            i = iss[l] + neigh[ne,0]
            j = jss[l] + neigh[ne,1]

            ## if the point is out of bounds skip this iteration
            if (i>ntx-1) or (i<0) or (j>nty-1) or (j<0) :
                continue
            
            if status[i,j]==0 : ## far

                ## add tt of point to binary heap and give handle
                tmptt = _calcttpt(ttime,ttlocmin,inittt,slowness,grd,cooa,coob,coovin,coovadj,i,j)
                # get handle
                ##han = sub2ind((ntx,nty),i,j)
                han = __NP.ravel_multi_index((i,j),(ntx,nty))
                # insert into heap
                bheap.insert_minheap(tmptt,han)
                # change status, add to narrow band
                status[i,j]=1

    #-------------------------------
    ## main FMM loop
    totnpts = ntx*nty
    for node in range(naccinit+1,totnpts+1): ## <<<<===| CHECK !!!!

        ## if no top left exit the game...
        if bheap.Nh<0 :
            break

        han,tmptt = bheap.pop_minheap()
        ia,ja = __NP.unravel_index(han,(ntx,nty)) #ind2sub((ntx,nty),han)
        #ja = div(han,ntx) +1
        #ia = han - ntx*(ja-1)
        # set status to accepted
        status[ia,ja] = 2 # 2=accepted
        # set traveltime of the new accepted point
        ttime[ia,ja] = tmptt

        ## try all neighbors of newly accepted point
        for ne in range(4): #=1:4 

            i = ia + neigh[ne,0]
            j = ja + neigh[ne,1]
                       
            ## if the point is out of bounds skip this iteration
            if (i>ntx-1) or (i<0) or (j>nty-1) or (j<0) :
                continue

            if status[i,j]==0 : ## far, active

                ## add tt of point to binary heap and give handle
                tmptt = _calcttpt(ttime,ttlocmin,inittt,slowness,grd,cooa,coob,coovin,coovadj,i,j)
                han = __NP.ravel_multi_index((i,j),(ntx,nty)) #sub2ind((ntx,nty),i,j)
                bheap.insert_minheap(tmptt,han)
                # change status, add to narrow band
                status[i,j]=1                

            elif status[i,j]==1 : ## narrow band                

                # update the traveltime for this point
                tmptt = _calcttpt(ttime,ttlocmin,inittt,slowness,grd,cooa,coob,coovin,coovadj,i,j)
                # get handle
                han = __NP.ravel_multi_index((i,j),(ntx,nty)) #sub2ind((nx,nty),i,j)
                # update the traveltime for this point in the heap
                bheap.update_node_minheap(tmptt,han)

    ##-------------------------------
    return ttime


##########################################################################

def _calcttpt(ttime,ttlocmin,inittt,slowness,grd,cooa,coob,coovin,coovadj,i,j):
    """
    Calculate the traveltime for a given point in the grid.
    """ 

    ##--------------------------------##
    distab2 = grd.hgrid**2
    distac = grd.hgrid
    distbc = __NP.sqrt(2.0)*grd.hgrid
    distab = grd.hgrid
    distAB2divBC = distab2/distbc
    #ttlocmin = Array{Float64,1}(8)
    ttlocmin[:] = inittt
    nvx = grd.nx
    nvy = grd.ny
    
    ## if on src, skip this iteration
    # onsrc[i,j]==true && continue
    
    ####################################################
    ##   Local solver (Podvin, Lecomte, 1991)         ##
    ####################################################
    
    ttlocmin[:] = inittt  
    ttdiffr = inittt  
    ttheadw = inittt
    ttc = inittt
    
    ## loop on triangles
    for cou in range(8): #1:8
        
        ipix = coovin[cou,0] + i
        jpix = coovin[cou,1] + j
        
        ## if the point is out of bounds skip this iteration
        if (ipix>nvx-1) or (ipix<0) or (jpix>nvy-1) or (jpix<0) :
            continue
        
        # ###################################################
        # ##  transmission
        tta = ttime[i+cooa[cou,0], j+cooa[cou,1]]
        ttb = ttime[i+coob[cou,0], j+coob[cou,1]]
        testas = distAB2divBC * slowness[ipix, jpix]
        ## stability condition
        if ((tta-ttb)>=0.0) and ( (tta-ttb)<=testas ) :
            ## distab = grd.hgrid #sqrt((xa-xb)**2+(ya-yb)**2)
            ttc = tta + __NP.sqrt( (distab*slowness[ipix,jpix])**2 - (tta-ttb)**2 )
     

        ## to avoid multiple calculations of the same arrivals...
        if (cou & 1)==1 :  # check if odd ##isodd(cou)  ## cou=[1,3,5,7]
            
            ###################################################
            ## diffraction
            ttdiffr = ttb + slowness[ipix, jpix] * distbc
            
            ###################################################
            ## head wave
            ## all these cases ARE necessary!
            ## ttheadw = inittt
            iadj = coovadj[cou,0] + i
            jadj = coovadj[cou,1] + j

            if (iadj>nvx-1) or (iadj<0) or (jadj>nvy-1) or (jadj<0) :
                ##ttheadw = inittt
                ttheadw = tta + distac * slowness[ipix,jpix]
            else :
                #slowadj = slowness[i+coovadj[cou,0], j+coovadj[cou,1]]
                #1.0/vel[i+coovadj[cou,0], j+coovadj[cou,1]]
                ttheadw = tta + distac * min( slowness[ipix,jpix],
                                              slowness[iadj,jadj] )
        
        ## minimum time
        ttlocmin[cou] = min(ttc,ttheadw,ttdiffr) # ,ttlocmin)!!            
    
    ##################################################
    #ttime[i,j] = min(ttimeold[i,j],minimum(ttlocmin))
    ttf = ttlocmin.min()
    return ttf
    ##################################################

##########################################################################

