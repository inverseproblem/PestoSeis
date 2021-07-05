"""Calculate rays in horizontally layered model

"""


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
.. module:: rayhorlayers.py
    :synopsis: Compute seismic ray paths in a horizontally layered model.
  
.. moduleauthor:: Andrea Zunino
"""

import numpy as __NP

""
def tracerayhorlay(laydep, vel, xystart,takeoffangle,maxnumiterations=20000) :
    """
    Trace rays in a horizontally layered model. 
      
    :param laydep: input depth of layers  
    :type mod: numpy.ndarray 
    :param vel: velocity for each layer 
    :type mod: numpy.ndarray 
    :param xystart: origin coordinates of the ray 
    :type mod: numpy.ndarray 
    :param takeoffangle: take off angles 
    :type mod: float 
    :param maxnumiterations: limit the number of ray segments to calculate, in case 
                             the ray never reaches the surface
    :type mod: float 

    :returns: coordinates of the the ray path, traveltime and distance covered 
    :rtype: ndarray,float,float  
   
    """

    #
    # v1,theta1
    #  
    # --------xpt1------------
    #           \
    # v2,theta2  \
    #             \  
    # ------------xpt2--------
    #        
    #  
    #   |\
    #   | \   theta
    #   |  \
    #   |-->\
    #   |    *
    #

    assert laydep.size == vel.size#+1
    assert xystart[1]>=0.0
    assert all(laydep>0.0) # first layer starts at 0 by definition

    laydep = __NP.append(0.0,__NP.asarray(laydep))
    nlay = laydep.size-1

    ids = __NP.where(laydep>=xystart[1]) ## >= !!
    ideplay = ids[0][0] #__NP.argmin(abs(pt[1]-laydep[ids]))
    raypar = __NP.sin(__NP.deg2rad(takeoffangle))/vel[ideplay]
    thetarad = __NP.deg2rad(takeoffangle)

    raycoo = __NP.zeros((1,2))
    raycoo[0,:] = __NP.array([xystart[0],xystart[1]])
    
    ##=================
    
    i=-1
    tt = 0.0
    dist = 0.0
    firstsegment=True
    
    while True :
        i+=1 # start from 0

        pt = __NP.array([raycoo[i,0],raycoo[i,1]])
        
        ## find closest layer below starting point (z)...
        ids = __NP.where(laydep>=pt[1]) ## >= !!
        ideplay = ids[0][0] #__NP.argmin(abs(pt[1]-laydep[ids]))
        #print(i,ids)
        
        if ideplay>nlay-1 :
            print(" timeray(): take off angle: {}.  Error, not enough layers!!!".format(takeoffangle))
            return raycoo,None,None
        else :
            if i==1 :
                print(" timeray(): take off angle: {}".format(takeoffangle))
        
        ##===============================
        
        if __NP.cos(thetarad)>=0.0 :
            direction="down"
        else :
            direction="up"
              
        ##===============================
        
        if direction=="down" and (not firstsegment) :
            ## arcsin domain goes [-1 1], so if
            ##   asarg>=0.0 we have a turning ray
            asarg = vel[ideplay] * raypar
            
            if abs(asarg)>=1.0 :
                # turning ray
                # get the new angle (Snell's law), vel of layer above        
                thetarad =  __NP.pi/2.0 - thetarad + __NP.pi/2.0
                zlay = laydep[ideplay-1]
                laythick = laydep[ideplay]-laydep[ideplay-1]
                vellay = vel[ideplay-1]
            else:
                # get the new angle (Snell's law), vel of layer below
                thetarad = __NP.arcsin( vel[ideplay] * raypar )
                zlay = laydep[ideplay+1]
                laythick = laydep[ideplay+1]-laydep[ideplay]
                vellay = vel[ideplay]

        elif direction=="up" and (not firstsegment) :
            # get the new angle (Snell's law), vel of layer above
            asarg = vel[ideplay-1] * raypar
            
            if abs(asarg)>=1.0:
                # turning ray
                thetarad = thetarad - __NP.pi/2.0 
                zlay = laydep[ideplay] #laydep[ideplay+1]
                laythick = laydep[ideplay]-laydep[ideplay-1] #laydep[ideplay+1]-laydep[ideplay]
                vellay = vel[ideplay] # vel[ideplay+1]
            else :           
                ## going up..
                thetarad = __NP.pi/2.0 -__NP.arcsin( vel[ideplay-1] * raypar ) + __NP.pi/2.0
                zlay = laydep[ideplay-1]
                laythick = laydep[ideplay] - laydep[ideplay-1]
                vellay = vel[ideplay-1]
                
        #############
        if i==0 :
            firstsegment=False # take off (so angle is fixed)
            if direction=="down":
                zlay = laydep[ideplay]
                vellay = vel[ideplay]
            elif direction=="up":
                zlay = laydep[ideplay]
                vellay = vel[ideplay]
                
        #####################################
    
        if abs(__NP.cos(thetarad))<=1e-12 :
            print(" timeray(): take off angle: {}. Horizontal ray.".format(__NP.rad2deg(thetarad)))
            return raycoo,None,None
        elif abs(__NP.sin(thetarad))<=1e-12 :
            xray = raycoo[i,0]
        else :
            # get the angular coefficient of ray segment
            m = __NP.cos(thetarad)/__NP.sin(thetarad)
            # find new intersection with layer at zl
            # intersection of two straight lines
            #  z = m*x +(z1-mx1)
            #  z = zlayer
            xray = (zlay + m*pt[0]-pt[1])/m
        
            
        ##-----------------------------------
        ## Do ray path calculations
        curdist = __NP.sqrt( (xray-pt[0])**2 + (zlay-pt[1])**2 )
        dist += curdist
        tt += curdist/vellay

        ##-----------------------------------
        raycoo = __NP.r_[raycoo, __NP.array([[xray,zlay]]) ]
       
        if (raycoo[-2,1]>0.0) and (abs(zlay-0.0)<=1e-6) : #direction=="up" and (abs(zlay-0.0)<=1e-6) :
            break
        elif (raycoo[-2,1]<0.0) and (abs(zlay-0.0)<=1e-6) :
            break

        elif i>maxnumiterations :
            break
        
    #############################
    #### return coordinates of the ray path, total traveltime and total length
    return raycoo,tt,dist
        

""
def _test() :
    
    import matplotlib.pyplot as PL
 
    Nlay = 5
    laydepth = __NP.linspace(0.0,2000.0,Nlay+1)[1:]
    vp = __NP.linspace(2000.0,3000.0,Nlay)
  
    xystart = __NP.array([0.0, 0.0])  ## ( x, y(depth) )
    takeoffangle = __NP.arange(0.0,90.0,10.0)

    nta = takeoffangle.size
    tt = __NP.zeros(nta)
    totdist = __NP.zeros(nta)
    xdist = __NP.zeros(nta)
    
    PL.figure(figsize=(10,8))
    PL.subplot2grid((3,4),(0,0),colspan=3,rowspan=2)
    PL.title("Rays")
    #PL.subplot(121)
    for i in range(nta):
        raycoo,tt[i],totdist[i] = tracerayhorlay(laydepth, vp, xystart,takeoffangle[i] )
        xdist[i] = raycoo[-1,0]
        #print takeoffangle[i], raycoo
        PL.plot(raycoo[:,0],raycoo[:,1],'-k',linewidth=0.6)
        PL.plot(raycoo[0,0],raycoo[0,1],'or',linewidth=0.6)
        PL.plot(raycoo[-1,0],raycoo[-1,1],'ob',linewidth=0.6)

    xmin,xmax = 0.0,xdist.max()
    for l in laydepth:
        PL.hlines(l,xmin,xmax,color='blue')

    PL.ylim([laydepth.max(),0.0])
    PL.xlabel("x [m]")
    PL.ylabel("y [m]")
    ##PL.xlimp([0.0,depthwe.max()])
    ##PL.gca().invert_yaxis()
    
    PL.subplot2grid((3,4),(0,3),rowspan=2)
    #PL.subplot(122)
    PL.title("Velocity")
    PL.step(__NP.append(vp[0],vp),__NP.append(0.0,laydepth), where='post', label='vp')
    #PL.plot(vp,depthvp,label="vp")a
    PL.gca().invert_yaxis()
    PL.legend()
    PL.ylim([laydepth.max(),0.0])
    PL.xlabel("Vp [m/s]")
    PL.ylabel("Depth [m]")

    PL.subplot2grid((3,4),(2,0),colspan=3)
    PL.title("Traveltimes at the surface")
    PL.plot(xdist,tt,'o')
    PL.xlabel("x [m]")
    PL.ylabel("Time [s]")

    PL.tight_layout()

    PL.show()


    




#############################################
# ############################################

if __name__ == "__main__" :

    _test()
