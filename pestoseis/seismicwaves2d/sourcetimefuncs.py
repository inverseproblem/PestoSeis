

import numpy as NP

#############################################################################

def gaussource( t, t0, f0 ) :
    """
      Derivative of a Gaussian source time function.
    """
    # boh = f0 .* (t-t0)
    # source = -8.0*boh.*exp( -boh.^2/(4.0*f0)^2 )
    a = (NP.pi*f0)**2
    source = - 8.0*a*(t-t0)*NP.exp( -a*(t-t0)**2 )
    
    return source

#############################################################################

def rickersource( t, t0, f0 ):
    """
      Ricker wavelet source time function.
    """
    b = (NP.pi*f0*(t-t0))**2
    w = (1.0-2.0*b)*NP.exp(-b)
    return w

#############################################################################
