

import numpy as np

#############################################################################

def gaussource( t, t0, f0 ) :
    """
      Derivative of a Gaussian source time function.
    """
    # boh = f0 .* (t-t0)
    # source = -8.0*boh.*exp( -boh.^2/(4.0*f0)^2 )
    a = (np.pi*f0)**2
    source = - 8.0*a*(t-t0)*np.exp( -a*(t-t0)**2 )
    
    return source

#############################################################################

def rickersource( t, t0, f0 ):
    """
      Ricker wavelet source time function.
    """
    b = (np.pi*f0*(t-t0))**2
    w = (1.0-2.0*b)*np.exp(-b)
    return w

#############################################################################

def ricker_1st_derivative_source( t, t0, f0 ):
  """
    First derivative of a wavelet source time function.
  """
  source = 2 * np.pi**2 * (t - t0) * f0**2 * np.exp(- np.pi**2 * (t - t0)**2 * f0**2) * (2 * np.pi**2 * (t - t0)**2 * f0**2 - 3)
  return source