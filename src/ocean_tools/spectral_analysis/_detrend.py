from scipy.signal import detrend
from scipy import linalg
import numpy as np
from typing import Tuple
import dask.array as da

def detrend_2d(data: np.ndarray, axis: Tuple[int]=(2, 0)):
    """ 
    This function is used to detrend the data along two axes.
    The common use case is to detrend the longitudes then the time.
    
    Note that the second detrend uses the detrended data
    
    Parameters
    ----------
    data : data to detrend
    axis: tuple specifying where to find the time and longitude
    
    Returns
    -------
    :
        detrended data 
    """
    if len(axis) != 2:
        raise Exception("Expected a tuple a two axes in 'axis' parameter")
    
    data2=detrend(data, axis=axis[0]) ## detrend longitudes
    
    return detrend(data2, axis=axis[1])  ## deterend temps


def linear_detrend(data, axis=-1):
    """ Scipy detrend but returning the coefficients, and working with dask array.
    
    The function has been simplified to support only linear detrending, 
    and no breakpoints can be specified.
    
    Parameters
    ----------
    data
        Input data to detrend
    axis
        Axis along which we must detrend
        
    Returns
    -------
    :
        The detrend coefficients
    """
    # From scipy detrend
    N = data.shape[axis]
    dtype = data.dtype

    # Restructure data so that axis is along first dimension and
    #  all other dimensions are collapsed into second dimension
    ndim = data.ndim
    if axis < 0:
        axis = axis + ndim
    newdims = np.r_[axis, 0:axis, axis + 1:ndim]
    newdata = np.reshape(
        np.transpose(data, tuple(newdims)),
        (N, data.size // N)
    )

    # Find leastsq fit and remove it
    A = da.ones((N, 2), dtype)
    A[:, 0] = np.cast[dtype](np.arange(1, N + 1) * 1.0 / N)
    coef, resids, rank, s = da.linalg.lstsq(A, newdata)
    
    # Reshape data with the proper form
    old_shape = [2, *np.array(data.shape)[np.r_[0:axis, axis + 1:ndim]]]
    old_order = tuple(np.r_[1:axis+1, 0, axis+1:ndim])
    return coef.reshape(old_shape).transpose(old_order)