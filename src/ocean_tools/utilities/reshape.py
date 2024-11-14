import xarray as xr
import numpy as np

def broadcast_arrays(z1: xr.DataArray, z2: xr.DataArray):
    # xarray broadcast is not able to align the chunks when broadcasting 
    dims = z1.dims
    additional_dims = [d for d in dims if d not in z2.dims]
    z2 = z2.expand_dims(additional_dims).transpose(*dims)
    return np.broadcast_arrays(z1.data, z2.data)


def slice_along_axis(array: np.ndarray, axis: int, slice_along_axis: slice):
    """Take a slice over a given axis.
    
    Similar to np.take but produce a view. The counter-part to this is that only
    slices are supported so no advanced indexing.

    Parameters
    ----------
    array
        Array to slice
    axis
        Axis along which the slice will be taken
    slice_along_axis
        Slice definition
    
    Returns
    -------
    :
        The sliced array
    """
    slices = [slice(None)] * array.ndim
    slices[axis] = slice_along_axis
    return array[tuple(slices)]
