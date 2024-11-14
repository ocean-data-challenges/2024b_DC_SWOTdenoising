import numpy as np
import dask.array as da
from typing import Union
from ocean_tools.utilities.reshape import slice_along_axis


def moving_window_cumsum(distances_along_axis, axis, h):
    """Cumsum over a moving window."""
    cumsum = np.cumsum(distances_along_axis, axis=axis)
    
    cumsum_after = slice_along_axis(
        cumsum,
        axis=axis,
        slice_along_axis=slice(2*h, None))
    
    cumsum_before = slice_along_axis(
        cumsum,
        axis=axis,
        slice_along_axis=slice(0, -2*h))
    
    moving_cumsum = cumsum_after - cumsum_before
    return moving_cumsum


def finite_differences(
    z: np.ndarray,
    distances_along_axis: np.ndarray,
    axis: int,
    h: int = 1):
    """ Compute the finite difference for the input array.

    The computation uses the central derivation scheme:
    (z[ii + h] - z[ii - h]) / 2*step

    Parameters
    ==========
    z
        Input vector to derive of shape (M, N, ...)
    distances_along_axis
        Distance between consecutive point along the given axis. The array is
        also of size (M, N, ...) but with nans for the last column. distances[0]
        should give the distances between point[0] and point[1]
    axis
        Direction of the derivation
    h
        Size in index to take for the central scheme. For ex., if h=1, we will
        take one point prior and one point later for each derivation.

    Returns
    =======
    deriv_z: np.ndarray
        The derivated array. It is the same size as the input array, but
        with nans on the borders
    """
    # Take slice over the dynamically given axis
    z_after = slice_along_axis(z, axis=axis, slice_along_axis=slice(2 * h, None))
    z_before = slice_along_axis(z, axis=axis, slice_along_axis=slice(0, - 2 * h))

    # Sum of distances in a rolling window
    delta_distances = moving_window_cumsum(distances_along_axis, axis=axis, h=h)

    # Compute the finite difference central scheme
    deriv_z = (z_after - z_before) / delta_distances

    # Fill with nans
    insert_shape = list(z.shape)
    insert_shape[axis] = h
    insertion = np.full(insert_shape, fill_value=np.nan)
    return np.concatenate([insertion, deriv_z, insertion], axis=axis)
