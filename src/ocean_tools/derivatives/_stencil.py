import numpy as np
import numba as nb
from sympy import finite_diff_weights
from scipy import ndimage
import logging

from ocean_tools.utilities.reshape import slice_along_axis

logger = logging.getLogger(__name__)

def stencil_derivation(
    z: np.ndarray,
    distances_along_axis: np.ndarray,
    axis: int,
    h: int = 1,
    handle_gaps: bool = True):
    """ Compute the first derivative using a stencil for the input array.

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
    gaps
        Whether to handle gaps or not by degrading the stencil order

    Returns
    =======
    deriv_z: np.ndarray
        The derivated array. It is the same size as the input array, but
        with nans on the borders
    """
    # Move axis and flatten as a 2D array
    z = np.swapaxes(z, 0, axis)
    distances_along_axis = np.swapaxes(distances_along_axis, 0, axis)

    old_shape = z.shape
    z = z.reshape(old_shape[0], -1)
    distances_along_axis = distances_along_axis.reshape(old_shape[0], -1)

    # Fix distances in the last row
    distances_along_axis = _stencil_distances(distances_along_axis)

    # Derivation
    derivative = _stencil_derivation0(z, distances_along_axis, h, handle_gaps)

    # Transform back
    derivative = derivative.reshape(old_shape)
    return np.swapaxes(derivative, 0, axis)


def _stencil_derivation0(
    z: np.ndarray,
    distances_along_axis: np.ndarray,
    h: int,
    handle_gaps: bool):
    # Produce the stencil weights (1D vector) for centered stenciu
    weights = _stencil_weights(half_size=h)

    # Get the proper shape for the stencil weights
    axes = list(np.arange(z.ndim))[1:]
    weights = np.expand_dims(weights, axis=axes)

    # Apply windowed product using correlation function
    windowed_product = ndimage.correlate(z, weights, mode="constant", cval=np.nan)

    if handle_gaps:
        windowed_product = _handle_gaps(z, windowed_product, h)

    derivative = windowed_product / distances_along_axis
    return derivative


def _stencil_distances(distances_along_axis: np.ndarray):
    # Distances along one axis should be a smooth function because it is a
    # pre-requisite for our stencil implementation. In fact, the stencil weights
    # we use suppose that the points are regularly spaced. In practice, this is
    # the case locally and we just need to change the distance locally. Just
    # filling the last row of distances by the previous valid row should be
    # sufficient.
    # This function could be improved by implementing smoothness detection to
    # check if the distances we use really are smoothed along an axis.
    if np.any(np.all(distances_along_axis[-1])):
        distances_along_axis = distances_along_axis.copy()
        distances_along_axis[-1] = distances_along_axis[-2]
    return distances_along_axis


def _handle_gaps(z, z_correlated, h):    
    # First step is to pad the input. This will allow us to look for windows
    # without having to take
    z = np.pad(z, ((h, h), (0, 0)), constant_values=np.nan)
    z_correlated = np.pad(z_correlated, ((h, h), (0, 0)), constant_values=np.nan)

    stencils = _generate_successive_stencils(h)[1:]
    for stencil in stencils:
        # Modification of z_correlated is done in-place
        _incremental_correlation(z, z_correlated, stencil)

    # Remove padding and return result
    return z_correlated[h:-h]


def _incremental_correlation(z, z_correlated, stencil):
    # Points to process. The function aims to update only the gap neighbours 
    # that have been set to nan
    z_mask_valids = ~np.isnan(z)
    z_correlated_invalids = np.isnan(z_correlated)
    mask = np.logical_and(z_mask_valids, z_correlated_invalids)
    indexes_x, indexes_y = mask.nonzero()

    # Crop decenter stencil if needed. A decentered stencil is given as an
    # symetrical array but with a 0 on one of the side (ex. [-1, 1, 0]). This
    # will make slicing a little more complex so we must define the half window
    # sizes on the left and right separately for each point
    h = stencil.size // 2
    if stencil[0] == 0:
        h_left, h_right = h - 1, h + 1
        stencil = stencil[1:]
    elif stencil[-1] == 0:
        h_left, h_right = h, h
        stencil = stencil[:-1]
    else:
        h_left, h_right = h, h + 1

    logger.debug(f"Using stencil of shape ({h_left}, {h_right - 1}) to process {len(indexes_x)} invalids")
    _correlation(z, z_correlated, stencil, h_left, h_right, indexes_x, indexes_y)


@nb.njit
def _correlation(z, z_correlated, stencil, h_left, h_right, indexes_x, indexes_y):
    # Correlation over a subset of the input array z.
    for ii, jj in zip(indexes_x, indexes_y):
        z_correlated[ii, jj] = np.sum(z[ii - h_left:ii + h_right, jj]*stencil)


def _stencil_weights(half_size=2, n=1, x0=0):
    # Compute weights for a centered stencil of a given size. The points are
    # assumed to be regularly spaced
    grid_points = np.arange(- half_size, half_size + 1).astype(int)
    return np.array(finite_diff_weights(n, grid_points, x0)[-1][-1]).astype(float) 


def _generate_successive_stencils(h: int):
    # Generate multiple centered stencils and two last decentered stencils of
    # size 2. This succession of stencil can be used to try and have a
    # derivative value where the high-order stencils do not have sufficient 
    # valid points
    stencils = []
    for half_size in np.arange(h, 0, -1):
        stencils.append(_stencil_weights(half_size))
    
    stencils.append(np.array([-1, 1, 0]))
    stencils.append(np.array([0, -1, 1]))
    
    return stencils
