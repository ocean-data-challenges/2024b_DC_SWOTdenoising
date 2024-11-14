from enum import Enum, auto
import numpy as np

from ._finite_differences_central import finite_differences
from ._stencil import stencil_derivation

class DerivationMethod(Enum):
    FINITE_DIFFERENCES = auto()
    STENCIL = auto()


def directional_derivative(
    z: np.ndarray,
    distances_along_axis: np.ndarray,
    axis: int,
    h: int = 1,
    method: DerivationMethod = DerivationMethod.STENCIL,
):
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
        Number of points taken for the derivation

    Returns
    =======
    deriv_z: np.ndarray
        The derivated array. It is the same size as the input array, but
        with nans on the borders
    """
    if method == DerivationMethod.STENCIL:
        return stencil_derivation(z, distances_along_axis, axis, h=h)
    elif method == DerivationMethod.FINITE_DIFFERENCES:
        return finite_differences(z, distances_along_axis, axis, h=h)
    else:
        raise Exception(f"Unknown method {method}")