"""Adapts the distance computation to multiple data shapes."""
from pyinterp.geodetic import System, coordinate_distances
import numpy as np
from ocean_tools.utilities.reshape import slice_along_axis

def distances_along_axis(
    longitudes: np.ndarray,
    latitudes: np.ndarray,
    axis: int = 0,
    return_full: bool = True,
    spherical_approximation: bool = True,
    **kwargs):
    """Compute the distances point to point along a given axis.

    In case the spherical approximation is used, the great circle distance is
    computed, with the earth radius being deduced from the spheroid model
    (mean_radius).

    Else, the distance will be computed using the ellipsoid model.

    Parameters
    ----------
    longitudes
        Longitudes in degrees
    latitudes
        Latitudes in degrees
    axis
        Axis along which the distance will be computed
    return_full
        True to return an array of the same shape as the input. Useful to
        conserve chunk sizes (Default to True). The last column of the array
        will be set to nan.
    spherical_approximation
        Whether to use a spherical earth or an ellipsoid earth model
    
    Returns
    -------
    distances_along_axis
        Distance between points along the given axis in meters. Last element or
        column is set to nan if return_full is set to True.
    """
    if spherical_approximation:
        distances_along_axis = _great_circle_distance_along_axis(
            longitudes,
            latitudes,
            axis=axis,
            **kwargs)
    else:
        distances_along_axis = _spheroid_distances_along_axis(
            longitudes,
            latitudes,
            axis=axis,
            **kwargs)

    # Add a nan at the end of the array to return an array with the same block
    # size. This will be easier for dask parallelization
    if return_full:
        append_shape = list(longitudes.shape)
        append_shape[axis] = 1

        distances_along_axis = np.append(
            distances_along_axis,
            np.full(append_shape, fill_value=np.nan),
            axis=axis)
    
    return distances_along_axis


def _spheroid_distances_along_axis(
    longitudes: np.ndarray,
    latitudes: np.ndarray,
    axis: int = 0,
    wgs: System = System()):
    # Slice along axis to compute distance point to point    
    lon0 = slice_along_axis(longitudes, axis, slice(0, -1))
    lon1 = slice_along_axis(longitudes, axis, slice(1, None))
    lat0 = slice_along_axis(latitudes, axis, slice(0, -1))
    lat1 = slice_along_axis(latitudes, axis, slice(1, None))

    # Compute distance on ellipsoid
    return coordinate_distances(
        lon0.ravel(),
        lat0.ravel(),
        lon1.ravel(),
        lat1.ravel(),
        wgs=wgs
    ).reshape(lon0.shape)
    

def _great_circle_distance_along_axis(
    longitudes: np.ndarray,
    latitudes: np.ndarray,
    axis: int = 0,
    wgs: System = System()):

    longitudes = np.radians(longitudes)
    latitudes = np.radians(latitudes)

    # Mean radius is in meters
    earth_radius = wgs.mean_radius()

    delta_lat = np.abs(np.diff(latitudes, axis=axis))
    tmp = np.abs(np.diff(longitudes, axis=axis))
    delta_lon = np.minimum(tmp, 2 * np.pi - tmp)

    lat0 = slice_along_axis(latitudes, axis=axis, slice_along_axis=slice(0, -1))
    lat1 = slice_along_axis(latitudes, axis=axis, slice_along_axis=slice(1, None))

    # Haversine formulae
    tmp = (
        np.sin(delta_lat / 2.0)**2 +
        np.sin(delta_lon / 2.0)**2 * np.cos(lat0) * np.cos(lat1)
    )

    return 2.0 * np.arcsin(np.sqrt(tmp)) * earth_radius

    