"""Make function dask compatible if needed."""
from ._track_orientation import track_orientation as track_orientation_numpy
from ._distances import (
    distances_along_axis as distances_along_axis_numpy,)
import numpy as np
import dask.array as da
from typing import Union

def distances_along_axis(
    longitudes: Union[da.array, np.ndarray],
    latitudes: Union[da.array, np.ndarray],
    axis: int = 0,
    **kwargs):
    try:
        if longitudes.chunks != latitudes.chunks:
            raise Exception(
                "Latitude and longitude arrays should have the same chunks.")
        
        overlap = np.zeros(latitudes.ndim)
        overlap[axis] = 1
        overlap = tuple(overlap)

        return da.map_overlap(
            distances_along_axis_numpy,
            longitudes,
            latitudes,
            depth=overlap,
            axis=axis,
            dtype=float,
            **kwargs)

    except AttributeError:
        # Arrays are numpy arrays, no need for adding dask overlay
        return distances_along_axis_numpy(
            longitudes,
            latitudes,
            axis=axis,
            **kwargs)


def track_orientation(
    latitude: Union[da.array, np.ndarray],
    longitude: Union[da.array, np.ndarray],
    along_track_axis: int,
    half_width: int = 1,
    **kwargs):
    try:
        if latitude.chunks != longitude.chunks:
            raise Exception(
                "Latitude and longitude arrays should have the same chunks.")
        
        overlap = np.zeros(latitude.ndim)
        overlap[along_track_axis] = half_width
        overlap = tuple(overlap)

        return da.map_overlap(
            track_orientation_numpy,
            latitude,
            longitude,
            depth=overlap,
            along_track_axis=along_track_axis,
            half_width=half_width,
            **kwargs)
    except AttributeError:
        # Arrays are numpy arrays, no need for adding dask overlay
        return track_orientation_numpy(
            latitude,
            longitude,
            along_track_axis=along_track_axis,
            half_width=half_width,
            **kwargs)


# Copy docstring for the functions
track_orientation.__doc__ = track_orientation_numpy.__doc__
distances_along_axis.__doc__ = distances_along_axis_numpy.__doc__