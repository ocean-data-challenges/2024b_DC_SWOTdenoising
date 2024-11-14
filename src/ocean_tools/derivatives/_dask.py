from ._dispatch import directional_derivative as directional_derivative_numpy
import numpy as np
import dask.array as da


def directional_derivative(z, distances_along_axis, axis=0, h=1, **kwargs):
    try:        
        if z.chunks != distances_along_axis.chunks:
            raise Exception(
                "z and distances_along_axis arrays should have the same chunks.")
        
        overlap = np.zeros(z.ndim)
        overlap[axis] = h
        overlap = tuple(overlap)

        return da.map_overlap(
            directional_derivative_numpy,
            z,
            distances_along_axis,
            depth=overlap,
            axis=axis,
            h=h,
            dtype=float,
            **kwargs)
    except AttributeError:
        # Arrays are numpy arrays, no need for adding dask overlay
        return directional_derivative_numpy(z,
            distances_along_axis,
            axis=axis,
            h=h,
            **kwargs)



directional_derivative.__doc__ = directional_derivative_numpy.__doc__