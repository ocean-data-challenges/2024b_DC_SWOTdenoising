import xarray as xr
import numpy as np

from ocean_tools.utilities.reshape import broadcast_arrays
from ._dask import directional_derivative as directional_derivative_dask

def directional_derivative(
    z: xr.DataArray,
    distances_along_dim: xr.DataArray,
    dim: str,
    **kwargs) -> xr.DataArray:

    dims = z.dims
    name = z.name
    z, distances_along_dim = broadcast_arrays(z, distances_along_dim)

    dz_along_dim = directional_derivative_dask(
        z,
        distances_along_dim,
        axis=dims.index(dim),
        **kwargs)

    return xr.DataArray(
        dz_along_dim,
        dims=dims,
        name=f"deriv_{name}_along_{dim}")

    