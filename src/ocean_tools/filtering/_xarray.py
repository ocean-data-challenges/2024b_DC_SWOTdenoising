import xarray as xr
from ._convolution2d import (
    Kernel2D,
    filter_convolution2d as filter_convolution2d_dask)
from ._butterworth import filter_butterworth as filter_butterworth_dask

def filter_butterworth(
    z: xr.DataArray,
    time: xr.DataArray,
    time_dim: str,
    **kwargs):

    z_filtered = filter_butterworth_dask(
        z.data, time.data, time_axis=list(z.dims).index(time_dim), **kwargs)

    return xr.DataArray(z_filtered, dims=z.dims, name="z_filtered")


def filter_convolution2d(
    z: xr.DataArray,
    kernel: Kernel2D,
    x_dim: int,
    y_dim: int,
    **kwargs):

    dims = list(z.dims)
    z_filtered = filter_convolution2d_dask(
        kernel,
        z.data,
        x_axis=dims.index(x_dim),
        y_axis=dims.index(y_dim),
        **kwargs)

    return xr.DataArray(z_filtered, dims=dims, name="z_filtered")