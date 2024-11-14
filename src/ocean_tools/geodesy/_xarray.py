import xarray as xr
from ._dask import (
    distances_along_axis as distances_along_axis_dask,
    track_orientation as track_orientation_dask,)
from ._track_orientation import projection_zonal_meridional as projection_zonal_meridional_dask

from ocean_tools.utilities.reshape import broadcast_arrays


def distances_along_axis(
    longitude: xr.DataArray,
    latitude: xr.DataArray,
    dim: str,
    **kwargs,
) -> xr.DataArray:
    """Compute distances in the along and across track directions.
    
    Parameters
    ----------
    longitude
        Swath longitudes
    latitude
        Swath latitudes
    kwargs
        Kwargs for the distance_along_axis method
    """
    dims = longitude.dims

    distances = distances_along_axis_dask(
        longitude.data,
        latitude.data,
        axis=dims.index(dim),
        **kwargs)

    return xr.DataArray(distances, dims=dims)


def track_orientation(
    latitude: xr.DataArray,
    longitude: xr.DataArray,
    along_track_dim: str,
    **kwargs):

    dims = list(latitude.dims)

    angles_zonal_along = track_orientation_dask(
        latitude.data,
        longitude.data,
        along_track_axis=dims.index(along_track_dim),
        **kwargs)

    return xr.DataArray(
        angles_zonal_along,
        dims=dims,
        name="angles_zonal_along")


def projection_zonal_meridional(
    v_along: xr.DataArray,
    v_across: xr.DataArray,
    angles_zonal_along: xr.DataArray,
    **kwargs):

    _, angles_zonal_along = broadcast_arrays(v_along, angles_zonal_along)

    out = projection_zonal_meridional_dask(
        v_along.data,
        v_across.data,
        angles_zonal_along,
        **kwargs)
    
    return (
        xr.DataArray(out["v_zonal"], dims=v_along.dims, name="v_zonal"),
        xr.DataArray(out["v_meridional"], dims=v_along.dims, name="v_meridional"))