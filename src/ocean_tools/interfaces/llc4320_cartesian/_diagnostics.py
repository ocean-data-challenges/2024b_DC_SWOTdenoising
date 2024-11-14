import numpy as np
import xarray as xr

from ocean_tools.geostrophy import geostrophic_surface_currents
from ocean_tools.derivatives import directional_derivative
from ocean_tools.geodesy import distances_along_axis


def currents(
    ds: xr.Dataset,
    target_field: str,
    longitude: str = "lon",
    latitude: str = "lat",
    **kwargs,
) -> xr.Dataset:
    """Compute distances in the along and across track directions.
    
    Parameters
    ----------
    longitude
        Swath longitudes
    latitude
        Swath latitudes
    kwargs
        Additional arguments for the derivation
    """
    longitude = ds[longitude]
    latitude = ds[latitude]
    ssh = ds[target_field]

    zonal_axis, meridional_axis = _find_zonal_meridional_axes(
        longitude,
        latitude)

    latitude, longitude = xr.broadcast(latitude, longitude)

    distances_zonal = distances_along_axis(
        longitude,
        latitude,
        dim=zonal_axis)

    distances_meridional = distances_along_axis(
        longitude,
        latitude,
        dim=meridional_axis)

    deriv_ssh_zonal = directional_derivative(
        ssh,
        distances_zonal,
        dim=zonal_axis,
        **kwargs)

    deriv_ssh_meridional = directional_derivative(
        ssh,
        distances_meridional,
        dim=meridional_axis,
        **kwargs)

    (
        speed_zonal,
        speed_meridional,) = geostrophic_surface_currents(
            deriv_ssh_zonal,
            deriv_ssh_meridional,
            latitude)
    speed_zonal.attrs["short_name"] = "U"
    speed_meridional.attrs["short_name"] = "V"

    return xr.Dataset({
        "distances_zonal": distances_zonal,
        "distances_meridional": distances_meridional,
        "speed_zonal": speed_zonal,
        "speed_meridional": speed_meridional,})


def speed_derivatives(
    ds: xr.Dataset,
    speed_zonal_name: str = "speed_zonal",
    speed_meridional_name: str = "speed_meridional",
    distances_zonal_name: str = "distances_zonal",
    distances_meridional_name: str = "distances_meridional",
    longitude: str = "lon",
    latitude: str = "lat",
    **kwargs) -> xr.Dataset:
    return xr.Dataset(dict(
        du_dlon=directional_derivative(ds[speed_zonal_name], ds[distances_zonal_name], dim=longitude, **kwargs),
        du_dlat=directional_derivative(ds[speed_zonal_name], ds[distances_meridional_name], dim=latitude, **kwargs),
        dv_dlon=directional_derivative(ds[speed_meridional_name], ds[distances_zonal_name], dim=longitude, **kwargs),
        dv_dlat=directional_derivative(ds[speed_meridional_name], ds[distances_meridional_name], dim=latitude, **kwargs)))


def _find_zonal_meridional_axes(longitude: xr.DataArray, latitude: xr.DataArray):
    # Retrieve meridional and zonal axes from a dataset containing data sampled
    # over a regular grid
    if len(longitude.shape) > 1 or len(latitude.shape) > 1:
        raise Exception("Expected longitudes and latitudes as 1D arrays" 
        "(cartesian grid)")

    zonal_axis = longitude.dims[0]
    meridional_axis = latitude.dims[0]

    return zonal_axis, meridional_axis
