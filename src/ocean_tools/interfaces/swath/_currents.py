import numpy as np
from typing import Dict

from ocean_tools.geostrophy import geostrophic_surface_currents
from ocean_tools.derivatives import directional_derivative
from ocean_tools.geodesy import (
    projection_zonal_meridional,
    track_orientation,
    distances_along_axis)


def currents(
    ds,
    target_field: str,
    longitude: str = "longitude",
    latitude: str = "latitude",
    **kwargs,
) -> Dict[str, np.ndarray]:
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

    angles_zonal_along = track_orientation(
        latitude,
        longitude,
        along_track_dim="num_lines")

    distances_across_track = distances_along_axis(
        longitude,
        latitude,
        dim="num_pixels")

    distances_along_track = distances_along_axis(
        longitude,
        latitude,
        dim="num_lines")

    deriv_ssh_along = directional_derivative(
        ssh,
        distances_along_track,
        dim="num_lines",
        **kwargs)

    deriv_ssh_across = -directional_derivative(
        ssh,
        distances_across_track,
        dim="num_pixels",
        **kwargs)

    (
        deriv_ssh_zonal,
        deriv_ssh_meridional) = projection_zonal_meridional(
            deriv_ssh_along,
            deriv_ssh_across,
            angles_zonal_along)

    (
        speed_along,
        speed_across,) = geostrophic_surface_currents(
            deriv_ssh_along,
            deriv_ssh_across,
            ds.latitude)

    (
        speed_zonal,
        speed_meridional,) = geostrophic_surface_currents(
            deriv_ssh_zonal,
            deriv_ssh_meridional,
            ds.latitude)
    speed_zonal.attrs["short_name"] = "U"
    speed_meridional.attrs["short_name"] = "V"

    ds.update({
        "distances_along_track": distances_along_track,
        "distances_across_track": distances_across_track,
        "speed_zonal": speed_zonal,
        "speed_meridional": speed_meridional,
        "speed_along": speed_along,
        "speed_across": speed_across,
        "deriv_ssh_across": deriv_ssh_across,
        "deriv_ssh_along": deriv_ssh_along,
        "deriv_ssh_zonal": deriv_ssh_zonal,
        "angles_zonal_along": angles_zonal_along,})

    return ds
