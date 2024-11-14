import numpy as np
from typing import Dict
from ocean_tools.geostrophy import geostrophic_surface_currents, relative_vorticity
from ocean_tools.derivatives import directional_derivative
from ocean_tools.geodesy import (
    projection_zonal_meridional,
    track_orientation,
    distances_along_axis, 
)


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

    # Angle de projection
    angles_zonal_along = track_orientation(
        latitude,
        longitude,
        along_track_dim="num_lines")

    #Distance along track et across track
    distances_across_track = distances_along_axis(
        longitude,
        latitude,
        dim="num_pixels")

    distances_along_track = distances_along_axis(
        longitude,
        latitude,
        dim="num_lines")

    #Dérivée dans le repère de la fauchée
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

    #Transformation dans le repère lon/lat
    (
        deriv_ssh_zonal,
        deriv_ssh_meridional) = projection_zonal_meridional(
            deriv_ssh_along,
            deriv_ssh_across,
            angles_zonal_along)

    #Vitesse dans le repère lon/lat
    (
        speed_zonal,
        speed_meridional,) = geostrophic_surface_currents(
            deriv_ssh_zonal,
            deriv_ssh_meridional,
            ds.latitude)
    speed_zonal.attrs["short_name"] = "U"
    speed_meridional.attrs["short_name"] = "V"
    
    #Vitesse dans le repère de la fauchée
    (
        speed_along,
        speed_across,) = geostrophic_surface_currents(
            deriv_ssh_along,
            deriv_ssh_across,
            ds.latitude)

    
    ds.update({
        "speed_across2" : speed_across,
        "speed_along2" : speed_along,
        "distance_across_track2" : distances_across_track,
        "distance_along_track2" : distances_along_track,
        "speed_zonal2": speed_zonal,
        "speed_meridional2": speed_meridional,
    })

    return ds