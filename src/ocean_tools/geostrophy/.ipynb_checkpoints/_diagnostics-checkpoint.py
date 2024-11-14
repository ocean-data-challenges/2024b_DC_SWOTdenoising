import numpy as np
# from typing import Union, Dict, Optional
# import dask.array as da
from ocean_tools.utilities.constants import g, Om_t as omega_earth


def geostrophic_surface_currents(deriv_ssh_to_x: np.ndarray,
                                 deriv_ssh_to_y: np.ndarray,
                                 latitudes: np.ndarray,
                                 time_axis: int = 2,
                                 x_axis: int=0):
    """
    Compute the geostrophic surface currents, their temporal mean and the associated anomalies
    ug - ug_mean, vg - vg_mean

    Parameters
    ----------
    deriv_ssh_to_x
        Derivative of sea surface height with respect to the X axis
    deriv_ssh_to_y
        Derivative of sea surface height with respect to the Y axis
    latitudes
        Point latitudes for coriolis factor computation (degrees)
    x_axis
        Axis index for the longitudes
    time_axis
        Temporal axis

    Returns
    -------
    ug: np.ndarray
        Longitude component of the geostrophic current derived from sea surface height
    vg: np.ndarray
        Latitude component of the geostrophic current derived from sea surface height field
    """
    # TODO: this should not be necessary. Input data should be compatible
    if len(deriv_ssh_to_x.shape) - len(latitudes.shape) > 1:
        latitudes = np.expand_dims(latitudes, [time_axis, x_axis])
    elif len(deriv_ssh_to_x.shape) - len(latitudes.shape) > 0:
        latitudes = np.expand_dims(latitudes, [time_axis])

    f = coriolis_factor(latitudes)
    ug = -(g / f) * deriv_ssh_to_y
    vg = (g / f) * deriv_ssh_to_x
    
    return dict(ug=ug, vg=vg)


# def anomalies(ug: np.ndarray, vg: np.ndarray, ug_mean: np.ndarray, vg_mean: np.ndarray, time_axis: int=2):
#     """Compute the geostrophic currents anomalies and their cross-products mean.

#     Parameters
#     ----------
#     ug: np.ndarray
#         Longitude component of the geostrophic current derived from sea surface height
#     vg: np.ndarray
#         Latitude component of the geostrophic current derived from sea surface height field
#     ug_mean: np.ndarray
#         Temporal mean of the longiude component of the geostrophic current derived from sea surface height field 
#     vg_mean: np.ndarray
#         Temporal mean of the latitude component of the geostrophic current derived from sea surface height field
#     time_axis
#         Temporal axis index

#     Returns
#     -------
#     anom_u: np.ndarray
#         Anomaly of the longitude component of the geostrophic current derived from sea surface height field: ug - ug_mean
#     anom_v: np.ndarray
#         Anomaly of the latitude component of the geostrophic current derived from sea surface height field: vg - vg_mean
#     anoms_uu_mean: np.ndarray
#         Temporal mean of anom_u**2
#     anoms_uv_mean: np.ndarray
#         Temporal mean of anom_u * anom_v
#     anoms_vv_mean: np.ndarray
#         Temporal mean of anom_v ** 2
#     anoms_uu_minus_vv_mean: np.ndarray
#         Temporal mean of (anom_u**2 - anom_v**2) / 2
#     eke_mean: np.ndarray
#         Temporal mean of (anom_u**2 + anom_v**2) / 2
#     """
#     anom_u = ug - np.expand_dims(ug_mean, time_axis)
#     anom_v = vg - np.expand_dims(vg_mean, time_axis)

#     return (
#         anom_u,
#         anom_v,
#         anom_u**2,
#         anom_u * anom_v,
#         anom_v ** 2,
#         (anom_u**2 - anom_v**2) / 2,
#         (anom_u**2 + anom_v**2) / 2)


# def strain_rate(
#     deriv_ug_to_x: np.ndarray,
#     deriv_ug_to_y: np.ndarray,
#     deriv_vg_to_x: np.ndarray,
#     deriv_vg_to_y: np.ndarray) -> np.ndarray:
#     """Compute strain rate. [1/s]

#     Params
#     ------
#     deriv_ug_to_x
#         Derivative of geostrophic surface current component ug with respect to the X axis
#     deriv_vg_to_x
#         Derivative of geostrophic surface current component vg with respect to the X axis
#     deriv_ug_to_y
#         Derivative of geostrophic surface current component ug with respect to the Y axis
#     deriv_vg_to_y
#         Derivative of geostrophic surface current component vg with respect to the Y axis

#     Returns
#     -------
#     sr: np.ndarray
#         strain rate
#     """
#     shear = deriv_ug_to_y + deriv_vg_to_x
#     norm = deriv_ug_to_x - deriv_vg_to_y

#     return np.sqrt(norm * norm + shear * shear)


def relative_vorticity(
    deriv_ug_to_y: np.ndarray,
    deriv_vg_to_x: np.ndarray,
    latitudes: np.ndarray,
    time_axis: int = 2,
    x_axis: int=0) -> np.ndarray:
    """Computes normalized relative vorticity from geostrophic velocities.
    It is defined as: (dvg/dx - dug/dy) / f
    Where f is the Coriolis factor. 

    Params
    ------
    latitudes
        Point latitudes for coriolis factor computation (degrees)
    deriv_vg_to_x
        Derivative of geostrophic surface current component vg with respect to the X axis
    deriv_ug_to_y
        Derivative of geostrophic surface current component ug with respect to the Y axis
    x_axis
        Axis index for the longitudes
    time_axis
        Temporal axis
        
    Returns
    -------
    rel_vort: np.ndarray
        strain rate
    """
    
    if len(deriv_ug_to_y.shape) - len(latitudes.shape) > 1:
        latitudes = np.expand_dims(latitudes, [time_axis, x_axis])
    elif len(deriv_ug_to_y.shape) - len(latitudes.shape) > 0:
        latitudes = np.expand_dims(latitudes, [time_axis])
    
    f = coriolis_factor(latitudes)

    return ( deriv_vg_to_x - deriv_ug_to_y ) / f


# def anisotropy_eddy_variability(
#     anoms_uu_minus_vv: np.ndarray,
#     anoms_uv: np.ndarray,
#     eke: np.ndarray,
# ) -> np.ndarray:
#     """Compute the anisotropy of the eddy variability.

#     Anisotropy mean can also be computed by given the temporal mean of the
#     inputs.

#     Parameters
#     ----------
#     anoms_uv: np.ndarray
#         anom_u * anom_v
#     anoms_uu_minus_vv: np.ndarray
#         (anom_u**2 - anom_v**2) / 2
#     eke: np.ndarray
#         (anom_u**2 + anom_v**2) / 2

#     Returns
#     -------
#     anisotropy: np.ndarray
#         anisotropy of eddy variabiliy for each time step. Ratio between eddy anisotropy
#         and eddy kinetic energy
#     """
#     L = np.sqrt(anoms_uu_minus_vv**2 + anoms_uv**2)
#     return L / eke

# def eke_transfer(
#     deriv_ug_mean_to_x: np.ndarray,
#     deriv_ug_mean_to_y: np.ndarray,
#     deriv_vg_mean_to_x: np.ndarray,
#     deriv_vg_mean_to_y: np.ndarray,
#     anoms_uu_mean: np.ndarray,
#     anoms_uv_mean: np.ndarray,
#     anoms_vv_mean: np.ndarray,
# ) -> np.ndarray:
#     """
#     EKE transfert.

#     Parameters
#     ----------
#     deriv_ug_mean_to_x
#         Derivative of geostrophic surface current component ug_mean with respect to the X axis
#     deriv_vg_mean_to_x
#         Derivative of geostrophic surface current component vg_mean with respect to the X axis
#     deriv_ug_mean_to_y
#         Derivative of geostrophic surface current component ug_mean with respect to the Y axis
#     deriv_vg_mean_to_y
#         Derivative of geostrophic surface current component vg_mean with respect to the Y axis

#     anoms_uu_mean: np.ndarray
#         Temporal mean of anom_u**2
#     anoms_uv_mean: np.ndarray
#         Temporal mean of anom_u * anom_v
#     anoms_vv_mean: np.ndarray
#         Temporal mean of anom_v ** 2

#     Returns
#     -------
#     eke_transfer: np.ndarray
#         Transfert of energy between mean flow and eddies. Positive (negative) values
#         for transfert from (to) eddies to (from) mean flow
#     """    
#     return (
#         anoms_uu_mean * deriv_ug_mean_to_x +
#         anoms_uv_mean * deriv_vg_mean_to_x + 
#         anoms_uv_mean * deriv_ug_mean_to_y + 
#         anoms_vv_mean * deriv_vg_mean_to_y
#     )


def coriolis_factor(latitudes: np.ndarray, degrees=True):
    """Compute the coriolis factor.

    Parameters
    ----------
    latitudes
        Latitudes over which to compute the coriolis force
    degrees
        Whether the input is in degrees (True) or in radians (False)
        (defaults to True)
    
    Returns
    -------
    :
        The coriolis factor 2 * omega * sin(lat)
    """
    if degrees:
        latitudes = np.deg2rad(latitudes)
    
    return 2 * omega_earth * np.sin(latitudes)
