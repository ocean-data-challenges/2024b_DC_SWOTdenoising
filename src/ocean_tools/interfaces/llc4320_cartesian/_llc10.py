import xarray as xr
import dask.array as da
import numpy as np
from typing import Union, Dict
from types import MethodType
import pyinterp.geodetic

from swot.filtering import filter_butterworth, filter_convolution2d
import swot.geostrophy as geostrophy
from swot.derivatives import directional_derivative
from swot.geodesy import distances_along_axis
from ._api import MesoscaleDiagnostics
from ._xarray_adapter import xarray_adapter

DATASET_DIMS = ("lon", "lat", "time")
zonal_axis = DATASET_DIMS.index("lon")
meridional_axis = DATASET_DIMS.index("lat")
time_axis = DATASET_DIMS.index("time")

@xr.register_dataset_accessor("llc10")
class LLC10(MesoscaleDiagnostics):
    """Accessor for supporting processing of LLC4320 downsampled data.
    
    This class relies on submodules to do the processing. It handles parallisation
    and fields storing on the dataset over which we work.

    The LLC4320 downsampled data are interpolated from the native grid to a regular
    grid of 1/10°.

    Parameters
    ----------
    xarray_obj
        Dataset associated with the accessor
    """

    def __init__(self, xarray_obj):
        
        super().__init__(
            # Transpose (lon, lat, time)
            xarray_obj.transpose(*DATASET_DIMS),
            dict(x_axis="lon", y_axis="lat", time_axis="time")
        )
        
        # Wrap functions using the adapter
        self._wrap_geostrophy_functions()
        
        self.distances = MethodType(
            xarray_adapter(
                "lon", "lat",
                distances_zonal=["lon", "lat"], distances_meridional=["lon", "lat"]
            )(_distances)
            , self)

        self.derivatives = MethodType(
            xarray_adapter(
                "z", "distances_zonal", "distances_meridional",
                dz_zonal=DATASET_DIMS, dz_meridional=DATASET_DIMS
            )(_derivatives),
            self)
        
        self.derivatives_2d = MethodType(
            xarray_adapter(
                "z", "distances_zonal", "distances_meridional",
                dz_zonal=['lon', 'lat'], dz_meridional=['lon', 'lat']
            )(_derivatives),
            self)

        self.filter_convolution2d = MethodType(
            xarray_adapter("z", z_filtered=DATASET_DIMS)(filter_convolution2d),
            self)

        self.filter_butterworth = MethodType(
            xarray_adapter("z", "time", z_filtered=DATASET_DIMS)(filter_butterworth),
            self)

    def _wrap_geostrophy_functions(self):
        self.geostrophic_currents = MethodType(
            xarray_adapter(
                "deriv_ssh_to_x", "deriv_ssh_to_y", "latitudes",
                ug=DATASET_DIMS,
                vg=DATASET_DIMS, 
                ug_mean=['lon','lat'],
                vg_mean=['lon','lat'], 
            )(geostrophy.geostrophic_surface_currents),
            self)

        self.anomalies = MethodType(
            xarray_adapter(
                "ug", "vg", "ug_mean", "vg_mean",
                anom_u=DATASET_DIMS,
                anom_v=DATASET_DIMS,
                anoms_uv=DATASET_DIMS,
                anoms_uu_mean=['lon', 'lat'],
                anoms_uv_mean=['lon', 'lat'],
                anoms_vv_mean=['lon', 'lat'],
                anoms_uu_minus_vv=DATASET_DIMS,
                anoms_uu_minus_vv_mean=['lon', 'lat'],
                eke=DATASET_DIMS,
                eke_mean=['lon', 'lat']
            )(geostrophy.anomalies),
            self)


        self.eddy_anisotropy = MethodType(
            xarray_adapter(
                "anoms_uu_minus_vv_mean", "anoms_uv_mean", "eke_mean","anoms_uu_minus_vv", "anoms_uv", "eke",
                anisotropy_mean=['lon', 'lat'],
                anisotropy=DATASET_DIMS
            )(geostrophy.anisotropy_eddy_variability),
            self)

        self.strain_rate = MethodType(
            xarray_adapter(
                "deriv_ug_to_x", "deriv_ug_to_y", "deriv_vg_to_x", "deriv_vg_to_y",
                sr=DATASET_DIMS
            )(geostrophy.strain_rate),
            self)
        
        self.relative_vorticity = MethodType(
            xarray_adapter(
                "deriv_ug_to_y", "deriv_vg_to_x", "latitudes", 
                rel_vor=DATASET_DIMS
            )(geostrophy.relative_vorticity),
            self)
        
        self.eke_transfer = MethodType(
            xarray_adapter(
                "deriv_ug_mean_to_x", "deriv_ug_mean_to_y", "deriv_vg_mean_to_x", "deriv_vg_mean_to_y",
                "anoms_uu_mean", "anoms_uv_mean", "anoms_vv_mean",
                eke_transfer=["lon", "lat"]
            )(geostrophy.eke_transfer),
            self)



def _derivatives(
    z: np.ndarray,
    distances_zonal: np.ndarray,
    distances_meridional: np.ndarray,
    **kwargs) -> Dict[str, np.ndarray]:
    """Compute the derivatives in the across and along track directions.
    
    Parameters
    ----------
    z
        Variable to derive
    distances_zonal
        Distances in the zonal direction
    distances_meridional
        Distances in the meridional direction
    kwargs
        Additional arguments for the underlying call to directional_derivative
    
    See Also
    --------
    swot.geodesy.directional_derivative

    Returns
    -------
    dz_zonal: np. ndarray
        Derivative in the zonal direction
    dz_meridional: np.ndarray
        Derivative in the meridional direction
    """
    # We need to adapt the shape of the distances to z
    
    if z.ndim > 2:
        distances_meridional = np.expand_dims(distances_meridional, axis=time_axis)
        distances_zonal = np.expand_dims(distances_zonal, axis=time_axis)

    return dict(
        dz_zonal=directional_derivative(
            z,
            distances_zonal,
            axis=zonal_axis,
            **kwargs),
        dz_meridional=directional_derivative(
            z,
            distances_meridional,
            axis=meridional_axis,
            **kwargs))


def _distances(lon, lat, **kwargs):
    """
    Wraps the pyinterp.geodetic.coordinates_distances to compute the dx
    and dy distances in meters.

    This function is adapted to a regular grid, meaning longitudes and latitudes
    are 1D arrays that must be combined to produce the full mesh grid.

    Parameters
    ----------
    lon
        Longitudes in degrees, 1D vector
    lat
        Latitudes in degrees, 1D vector
    kwargs
        Kwargs for the distance_along_axis method

    Returns
    -------
    distances_zonal: np.ndarray
        Distance between points along the longitude axis in meters
    distances_meridional: np.ndarray
        Distance between points along the latitude axis in meters
    """
    lats, lons = np.meshgrid(lat, lon)
    dx = distances_along_axis(lons, lats, axis=zonal_axis, **kwargs)
    dy = distances_along_axis(lons, lats, axis=meridional_axis, **kwargs)

    return dict(distances_zonal=dx, distances_meridional=dy)