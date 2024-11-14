import dataclasses as dc
from heapq import merge
import numpy as np
from typing import Dict, List, Tuple, Union
import dask.array as da
from astropy.convolution import convolve, Tophat2DKernel
import scipy.ndimage
import xarray as xr
import logging

logger = logging.getLogger(__name__)

def coarse_graining(
    ug: xr.DataArray,
    vg: xr.DataArray,
    dug_dlon: xr.DataArray,
    dug_dlat: xr.DataArray,
    dvg_dlon: xr.DataArray,
    dvg_dlat: xr.DataArray,
    distances_zonal: xr.DataArray,
    distances_meridional: xr.DataArray,
    scales: Union[float, List[float]] = 150000,
    sea_water_density: float = 1027.4,
    unit_conversion_factor: float = 1e6) -> xr.Dataset:
    """Coarse graining computation over a xarray dataset.
    
    We assume the dataset already has the current u and v, their respective
    derivatives in the X and Y directions and the distances.
    
    Parameters
    ----------
    ds
        Dataset that contains the geostrophic currents, its derivatives and the
        distances along and across track.
    spatial_dims
        Spatial dimensions of the dataset. Defaults to 'lon' and 'lat'
    sea_water_density
        The sea water density in kg/m3
    unit_conversion_factor
        By default, the EKE flux is returned in [mW/m2/km], which give a unit
        conversion factor of 1e6
    variables
        Name of the variables used for 

    Parameters
    ----------
    scale
        Scale to study in meters (defaults to 150000m = 150km)

    Returns
    -------
    eke_flux:
        Energy flux at a given scale L [mW/m2/km]


    """
    if np.isscalar(scales):
        scales = [scales]

    # Prepare dataset
    logger.debug("Renaming variables with canonical names")
    ds = xr.merge([
        ug.rename("ug"),
        vg.rename("vg"),
        dug_dlon.rename("deriv_ug_to_x"),
        dug_dlat.rename("deriv_ug_to_y"),
        dvg_dlon.rename("deriv_vg_to_x"),
        dvg_dlat.rename("deriv_vg_to_y"),])
    ds = _add_variables(ds)

    # Prepare anomalies
    logger.debug("Computing mean and persisting in local process memory")
    spatial_dims = distances_zonal.dims
    ds_mean = ds.mean(dim=spatial_dims, keepdims=True)
    ds_mean = ds_mean.persist()

    step = _compute_step(distances_zonal, distances_meridional)
    eke_fluxes = []
    for scale in scales:
        # Compute kernel size
        kernel, half_size = _compute_kernel(step, scale)

        depth=tuple([
            half_size if dim in spatial_dims else 0 
            for dim in ds.dims])

        eke_flux = _eke_flux_map(ds, ds_mean, depth, kernel, sea_water_density, unit_conversion_factor)
        eke_fluxes.append(eke_flux.rename(f"eke_flux_{int(scale / 1e3)}"))
            
    return xr.merge(eke_fluxes)


def _eke_flux_map(ds, ds_mean, depth, kernel, sea_water_density, unit_conversion_factor):

    if isinstance(ds.ug.data, da.Array):
        logger.debug("Detected dask array, eke flux will be map on dask")       
        logger.debug(f"Depth for overlap: {depth}")
        shape = ds.ug.shape
        chunks = ds.ug.chunks
        
        eke_flux = da.map_overlap(
            _eke_flux,
            ds.ug.data, ds.ug2.data, ds.vg.data, ds.vg2.data, ds.ugvg.data,
            ds.deriv_ug_to_x.data, ds.deriv_ug_to_y.data, ds.deriv_vg_to_x.data, ds.deriv_vg_to_y.data,
            da.broadcast_to(ds_mean.ug.data, shape, chunks),
            da.broadcast_to(ds_mean.ug2.data, shape, chunks),
            da.broadcast_to(ds_mean.vg.data, shape, chunks),
            da.broadcast_to(ds_mean.vg2.data, shape, chunks),
            da.broadcast_to(ds_mean.ugvg.data, shape, chunks),
            da.broadcast_to(ds_mean.deriv_ug_to_x.data, shape, chunks),
            da.broadcast_to(ds_mean.deriv_ug_to_y.data, shape, chunks),
            da.broadcast_to(ds_mean.deriv_vg_to_x.data, shape, chunks),
            da.broadcast_to(ds_mean.deriv_vg_to_y.data, shape, chunks),
            depth=depth,
            dtype=float,
            kernel=kernel,
            sea_water_density=sea_water_density,
            unit_conversion_factor=unit_conversion_factor,
        )

    else:
        eke_flux = _eke_flux(
            ds.ug.data,
            ds.ug2.data,
            ds.vg.data,
            ds.vg2.data,
            ds.ugvg.data,
            ds.deriv_ug_to_x.data,
            ds.deriv_ug_to_y.data,
            ds.deriv_vg_to_x.data,
            ds.deriv_vg_to_y.data,
            ds_mean.ug.data,
            ds_mean.ug2.data,
            ds_mean.vg.data,
            ds_mean.vg2.data,
            ds_mean.ugvg.data,
            ds_mean.deriv_ug_to_x.data,
            ds_mean.deriv_ug_to_y.data,
            ds_mean.deriv_vg_to_x.data,
            ds_mean.deriv_vg_to_y.data,
            kernel=kernel,
            sea_water_density=sea_water_density,
            unit_conversion_factor=unit_conversion_factor)

    return xr.DataArray(
        eke_flux,
        dims=ds.ug.dims,
        coords=ds.ug.coords)


def _add_variables(ds):
    # Create variables of interest
    logger.debug("Adding ug*2, vg*2 and ug*vg in dataset")
    ds["ug2"] = ds.ug**2
    ds["vg2"] = ds.vg**2
    ds["ugvg"] = ds.ug * ds.vg
    return ds

def _compute_step(distances_zonal, distances_meridional):
    # Kernel computation. We need to get the step between two points in
    # meters to estimate the window size in terms of indexes
    step = np.nanmean([
        np.nanmedian(distances_zonal.values),
        np.nanmedian(distances_meridional.values)])
    logger.debug(f"Compute step for kernel: {step} m")
    return step

def _compute_kernel(step: float, scale: float) -> Tuple[np.ndarray, int]:
    kernel_size = int(scale / step)
    logger.debug(f"Kernel size deduced from step: {kernel_size}")
    kernel = Tophat2DKernel(kernel_size)
    kernel = np.expand_dims(kernel, axis=2)
    return kernel, kernel_size

def anomalies_convolutions(data, mean, kernel: np.ndarray) -> List[np.ndarray]:
    """Compute convolutions over an anomaly and restore the mean afterwards."""
    mask = np.isnan(data)
    data[mask] = 0
    conv = scipy.ndimage.correlate(
        data - mean,
        kernel,
        mode="constant",
        cval=0) + mean
    conv[mask] = np.nan
    return conv


def _eke_flux(
    ug: np.ndarray,
    ug2: np.ndarray,
    vg: np.ndarray,
    vg2: np.ndarray,
    ugvg: np.ndarray,
    deriv_ug_to_x: np.ndarray,
    deriv_ug_to_y: np.ndarray,
    deriv_vg_to_x: np.ndarray,
    deriv_vg_to_y: np.ndarray,
    ug_mean: np.ndarray,
    ug2_mean: np.ndarray,
    vg_mean: np.ndarray,
    vg2_mean: np.ndarray,
    ugvg_mean: np.ndarray,
    deriv_ug_to_x_mean: np.ndarray,
    deriv_ug_to_y_mean: np.ndarray,
    deriv_vg_to_x_mean: np.ndarray,
    deriv_vg_to_y_mean: np.ndarray,
    kernel: np.ndarray,
    sea_water_density: float,
    unit_conversion_factor: float):

    # Convolutions
    ug_conv = anomalies_convolutions(ug, ug_mean, kernel)
    ug2_conv = anomalies_convolutions(ug2, ug2_mean, kernel)
    vg_conv = anomalies_convolutions(vg, vg_mean, kernel)
    vg2_conv = anomalies_convolutions(vg2, vg2_mean, kernel)
    ugvg_conv = anomalies_convolutions(ugvg, ugvg_mean, kernel)
    deriv_ug_to_x_conv = anomalies_convolutions(deriv_ug_to_x, deriv_ug_to_x_mean, kernel)
    deriv_ug_to_y_conv = anomalies_convolutions(deriv_ug_to_y, deriv_ug_to_y_mean, kernel)
    deriv_vg_to_x_conv = anomalies_convolutions(deriv_vg_to_x, deriv_vg_to_x_mean, kernel)
    deriv_vg_to_y_conv = anomalies_convolutions(deriv_vg_to_y, deriv_vg_to_y_mean, kernel)

    # EKE flux term computation
    return - (
        (ug2_conv - ug_conv ** 2) * deriv_ug_to_x_conv +
        (ugvg_conv - ug_conv * vg_conv) * (deriv_ug_to_y_conv + deriv_vg_to_x_conv) +
        (vg2_conv - vg_conv ** 2) * deriv_vg_to_y_conv
    ) * sea_water_density * unit_conversion_factor
