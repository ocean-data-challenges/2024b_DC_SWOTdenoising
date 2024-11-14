import dataclasses as dc
import numpy as np
from typing import Dict, List, Tuple
import dask.array as da
from astropy.convolution import convolve, Tophat2DKernel
import xarray as xr

@dc.dataclass
class CoarseGraining:
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
    """
    ds: xr.Dataset
    spatial_dims: Tuple[str, str] = ("lon", "lat")
    sea_water_density: float = 1027.4
    unit_conversion_factor: float = 1e6
    variables: List = dc.field(default_factory=lambda: {
        "ug": "Geostrophic current speed in zonal direction (m/s)",
        "vg": "Geostrophic current speed in meridional direction (m/s)",
        "deriv_ug_to_x": "Derivative of the zonal geostrophic current speed (U) in the zonal direction",
        "deriv_ug_to_y": "Derivative of the zonal geostrophic current speed (U) in the meridional direction",
        "deriv_vg_to_x": "Derivative of the zonal geostrophic current speed (U) in the zonal direction",
        "deriv_vg_to_y": "Derivative of the zonal geostrophic current speed (U) in the meridional direction",
        "distances_zonal": "Distances between two points in the zonal direction (m)",
        "distances_meridional": "Distances between two points in the meridional direction (m)",})
    
    def __post_init__(self):
        # Check we have all the fields we need for coarse graining.
        check_all_fields_in_dataset(self.ds, self.variables)

        # Don't keep the unused variables. For now this step is optional as the
        # variable loading should be done in a lazy way
        self.ds = self.ds[self.variables.keys()]

        # We will need the median distance between two points in order to
        # estimate the kernel size in indices
        self._compute_step()

        # Add variables needed for coarse graining computation
        self._add_variables()
       
    def _compute_step(self):
        # Kernel computation. We need to get the step between two points in
        # meters to estimate the window size in terms of indexes
        self.step = np.nanmedian([
            np.nanmedian(self.ds.distances_zonal.values),
            np.nanmedian(self.ds.distances_meridional.values)
        ])

    def _compute_kernel(self, scale: float) -> Tuple[np.ndarray, int]:
        kernel_size = int(scale / self.step)
        kernel = Tophat2DKernel(kernel_size)
        kernel = np.expand_dims(kernel, axis=2)
        return kernel, kernel_size
    
    def _add_variables(self):
        # Create variables of interest
        self.ds["ug2"] = self.ds.ug**2
        self.ds["vg2"] = self.ds.vg**2
        self.ds["ugvg"] = self.ds.ug * self.ds.vg

    def run(self, scale: float = 150000) -> xr.DataArray:
        """Compute the EKE flux using the coarse graining method.

        Parameters
        ----------
        scale
            Scale to study in meters (defaults to 150000m = 150km)

        Returns
        -------
        eke_flux:
            Energy flux at a given scale L [mW/m2/km]
        """
        ds = self.ds
        ds_mean = ds.mean(dim=self.spatial_dims, keepdims=True)

        kernel, half_size = self._compute_kernel(scale)

        depth=tuple([
            half_size if dim in self.spatial_dims else 0 
            for dim in ds.dims])
        
        if isinstance(ds.ug.data, da.Array):
            # ds_mean = ds_mean.persist()
            
            shape = ds.ug.shape
            chunks = ds.ug.chunks
            
            eke_flux = da.map_overlap(
                self._eke_flux,
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
            )
        else:
            eke_flux = self._eke_flux(
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
                kernel=kernel)

        return xr.DataArray(
            eke_flux,
            dims=ds.ug.dims,
            coords=ds.ug.coords)

    def _eke_flux(self,
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
        kernel: np.ndarray):

        # Convolutions
        convolutions = anomalies_convolutions(
            (ug, ug_mean),
            (ug2, ug2_mean),
            (vg, vg_mean),
            (vg2, vg2_mean),
            (ugvg, ugvg_mean),
            (deriv_ug_to_x, deriv_ug_to_x_mean),
            (deriv_ug_to_y, deriv_ug_to_y_mean),
            (deriv_vg_to_x, deriv_vg_to_x_mean),
            (deriv_vg_to_y, deriv_vg_to_y_mean),
            kernel=kernel)

        # Eke flux computation
        return eke_flux(
            *convolutions,
            sea_water_density=self.sea_water_density,
            unit_conversion_factor=self.unit_conversion_factor)


def check_all_fields_in_dataset(ds: xr.Dataset, variables: Dict[str, str]):
    """Check that all expected fields are present in the dataset.
    
    For more clarity, a description of the expected variable is given in the
    error message.
    
    Parameters
    ----------
    ds
        The dataset we need to check
    fields
        Dictionnary
    
    Raises
    ------
    Exception
        In case there is a missing field in the dataset
    """
    try:
        for field, description in variables.items():
            ds[field]
    except KeyError:
        raise Exception(f"""
        Dataset misses the field name='{field}', description='{description}'. If
        this variable exists in your dataset, rename it with its expected name
        """)


def anomalies_convolutions(
    *args: List[Tuple[np.ndarray, np.ndarray]],
    kernel: np.ndarray) -> List[np.ndarray]:
    """Compute convolutions over an anomaly.
    
    This method groups the convolutions over multiple variables anomalies. It
    takes the variable data and removes its mean, compute the convolution and
    add the mean to the result.
    """
    convolutions = []
    for data, mean in args:
        convolutions.append(
            convolve(
                data - mean,
                kernel=kernel,
                boundary='extend', 
                fill_value = 0,
                nan_treatment='fill',
                preserve_nan = True
            ) + mean)
    return convolutions


def eke_flux(
    ug_conv: np.ndarray,
    ug2_conv: np.ndarray,
    vg_conv: np.ndarray,
    vg2_conv: np.ndarray,
    ugvg_conv: np.ndarray,
    deriv_ug_to_x_conv: np.ndarray,
    deriv_ug_to_y_conv: np.ndarray,
    deriv_vg_to_x_conv: np.ndarray,
    deriv_vg_to_y_conv: np.ndarray,
    sea_water_density: float,
    unit_conversion_factor: float) -> np.ndarray:
    """Compute the EKE flux from the convolutions results."""
    return - (
        (ug2_conv - ug_conv ** 2) * deriv_ug_to_x_conv +
        (ugvg_conv - ug_conv * vg_conv) * (deriv_ug_to_y_conv + deriv_vg_to_x_conv) +
        (vg2_conv - vg_conv ** 2) * deriv_vg_to_y_conv
    ) * sea_water_density * unit_conversion_factor
