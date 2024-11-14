import xarray as xr
import numpy as np
import xgcm
from types import MethodType
import zarr
import shutil
import os
import dask.array as da

from swot.filtering import filter_butterworth, filter_convolution2d
import swot.geostrophy as geostrophy
from ._api import MesoscaleDiagnostics
from ._xarray_adapter import xarray_adapter

@xr.register_dataset_accessor("llc")
class LLC(xgcm.Grid, MesoscaleDiagnostics):
 
    def __init__(self, xarray_obj):
        self._wrap_geostrophy_functions()
        xgcm.Grid.__init__(self, xarray_obj, periodic=False, metrics = {
            ('X',): ['dxC', 'dxG'], # X distances
            ('Y',): ['dyC', 'dyG'], # Y distances
        })
        MesoscaleDiagnostics.__init__(self, xarray_obj, dict(x_axis=["i", "i_g"], y_axis=["j", "j_g"], time_axis="time"))
    
        self.filter_convolution2d = MethodType(
            xarray_adapter("z", z_filtered=["time", "face", "j", "i"])(filter_convolution2d),
            self)

        self.filter_butterworth = MethodType(
            xarray_adapter("z", "time", z_filtered=['time','face', 'j', 'i'])(filter_butterworth),
            self)
    
    def _wrap_geostrophy_functions(self):       
        self.anomalies = MethodType(
            xarray_adapter(
                "ug", "vg", "ug_mean", "vg_mean",
                anom_u=['time','face', 'j', 'i'],
                anom_v=['time','face', 'j', 'i'],
                anoms_uu_mean=['face', 'j', 'i'],
                anoms_uv_mean=['face','j', 'i'],
                anoms_vv_mean=['face','j', 'i'],
                anoms_uu_minus_vv_mean=['face','j', 'i'],
                eke_mean=['face','j', 'i']
            )(geostrophy.anomalies),
            self)

        self.eddy_anisotropy = MethodType(
            xarray_adapter(
                "anoms_uu_minus_vv_mean", "anoms_uv_mean", "eke_mean",
                anisotropy=['face','j', 'i']
            )(geostrophy.anisotropy_eddy_variability),
            self)
        
        self.geostrophic_surface_currents = MethodType(
            xarray_adapter(
                "deriv_ssh_to_x", "deriv_ssh_to_y", "latitudes",
                ug=['time','face','j', 'i_g'],
                vg=['time','face','j_g', 'i'],
                ug_mean=['face','j', 'i_g'],
                vg_mean=['face','j_g', 'i'],
            )(geostrophy.geostrophic_surface_currents),
            self)

        self.strain_rate = MethodType(
            xarray_adapter(
                "deriv_ug_to_x", "deriv_ug_to_y", "deriv_vg_to_x", "deriv_vg_to_y",
                sr=['time','face','j', 'i']
            )(geostrophy.strain_rate),
            self)

        self.eke_transfer = MethodType(
            xarray_adapter(
                "deriv_ug_mean_to_x", "deriv_ug_mean_to_y", "deriv_vg_mean_to_x", "deriv_vg_mean_to_y",
                "anoms_uu_mean", "anoms_uv_mean", "anoms_vv_mean",
                eke_transfer=["face", "j", "i"]
            )(geostrophy.eke_transfer),
            self
        )


    def derivatives(self, variables):
        for variable in variables:
            da = self.derivative(self._obj[variable], "X")
            new_name = f"deriv_{variable}_to_x"

            self._obj = self._obj.update(
                da.to_dataset(name=new_name).reset_coords()[[new_name]]
            )

            da = self.derivative(self._obj[variable], "Y")
            new_name = f"deriv_{variable}_to_y"
            self._obj = self._obj.update(
                da.to_dataset(name=new_name).reset_coords()[[new_name]]
            )

        return self


    def coriolis(self, /, YC=None):
        """Compute coriolis force for the grid."""
        f = geostrophy.coriolis_factor(self._obj.YC) # at center points
        f_xi = geostrophy.coriolis_factor(self._obj.YG) # at vorticity points
        f_i = self.interp(f_xi, "X").chunk({"i":None}) # at v points
        f_j = self.interp(f_xi, "Y").chunk({"j":None}) # at u points
        
        self._obj = self._obj.assign_coords(f=f, f_xi=f_xi, f_i=f_i, f_j=f_j)
        return self

    def vorticity(self, ij_relabelling=True):
        zeta = (
            self.diff(self._obj.SSV * self._obj.dyC, "X") - 
            self.diff(self._obj.SSU * self._obj.dxC, "Y")
        )
        zeta = zeta / self._obj.rAz

        if ij_relabelling:
            zeta = zeta.rename({"i_g": "i", "j_g": "j"}) # artificially move to center points
        
        self._obj["vorticity"] = zeta
        return self

    def hdivergence(self):
        div = self.diff(self._obj.SSU*self._obj.dyG, "X") + self.diff(self._obj.SSV *self._obj.dxG, "Y")
        div = div / self._obj.rA

        self._obj["hdivergence"] = div
        return self


def faces_covering_area(ds, area, lon_axis="XC", lat_axis="YC"):
    """ Return faces covering the input area.

    Parameters
    ----------
    ds: xr.Dataset
        Dataset containing the LLC4320 grid
    area: Tuple[float, float, float, float]
        Area for which we want to know the face. It consists in
        a tuple with (lon_min, lat_min, lon_max, lat_max). The
        longitudes must be between [-180, 180[
    lon_axis: str
        Name of the longitude axis (default is XC)
    lat_axis: str
        Name of the latitude axis (default is YC)

    Returns
    -------
    selected_faces : dask.delayed(List[int])
        The faces that covers the area. It is a delayed object 
        and the computation can be done over a dask cluster
    """
    selected_faces = []

    for face in ds.face.data.astype(int):
        mask= np.logical_and.reduce([
            ds.XC[face].data >= area[0],
            ds.XC[face].data <= area[2],
            ds.YC[face].data >= area[1],
            ds.YC[face].data <= area[3]
        ])
        if np.sum(mask) > 0:
            selected_faces.append(face)
        
    return selected_faces
