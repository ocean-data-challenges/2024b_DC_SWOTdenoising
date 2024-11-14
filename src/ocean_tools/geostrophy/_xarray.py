import xarray as xr

from ._diagnostics import (
    #anomalies as _anomalies,
    geostrophic_surface_currents as _geostrophic_surface_currents,
    #anisotropy_eddy_variability as _anisotropy,
    #eke_transfer as _eke_transfer,
    strain_rate as _strain_rate,
    relative_vorticity as _relative_vorticity
    )

from ocean_tools.utilities.reshape import broadcast_arrays
from ocean_tools.utilities.decorators import stamp_processing_info

"""
@stamp_processing_info
def anomalies(
    ug: xr.DataArray,
    vg: xr.DataArray,
    ug_mean: xr.DataArray,
    vg_mean: xr.DataArray,
    time_dim: str) -> xr.Dataset:
    
    dims = list(ug.dims)
    dims_mean = list(ug_mean.dims)

    out = _anomalies(
        ug.data,
        vg.data,
        ug_mean.data,
        vg_mean.data,
        time_axis=dims.index(time_dim)
    )

    return xr.Dataset(dict(
            anom_u=xr.DataArray(out[0], dims=dims),
            anom_v=xr.DataArray(out[1], dims=dims),
            anoms_uu=xr.DataArray(out[2], dims=dims),
            anoms_uv=xr.DataArray(out[3], dims=dims),
            anoms_vv=xr.DataArray(out[4], dims=dims),
            anoms_uu_minus_vv=xr.DataArray(out[5], dims=dims),
            anoms_uu_plus_vv=xr.DataArray(out[6], dims=dims) ))
"""
@stamp_processing_info
def geostrophic_surface_currents(
    deriv_ssh_to_x: xr.DataArray,
    deriv_ssh_to_y: xr.DataArray,
    latitudes: xr.DataArray):
    
    dims = list(deriv_ssh_to_x.dims)
    _, latitudes = broadcast_arrays(deriv_ssh_to_x, latitudes)

    # x_axis is set to None because our arrays are already broadcast to the same
    # shape 
    out = _geostrophic_surface_currents(
        deriv_ssh_to_x.data,
        deriv_ssh_to_y.data,
        latitudes,
        x_axis=None)

    return (
        xr.DataArray(out["ug"], dims=dims, name="speed_x", attrs=dict(units="m/s")),
        xr.DataArray(out["vg"], dims=dims, name="speed_y", attrs=dict(units="m/s")),)


@stamp_processing_info
def strain_rate(
    deriv_ug_to_x: xr.DataArray,
    deriv_ug_to_y: xr.DataArray,
    deriv_vg_to_x: xr.DataArray,
    deriv_vg_to_y: xr.DataArray,) -> xr.DataArray:

    dims = list(deriv_ug_to_x.dims)

    # x_axis is set to None because our arrays are already broadcast to the same
    # shape 
    out = _strain_rate(
        deriv_ug_to_x.data,
        deriv_ug_to_y.data,
        deriv_vg_to_x.data,
        deriv_vg_to_y.data)

    return xr.DataArray(out, dims=dims, name="strain_rate", attrs=dict(units="s-1"))


@stamp_processing_info
def relative_vorticity(
    deriv_ug_to_y: xr.DataArray,
    deriv_vg_to_x: xr.DataArray,
    latitudes: xr.DataArray) -> xr.DataArray:
    
    dims = list(deriv_ug_to_y.dims)
    _, latitudes = broadcast_arrays(deriv_ug_to_y, latitudes)

    out = _relative_vorticity(
        deriv_ug_to_y.data,
        deriv_vg_to_x.data,
        latitudes)
    
    return xr.DataArray(out, dims=dims, name="vorticity")


# @stamp_processing_info
# def eke_transfer(
#     deriv_ug_to_x: xr.DataArray,
#     deriv_ug_to_y: xr.DataArray,
#     deriv_vg_to_x: xr.DataArray,
#     deriv_vg_to_y: xr.DataArray,
#     anoms_uu: xr.DataArray,
#     anoms_uv: xr.DataArray,
#     anoms_vv: xr.DataArray) -> xr.DataArray:
    
#     out = _eke_transfer(
#         deriv_ug_to_x.data,
#         deriv_ug_to_y.data,
#         deriv_vg_to_x.data,
#         deriv_vg_to_y.data,
#         anoms_uu.data,
#         anoms_uv.data,
#         anoms_vv.data)
    
#     return xr.DataArray(out, dims=deriv_ug_to_x.dims, name="eke_transfer")


# @stamp_processing_info
# def anisotropy(
#     anoms_uu_minus_vv: xr.DataArray,
#     anoms_uv: xr.DataArray,
#     eke: xr.DataArray) -> xr.DataArray:
    
#     dims = list(eke.dims)

#     out = _anisotropy(
#         anoms_uu_minus_vv.data,
#         anoms_uv.data,
#         eke.data)
    
#     return xr.DataArray(out, dims=dims, name="anisotropy")
