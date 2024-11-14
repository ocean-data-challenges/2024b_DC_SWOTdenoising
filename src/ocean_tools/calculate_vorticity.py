from ocean_tools.geostrophy import strain_rate, relative_vorticity
from ocean_tools.derivatives import DerivationMethod
from ocean_tools.derivatives import directional_derivative
import xarray as xr

def speed_derivatives(
    ds: xr.Dataset,
    method,
    h,
    speed_across_name: str = "speed_across2",
    speed_along_name: str = "speed_along2",
    distances_across_name: str = "distance_across_track2",
    distances_along_name: str = "distance_along_track2",
    num_pixels: str = "num_pixels",
    num_lines: str = "num_lines",
    **kwargs) -> xr.Dataset:
    
    d_along_across=-directional_derivative(ds[speed_across_name], ds[distances_across_name], dim=num_pixels, method=method, h=h)
    d_along_along=directional_derivative(ds[speed_across_name], ds[distances_along_name], dim=num_lines, method=method, h=h)
    d_across_across=-directional_derivative(ds[speed_along_name], ds[distances_across_name], dim=num_pixels, method=method, h=h)
    d_across_along=directional_derivative(ds[speed_along_name], ds[distances_along_name], dim=num_lines, method=method, h=h)
    
    _strain_rate = strain_rate(deriv_ug_to_x = d_across_along,
                               deriv_ug_to_y = d_across_across,
                               deriv_vg_to_x = d_along_along,
                               deriv_vg_to_y = d_along_across,)

    _relative_vorticity = relative_vorticity(deriv_ug_to_y = d_across_across,
                                             deriv_vg_to_x = d_along_along,
                                             latitudes = ds.latitude,)

    ds = ds.update({"relative_vorticity" : _relative_vorticity, "strain_rate" : _strain_rate}) 
    return ds