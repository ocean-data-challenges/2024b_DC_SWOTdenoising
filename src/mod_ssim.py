import numpy as np
import xarray as xr 
import datetime
from scipy import signal
#from ocean_tools.calculate_currents import currents
#from ocean_tools.calculate_vorticity import speed_derivatives
#from ocean_tools.derivatives import DerivationMethod



def calcul_derivatives(ds_inputs, var, name = 'true'):
    ds = ds_inputs.copy()
    ds = currents(ds, var, longitude = "longitude", latitude = "latitude", h=4).load()
    ds = speed_derivatives(ds, method=DerivationMethod.STENCIL, h=4).load()

    velocity = np.sqrt(ds.speed_meridional2**2 + ds.speed_zonal2**2).values
    relative_vorticity = ds['relative_vorticity'].values
    
    ds_inputs['velocity_{}'.format(name)] = (['num_lines', 'num_pixels'], velocity)
    ds_inputs['relative_vorticity_{}'.format(name)] = (['num_lines', 'num_pixels'], relative_vorticity)
    return ds_inputs



def structural_similarity(da1,da2):
    num_lines1 = da1.num_lines.shape[0]
    num_lines2 = da2.num_lines.shape[0]
    if num_lines1 > num_lines2:
        print('WARNING')
        da1 = da1.isel(num_lines = slice(0, num_lines2))
    elif num_lines1 < num_lines2:
        print('WARNING')
        da2 = da2.isel(num_lines = slice(0, num_lines1))
        
    K1 = 0.01
    K2 = 0.03
    data_range = float((da2.max() - da2.min()).values)
    
    mean1 = np.abs(float(da1.mean().values))
    mean2 = np.abs(float(da2.mean().values))
    
    std1 = float(da1.std().values)
    std2 = float(da2.std().values)
    
    cov12 = float(xr.cov(da1, da2).values)
    
    c1 = (K1*data_range)**2
    c2 = (K2*data_range)**2
    c3 = c2/2
    
    a = (2*mean1*mean2 + c1) * (2*std1*std2 + c2) * (cov12 + c3)
    b = (mean1**2 + mean2**2 + c1) * (std1**2 + std2**2 + c2) * (std1*std2 + c3)
    
    ssim = a/b
    return ssim