""" 
mod_switchvar.py - Ocean Variable Conversion Module

This module provides functions for converting ocean variables:
 
ssh_to_currents(ssh, lon, lat)
    Compute ocean current components from sea surface height (SSH).

currents_to_potential_vorticity(u, v, lon, lat, h)
    Compute oceanic potential vorticity (PV) from ocean current components.

currents_to_relative_vorticity(u, v, lon, lat)
    Compute oceanic relative vorticity (zeta) from ocean current components.

lonlat2dxdy(lon, lat)
    Compute spatial grid spacings (dx, dy) from longitude and latitude arrays.
    
"""

import numpy as np 
import xarray as xr
from src.mod_compare import *
from src.mod_interp import *

  

def ssh_to_currents(ssh, lon, lat):
    """
    Compute ocean current components from ssh.

    Parameters
    ----------
    u : xarray.dataset
        3D array representing sea surface height (time, x, y). 
    lon : numpy.ndarray
        1D or 2D array of longitudes corresponding to the x-dimension of u and v.
    lat : numpy.ndarray
        1D or 2D array of latitudes corresponding to the y-dimension of u and v.

    Returns
    -------
    u : numpy.ndarray
        3D array of oceanic zonal current component (time, x, y).
    v : numpy.ndarray
        3D array of oceanic meridional current component (time, x, y).
    """
    
    if len(np.shape(lon))==1: 
        lon, lat = np.meshgrid(lon,lat)
    
    # Compute geostrophic current velocities
    g = 9.81  
    ff = 4*np.pi/86164*np.sin(lat*np.pi/180)
    DX,DY = lonlat2dxdy(lon,lat)
    
    dlon = np.gradient(lon, axis=1)
    dlat = np.gradient(lat, axis=0)
    
    u = ssh.differentiate('latitude') 
    u.data *= -g/ff/DY
    v = ssh.differentiate('longitude')
    v.data *= g/ff/DX 
    uu = np.array( u.data * dlat )
    vv = np.array( v.data * dlon )

    return uu, vv


def currents_to_potential_vorticity(u, v, lon, lat, h):
    """
    Compute oceanic potential vorticity (PV) from ocean current components.

    Parameters
    ----------
    u : numpy.ndarray
        3D array representing eastward ocean current velocity (time, x, y).
    v : numpy.ndarray
        3D array representing northward ocean current velocity (time, x, y).
    lon : numpy.ndarray
        1D or 2D array of longitudes corresponding to the x-dimension of u and v.
    lat : numpy.ndarray
        1D or 2D array of latitudes corresponding to the y-dimension of u and v. 
    h : numpy.ndarray
        2D array of water depth.

    Returns
    -------
    pv : numpy.ndarray
        3D array of oceanic potential vorticity (time, x, y).
    """
    
    
    if len(np.shape(lon))==1: 
        lon, lat = np.meshgrid(lon,lat)
        
    # Calculate the Coriolis parameter (f) using the latitude
    f = 2 * np.sin(np.deg2rad(lat))
    
    ff = 4*np.pi/86164*np.sin(lat*np.pi/180) 
    

    # Compute grid spacing in longitude and latitude
    dlon = np.gradient(lon, axis=1)
    dlat = np.gradient(lat, axis=0)

    # Calculate the derivatives of u and v with respect to longitude and latitude
    du_dlon = np.gradient(u, axis=2) / (dlon * 111.32e3 * np.cos(np.deg2rad(lat)))
    dv_dlat = np.gradient(v, axis=1) / (dlat * 111.32e3)

    # Compute oceanic relative vorticity (zeta)
    zeta = dv_dlat - du_dlon + ff[None, :, :] 

    # Compute potential vorticity (PV)
    pv = zeta / h[None, :, :]

    return pv



def currents_to_relative_vorticity_old(u, v, lon, lat):
    """
    Compute oceanic relative vorticity (zeta) from ocean current components.

    Parameters
    ----------
    u : numpy.ndarray
        3D array representing eastward ocean current velocity (time, x, y).
    v : numpy.ndarray
        3D array representing northward ocean current velocity (time, x, y).
    lon : numpy.ndarray
        1D or 2D array of longitudes corresponding to the x-dimension of u and v.
    lat : numpy.ndarray
        1D or 2D array of latitudes corresponding to the y-dimension of u and v.

    Returns
    -------
    zeta : numpy.ndarray
        3D array of oceanic relative vorticity (time, x, y).
    """
    
    
    if len(np.shape(lon))==1: 
        lon, lat = np.meshgrid(lon,lat)
        
    # Calculate the Coriolis parameter (f) using the latitude 
    ff = 4*np.pi/86164*np.sin(lat*np.pi/180)

    # Compute grid spacing in longitude and latitude
    dlon = np.gradient(lon, axis=1)
    dlat = np.gradient(lat, axis=0)

    # Calculate the derivatives of u and v with respect to longitude and latitude
    du_dlon = np.gradient(u, axis=2) / (dlon * 111.32e3 * np.cos(np.deg2rad(lat)))
    dv_dlat = np.gradient(v, axis=1) / (dlat * 111.32e3)

    # Compute oceanic relative vorticity (zeta)
    zeta = dv_dlat - du_dlon  

    return zeta / ff[None, :, :]


def currents_to_relative_vorticity(u, v, lon, lat):
    """
    Compute oceanic relative vorticity (zeta) from ocean current components.

    Parameters
    ----------
    u : numpy.ndarray
        3D array representing eastward ocean current velocity (time, x, y).
    v : numpy.ndarray
        3D array representing northward ocean current velocity (time, x, y).
    lon : numpy.ndarray
        1D or 2D array of longitudes corresponding to the x-dimension of u and v.
    lat : numpy.ndarray
        1D or 2D array of latitudes corresponding to the y-dimension of u and v.

    Returns
    -------
    rv : numpy.ndarray
        3D array of oceanic relative vorticity (time, x, y).
    """
     

    if len(np.shape(lon))==1: 
        lon, lat = np.meshgrid(lon, lat)

    dx,dy = lonlat2dxdy(lon, lat)

    u = np.moveaxis(u, 0, -1)
    v = np.moveaxis(v, 0, -1)

    _dx = dx[:,:,np.newaxis]
    _dy = dy[:,:,np.newaxis] 


    # Calculate the Coriolis parameter (f) using the latitude 
    ff = 4*np.pi/86164*np.sin(lat*np.pi/180)

    rv = np.zeros_like(u)

    rv = np.gradient(v, axis=1)/_dx - np.gradient(u, axis=0)/_dy

    rv = np.moveaxis(rv, -1, 0)

    f = 4*np.pi/86164*np.sin(lat*np.pi/180)

    rv = rv / f[None,:,:]

    return rv


def lonlat2dxdy(lon,lat):
    """
    Compute spatial grid spacings (dx, dy) from longitude and latitude arrays.

    Parameters
    ----------
    lon : numpy.ndarray
        2D array of longitudes.
    lat : numpy.ndarray
        2D array of latitudes.

    Returns
    -------
    dx : numpy.ndarray
        2D array of grid spacings in the longitudinal direction.
    dy : numpy.ndarray
        2D array of grid spacings in the latitudinal direction.
    """
    
    dlon = np.gradient(lon) 
    dlat = np.gradient(lat)
    dx = np.sqrt((dlon[1]*111000*np.cos(np.deg2rad(lat)))**2
                 + (dlat[1]*111000)**2)
    dy = np.sqrt((dlon[0]*111000*np.cos(np.deg2rad(lat)))**2
                 + (dlat[0]*111000)**2)
    dx[0,:] = dx[1,:]
    dx[-1,: ]= dx[-2,:] 
    dx[:,0] = dx[:,1]
    dx[:,-1] = dx[:,-2]
    dy[0,:] = dy[1,:]
    dy[-1,:] = dy[-2,:] 
    dy[:,0] = dy[:,1]
    dy[:,-1] = dy[:,-2]
    
    return dx,dy
 

def sla_to_ssh(ds, path_mdt='../data/cnes_obs-sl_glo_phy-mdt_my_0.125deg_P20Y_1695393893725_1mdt.nc', ds_mdt=None, name_lonlat=['longitude','latitude']):
    """ 
    Compute sea surface height (SSH) from sea level anomaly (SLA) and mean dynamic topography (MDT).

    Parameters
    ----------
    ds : xarray.Dataset
        3D array representing sea level anomaly (SLA) (time, x, y).
    path_mdt : str, optional
        Path to the MDT dataset file.
    ds_mdt : xarray.Dataset, optional
        2D array representing mean dynamic topography (MDT) (x, y).
    name_lonlat : list of str, optional
        List containing the names of longitude and latitude variables in the dataset.

    Returns
    -------
    ds : xarray.Dataset
        3D array of sea surface height (SSH) (time, x, y).

    """
    if 'ssh' in ds.data_vars:
        print('Warning: ssh variable already exists in provided dataset. Naming new ssh: ssh2')
        
        
    if ds_mdt is None: 
        
        try:  
            ds_mdt = xr.open_dataset(path_mdt)
            
        except:
            print('No mdt provided.')
            raise
            
    lon_min = np.min(ds[name_lonlat[0]])
    lon_max = np.max(ds[name_lonlat[0]])
    lat_min = np.min(ds[name_lonlat[1]])
    lat_max = np.max(ds[name_lonlat[1]])
            
    ds_mdt = regional_zoom(ds_mdt, [lon_min-1,lon_max+1], [lat_min-1,lat_max+1], namelon='longitude', namelat='latitude', change_lon=True)
        

    lon,lat = np.meshgrid(ds[name_lonlat[0]],ds[name_lonlat[1]])

    mdt_interp = interp2d(ds_mdt,{'lon':'longitude','lat':'latitude','var':'mdt'},lon,lat)

    mdt_interp3d = np.repeat(mdt_interp.transpose()[:, :, np.newaxis], ds.time.size, axis=2).transpose()

    if 'ssh' in ds.data_vars:
        ds = ds.assign(ssh2=lambda ds: ds.sla  + mdt_interp3d)
    else:    
        ds = ds.assign(ssh=lambda ds: ds.sla  + mdt_interp3d)

    return ds