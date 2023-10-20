import cartopy.crs as ccrs
import hvplot.xarray
import pandas as pd
import xarray as xr
import numpy as np
import warnings





def regional_zoom(ds_glob, boxlon, boxlat, namelon='lon', namelat='lat', change_lon = True):

    """
    Select a geographical region in a xarray.dataset using longitude and latitude information.

    Parameters
    ----------
    ds_glob : xarray.Dataset
        Input dataset organized in longitudes and latitudes coordinates.
    box_lon : Array of float
        Min and max longitudes of the requested regional zoom.
    box_lat : Array of float
        Min and max latitudes of the requested regional zoom. 
    namelon : str, optional
        Name of the longitude coordinate.
    namelat : str, optional
        Name of the latitude coordinate.
    change_lon : Boolean, optional
        Switches from (0,360) longitude to (-180,180) if True. Default True.

    Returns
    -------
    xarray.Dataset
        The input dataset restricted in the requested geographical region.
    """

    min_lon = boxlon[0] 
    max_lon = boxlon[1] 
    min_lat = boxlat[0]
    max_lat = boxlat[1] 

    import copy
    ds_reg = copy.deepcopy(ds_glob)

    long_attrs = ds_reg[namelon].attrs 
    
    if change_lon: 
        ds_reg[namelon] = ds_reg[namelon].where(ds_reg[namelon]<180,ds_reg[namelon]-360).values
        ds_reg = ds_reg.sortby(ds_reg[namelon])
      
  
    ds_reg = ds_reg.where(ds_reg[namelon]<max_lon,drop=True)
    ds_reg = ds_reg.where(ds_reg[namelon]>min_lon,drop=True)
    ds_reg = ds_reg.where(ds_reg[namelat]<max_lat,drop=True)
    ds_reg = ds_reg.where(ds_reg[namelat]>min_lat,drop=True) 
    
    
    #ds_reg[namelon] = ds_reg[namelon].values%360
    #ds_reg = ds_reg.sortby(ds_reg[namelon])
    ds_reg[namelon].attrs = long_attrs 
    #ds_reg = ds_reg.sortby(ds_reg['time'])

    return ds_reg


def convert_longitude(ds, lon_name):
    """
    Converts longitude coordinates ranging from 0 to 360 to -180 to 180.
    
    Parameters
    ----------
    ds : xarray.Dataset
        The input dataset containing longitude coordinates ranging from 0 to 360.
    lon_name : str
        The name of the longitude coordinate in the dataset.
    
    Returns
    -------
    xarray.Dataset
        The input dataset with longitude coordinates ranging from -180 to 180.
    """
    lon = ds[lon_name]
    if lon.min() < 0 or lon.max() > 360:
        lon = ((lon + 180) % 360) - 180
        ds = ds.assign_coords(**{lon_name: lon})
        message = f"Converted {lon_name} from 0-360 to -180-180"
        warnings.warn(message)
    else:
        message = f"{lon_name} already in -180-180"
        warnings.warn(message, category=UserWarning)
    return ds



def compare_stat_score_map(study_filename, ref_filename,boxlon = [-180, 180],boxlat = [-90, 90], vartype='err_var'):
    """
    Plot statistical score maps.

    Parameters
    ----------
    study_filename : str
        Path to the NetCDF file containing the studied data. 
    ref_filename : str
        Path to the NetCDF file containing the reference data.
    boxlon : Array of floats, optional
        Array containing min and max longitude to plot.
    boxlar : Array of floats, optional
        Array containing min and max latitude to plot.
    vartype : str, optional
        Name of the variable to plot, 'err_var' by default or 'exp_var'.

    Returns
    -------
    holoviews.Layout
        A HoloViews layout containing six plots.
    """

    
    ds_ref_binning_allscale = xr.open_dataset(ref_filename, group='all_scale')
    ds_ref_binning_allscale = regional_zoom(ds_ref_binning_allscale,boxlon, boxlat)
    ds_ref_binning_filtered = xr.open_dataset(ref_filename, group='filtered')
    ds_ref_binning_filtered = regional_zoom(ds_ref_binning_filtered,boxlon, boxlat)
    
    explained_variance_ref_all_scale = 1 - ds_ref_binning_allscale['variance_mapping_err']/ds_ref_binning_allscale['variance_track']
    explained_variance_ref_filtered = 1 - ds_ref_binning_filtered['variance_mapping_err']/ds_ref_binning_filtered['variance_track']
    
    ds_study_binning_allscale = xr.open_dataset(study_filename, group='all_scale')
    ds_study_binning_allscale = regional_zoom(ds_study_binning_allscale,boxlon, boxlat)
    ds_study_binning_filtered = xr.open_dataset(study_filename, group='filtered')
    ds_study_binning_filtered = regional_zoom(ds_study_binning_filtered,boxlon, boxlat)
    
    explained_variance_study_all_scale = 1 - ds_study_binning_allscale['variance_mapping_err']/ds_study_binning_allscale['variance_track']
    explained_variance_study_filtered = 1 - ds_study_binning_filtered['variance_mapping_err']/ds_study_binning_filtered['variance_track']
    
    
    if vartype == 'err_var':
    
    
        fig1 = ds_ref_binning_allscale['variance_mapping_err'].hvplot.quadmesh(x='lon',
                                                                  y='lat',
                                                                  clim=(0, 0.002),
                                                                  cmap='Reds',
                                                                  rasterize=True,
                                                                  title='Error variance [All scale]: '+ds_ref_binning_allscale.attrs['method'])

        fig2 = ds_ref_binning_filtered['variance_mapping_err'].hvplot.quadmesh(x='lon',
                                                                  y='lat',
                                                                  clim=(0, 0.002),
                                                                  cmap='Reds',
                                                                  rasterize=True,
                                                                  title='Error variance [65:500km]: '+ds_ref_binning_filtered.attrs['method'])
        
    
        fig3 = ds_study_binning_allscale['variance_mapping_err'].hvplot.quadmesh(x='lon',
                                                                  y='lat',
                                                                  clim=(0, 0.002),
                                                                  cmap='Reds',
                                                                  rasterize=True,
                                                                  title='Error variance [All scale]: '+ds_study_binning_allscale.attrs['method'])

        fig4 = ds_study_binning_filtered['variance_mapping_err'].hvplot.quadmesh(x='lon',
                                                                  y='lat',
                                                                  clim=(0, 0.002),
                                                                  cmap='Reds',
                                                                  rasterize=True,
                                                                  title='Error variance [65:500km]: '+ds_study_binning_filtered.attrs['method'])

        fig5 = (100*(ds_study_binning_allscale['variance_mapping_err'] - ds_ref_binning_allscale['variance_mapping_err'])/ds_ref_binning_allscale['variance_mapping_err']).hvplot.quadmesh(x='lon',
                                                                  y='lat',
                                                                  clim=(-20, 20),
                                                                  cmap='coolwarm',
                                                                  rasterize=True,
                                                                  title='Error variance [All scale]')

        fig6 = (100*(ds_study_binning_filtered['variance_mapping_err'] - ds_ref_binning_filtered['variance_mapping_err'])/ds_ref_binning_filtered['variance_mapping_err']).hvplot.quadmesh(x='lon',
                                                                  y='lat',
                                                                  clim=(-20, 20),
                                                                  cmap='coolwarm',
                                                                  rasterize=True,
                                                                  title='Error variance [65:500km]')
    
    if vartype == 'expl_var':
    
        fig1 = explained_variance_ref_all_scale.hvplot.quadmesh(x='lon',
                                                                  y='lat',
                                                                  clim=(0, 1),
                                                                  cmap='RdYlGn',
                                                                  rasterize=True,
                                                                  title='Explained variance [All scale]: '+ds_ref_binning_allscale.attrs['method']
                                                               )

        fig2 = explained_variance_ref_filtered.hvplot.quadmesh(x='lon',
                                                                  y='lat',
                                                                  clim=(0, 1),
                                                                  cmap='RdYlGn',
                                                                  rasterize=True,
                                                                  title='Explained variance [65:500km]: '+ds_ref_binning_allscale.attrs['method']
                                                               )
        
    
        fig3 = explained_variance_study_all_scale.hvplot.quadmesh(x='lon',
                                                                  y='lat',
                                                                  clim=(0, 1),
                                                                  cmap='RdYlGn',
                                                                  rasterize=True,
                                                                  title='Explained variance [All scale]: '+ds_study_binning_allscale.attrs['method']
                                                               )

        fig4 = explained_variance_study_filtered.hvplot.quadmesh(x='lon',
                                                                  y='lat',
                                                                  clim=(0, 1),
                                                                  cmap='RdYlGn',
                                                                  rasterize=True,
                                                                  title='Explained variance [65:500km]: '+ds_study_binning_allscale.attrs['method']
                                                               )

        fig5 = (explained_variance_study_all_scale - explained_variance_ref_all_scale).hvplot.quadmesh(x='lon',
                                                                  y='lat',
                                                                  clim=(-0.2, 0.2),
                                                                  cmap='coolwarm_r',
                                                                  rasterize=True,
                                                                  title='Gain(+)/Loss(-) Explained variance [All scale]')

        fig6 = (explained_variance_study_filtered - explained_variance_ref_filtered).hvplot.quadmesh(x='lon',
                                                                  y='lat',
                                                                  clim=(-0.2, 0.2),
                                                                  cmap='coolwarm_r',
                                                                  rasterize=True,
                                                                  title='Gain(+)/Loss(-) Explained variance [65:500km]')
    
#     fig5 = ds_binning_allscale['rmse'].hvplot.quadmesh(x='lon',
#                                                               y='lat',
#                                                               clim=(0, 0.1),
#                                                               cmap='Reds',
#                                                               rasterize=True,
#                                                               title='RMSE [All scale]')
    
#     fig6 = ds_binning_filtered['rmse'].hvplot.quadmesh(x='lon',
#                                                               y='lat',
#                                                               clim=(0, 0.1),
#                                                               cmap='Reds',
#                                                               rasterize=True,
#                                                               title='RMSE [65:500km]')
    
    return (fig1 + fig2 + fig3 + fig4 + fig5 + fig6).cols(2)



def compare_psd_score(study_filename, ref_filename,boxlon = [-180, 180],boxlat = [-90, 90]):
    """
    Compare Power Spectral Density (PSD) scores between a study dataset and a reference dataset.
 
    Parameters
    ----------
    study_filename : str
        Filepath to the study dataset NetCDF file for comparison.
    ref_filename : str
        Filepath to the reference dataset NetCDF file for comparison.
    boxlon : list of float, optional
        Longitude bounding box for regional zoom, by default [-180, 180].
    boxlat : list of float, optional
        Latitude bounding box for regional zoom, by default [-90, 90].

    Returns
    -------
    holoviews.Layout
        A HoloViews layout containing three plots:
        - Effective resolution of the reference dataset.
        - Effective resolution of the study dataset.
        - Gain(-)/loss(+) of effective resolution compared to the reference dataset.

    Examples
    --------
    >>> study_file = "study_data.nc"
    >>> ref_file = "reference_data.nc"
    >>> compare_psd_score(study_file, ref_file, boxlon=[-120, -70], boxlat=[20, 50])
    """
    
    ds_ref = xr.open_dataset(ref_filename)
    ds_ref = regional_zoom(ds_ref,boxlon, boxlat)
    ds_study = xr.open_dataset(study_filename)
    ds_study = regional_zoom(ds_study,boxlon, boxlat)
    
    fig0 = ds_ref.effective_resolution.hvplot.quadmesh(x='lon', 
                                                   y='lat', 
                                                   cmap='Spectral_r', 
                                                   clim=(100, 500), 
                                                   title='Effective resolution [km]: '+ds_ref.attrs['method'],
                                                   rasterize=True, 
                                                   projection=ccrs.PlateCarree(), 
                                                   project=True, 
                                                   geo=True, 
                                                   coastline=True)
    
    
    fig1 = ds_study.effective_resolution.hvplot.quadmesh(x='lon', 
                                                   y='lat', 
                                                   cmap='Spectral_r', 
                                                   clim=(100, 500), 
                                                   title='Effective resolution [km]: '+ds_study.attrs['method'],
                                                   rasterize=True, 
                                                   projection=ccrs.PlateCarree(), 
                                                   project=True, 
                                                   geo=True, 
                                                   coastline=True)
    
    fig2 = (100*(ds_study.effective_resolution - ds_ref.effective_resolution)/ds_ref.effective_resolution).hvplot.quadmesh(x='lon', 
                                                   y='lat', 
                                                   cmap='coolwarm', 
                                                   clim=(-20, 20), 
                                                   title='Gain(-)/loss(+) Effective resolution [%]',
                                                   rasterize=True, 
                                                   projection=ccrs.PlateCarree(), 
                                                   project=True, 
                                                   geo=True, 
                                                   coastline=True)
    
    return (fig0 + fig1 + fig2).cols(1)


