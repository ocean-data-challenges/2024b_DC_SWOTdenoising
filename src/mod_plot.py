import cartopy.crs as ccrs
import hvplot.xarray
import pandas as pd
import xarray as xr
import numpy as np
import cartopy
import matplotlib.axes
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import cartopy.feature as cfeature 
from src.mod_compare import *
 

def plot_stat_score_map(filename):
    """
    Generate and display error and explained variance maps.

    Parameters
    ----------
    filename : str
        Path to the input NetCDF file containing required data.

    Returns
    -------
    holoviews layout
        Layout containing quadmesh plots for error and explained variance maps.
    """
    
    ds_binning_allscale = xr.open_dataset(filename, group='all_scale')
    ds_binning_filtered = xr.open_dataset(filename, group='filtered')
    
    fig1 = ds_binning_allscale['variance_mapping_err'].hvplot.quadmesh(x='lon',
                                                              y='lat',
                                                              clim=(0, 0.002),
                                                              cmap='Reds',
                                                              rasterize=True,
                                                              title='Error variance [All scale]')
    
    fig2 = ds_binning_filtered['variance_mapping_err'].hvplot.quadmesh(x='lon',
                                                              y='lat',
                                                              clim=(0, 0.002),
                                                              cmap='Reds',
                                                              rasterize=True,
                                                              title='Error variance [65:500km]')
    
    fig3 = (1 - ds_binning_allscale['variance_mapping_err']/ds_binning_allscale['variance_track']).hvplot.quadmesh(x='lon',
                                                              y='lat',
                                                              clim=(0, 1),
                                                              cmap='RdYlGn',
                                                              rasterize=True,
                                                              title='Explained variance [All scale]')
    
    fig4 = (1 - ds_binning_filtered['variance_mapping_err']/ds_binning_filtered['variance_track']).hvplot.quadmesh(x='lon',
                                                              y='lat',
                                                              clim=(0, 1),
                                                              cmap='RdYlGn',
                                                              rasterize=True,
                                                              title='Explained variance [65:500km]')
    
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
    
    return (fig1 + fig2 + fig3 + fig4).cols(2)


def plot_stat_score_timeseries(filename):
    """
    Generate and display timeseries of error and explained variance.

    Parameters
    ----------
    filename : str
        Path to the input NetCDF file containing required data.

    Returns
    -------
    holoviews layout
        Layout containing line plots for error and explained variance timeseries.
    """
    
    ds_binning_allscale = xr.open_dataset(filename, group='all_scale')
    ds_binning_filtered = xr.open_dataset(filename, group='filtered')
    
    fig0 = ds_binning_allscale['timeserie_variance_mapping_err'].hvplot.line(x='time', 
                                                                         y='timeserie_variance_mapping_err',  
                                                                         label='All scale',
                                                                             grid=True,
                                                                         title='Daily averaged Error variance')*ds_binning_filtered['timeserie_variance_mapping_err'].hvplot.line(x='time', 
                                                                                                                                                                y='timeserie_variance_mapping_err', 
                                                                         label='Filtered',
                                                                         grid=True,
                                                                         title='Daily averaged Error variance')
    
    
    ds_binning_allscale['explained_variance_score'] =  1. - (ds_binning_allscale['timeserie_variance_mapping_err']/ds_binning_allscale['timeserie_variance_track'])
    ds_binning_filtered['explained_variance_score'] =  1. - (ds_binning_filtered['timeserie_variance_mapping_err']/ds_binning_filtered['timeserie_variance_track'])
    
    fig1 = ds_binning_allscale['explained_variance_score'].hvplot.line(x='time', 
                                                                       y='explained_variance_score',
                                                                       label='All scale',
                                                                       grid=True,
                                                                       title='Explained variance score', 
                                                                       ylim=(0, 1))*ds_binning_filtered['explained_variance_score'].hvplot.line(x='time', 
                                                                                                                                                y='explained_variance_score',
                                                                                                                                                label='Filtered',
                                                                                                                                                grid=True,
                                                                                                                                                title='Explained variance score', 
                                                                                                                                                ylim=(0, 1))
    
    return (fig0 + fig1).cols(1)


def plot_stat_by_regimes(stat_output_filename):
    """
    Generate and display statistical summary by geographical regions.

    Parameters
    ----------
    stat_output_filename : str
        Path to the input NetCDF file containing statistical data.

    Returns
    -------
    pandas DataFrame
        DataFrame containing the statistical summary by regions and variables.
    """
    
    my_dictionary = {}
    for region in ['coastal', 'offshore_highvar', 'offshore_lowvar', 'equatorial_band', 'arctic', 'antarctic']:
        my_dictionary[f'{region}'] = {}
        for var_name in ['mapping_err', 'sla_unfiltered', 'mapping_err_filtered', 'sla_filtered']:
        
            ds = xr.open_dataset(stat_output_filename, group=f'{region}_{var_name}')

            my_dictionary[f'{region}'][f'{var_name}_var [m²]'] =  ds['variance'].values[0]
            #my_dictionary[f'{region}'][f'{var_name}_rms'] =  ds['rmse'].values[0]
    
    for region in ['coastal', 'offshore_highvar', 'offshore_lowvar', 'equatorial_band', 'arctic', 'antarctic']:
        my_dictionary[region]['var_score_allscale'] = 1. - my_dictionary[region]['mapping_err_var [m²]']/my_dictionary[region]['sla_unfiltered_var [m²]']
        my_dictionary[region]['var_score_filtered'] = 1. - my_dictionary[region]['mapping_err_filtered_var [m²]']/my_dictionary[region]['sla_filtered_var [m²]']
    
    return pd.DataFrame(my_dictionary.values(), index=my_dictionary.keys())



def plot_stat_uv_by_regimes(stat_output_filename):
    """
    Generate and display statistical summary of velocity components by geographical regions.

    Parameters
    ----------
    stat_output_filename : str
        Path to the input NetCDF file containing statistical data.

    Returns
    -------
    pandas DataFrame
        DataFrame containing the statistical summary by regions and velocity components.
    """
    my_dictionary = {}
    for region in ['coastal', 'offshore_highvar', 'offshore_lowvar', 'equatorial_band', 'arctic', 'antarctic']:
        my_dictionary[f'{region}'] = {}
        for var_name in ['mapping_err_u', 'mapping_err_v', 'ugos_interpolated', 'EWCT', 'vgos_interpolated', 'NSCT']:
        
            ds = xr.open_dataset(stat_output_filename, group=f'{region}_{var_name}')

            my_dictionary[f'{region}'][f'{var_name}_var [m²/s²]'] =  ds['variance'].values[0]
            #my_dictionary[f'{region}'][f'{var_name}_rms'] =  ds['rmse'].values[0]
    
    for region in ['coastal', 'offshore_highvar', 'offshore_lowvar', 'equatorial_band', 'arctic', 'antarctic']:
        my_dictionary[region]['var_score_u_allscale'] = 1. - my_dictionary[region]['mapping_err_u_var [m²/s²]']/my_dictionary[region]['EWCT_var [m²/s²]']
        my_dictionary[region]['var_score_v_allscale'] = 1. - my_dictionary[region]['mapping_err_v_var [m²/s²]']/my_dictionary[region]['NSCT_var [m²/s²]']
    
    return pd.DataFrame(my_dictionary.values(), index=my_dictionary.keys())
    


def plot_effective_resolution(filename):
    """
    Generate and display a quadmesh plot of effective resolution.

    Parameters
    ----------
    filename : str
        Path to the input NetCDF file containing required data.

    Returns
    -------
    holoviews quadmesh
        Quadmesh plot of effective resolution.
    """
    
    ds = xr.open_dataset(filename)
    
    fig0 = ds.effective_resolution.hvplot.quadmesh(x='lon', 
                                                   y='lat', 
                                                   cmap='Spectral_r', 
                                                   clim=(100, 500), 
                                                   title='Effective resolution [km]',
                                                   rasterize=True, 
                                                   projection=ccrs.PlateCarree(), 
                                                   project=True, 
                                                   geo=True, 
                                                   coastline=True)
    
    return fig0

 


def plot_effective_resolution_png(filename,region='glob',box_lonlat=None, change_lon=True):
    """
    Generate and save a PNG image of effective resolution.

    Parameters
    ----------
    filename : str
        Path to the input NetCDF file containing required data.
    region : str, optional
        Region name for the image, by default 'glob'.
    box_lonlat : dict, optional
        Dictionary containing 'lon_min', 'lon_max', 'lat_min', and 'lat_max' values,
        defining the bounding box for the image.

    Returns
    -------
    None
    """

    ds = xr.open_dataset(filename)
    
    
    try : 
        method_name = ds.attrs['method']
    except:
        method_name = ''
    
    if box_lonlat is not None:
        lon_min = box_lonlat['lon_min']
        lon_max = box_lonlat['lon_max']
        lat_min = box_lonlat['lat_min']
        lat_max = box_lonlat['lat_max']
        #ds = ds.sel({'lon':slice(lon_min,lon_max),'lat':slice(lat_min,lat_max)}) 
        
        ds = regional_zoom(ds, [lon_min,lon_max], [lat_min,lat_max], namelon='lon', namelat='lat', change_lon=change_lon)
    
    
    fig, axs = plt.subplots(nrows=1,ncols=1,
                        subplot_kw={'projection': ccrs.PlateCarree()},
                        figsize=(5.5,3.25))

    
    vmin = 100.
    vmax= 500.
    p0 = axs.pcolormesh(ds.lon, ds.lat, ds.effective_resolution, vmin=vmin, vmax=vmax, cmap='Spectral_r')
    axs.set_title('SSH Map Effective resolution')
    axs.add_feature(cfeature.LAND, color='w', zorder=12)
    axs.coastlines(resolution='10m', lw=0.5, zorder=13)
    # optional add grid lines
    p0.axes.gridlines(color='black', alpha=0., linestyle='--')
    # draw parallels/meridiens and write labels
    gl = p0.axes.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                          linewidth=0.1, color='black', alpha=0.5, linestyle='--')
    # adjust labels to taste
    gl.xlabels_top = False
    gl.ylabels_right = False
    gl.ylocator = mticker.FixedLocator([-90, -60, -30, 0, 30, 60, 90])
    gl.xformatter = LONGITUDE_FORMATTER
    gl.yformatter = LATITUDE_FORMATTER
    gl.xlabel_style = {'size': 10, 'color': 'black'}
    gl.ylabel_style = {'size': 10, 'color': 'black'}
    
    
    cax = fig.add_axes([0.92, 0.25, 0.04, 0.6])
    fig.colorbar(p0, cax=cax, orientation='vertical')
    cax.set_ylabel('Effective resolution [km]', fontweight='bold')
    
    fig.subplots_adjust(bottom=0.2, top=0.9, left=0.1, right=0.9,
                    wspace=0.02, hspace=0.01)
    
    plt.savefig("../figures/Maps_"+str(method_name)+"_effres_"+region+".png", bbox_inches='tight')
    
    print('Averaged effective resolution:',np.nanmean(ds.effective_resolution),'km')



def plot_psd_scores(filename):
    """
    Generate and display plots related to Power Spectral Density (PSD) scores.

    Parameters
    ----------
    filename : str
        Path to the input NetCDF file containing required data.

    Returns
    -------
    holoviews layout
        Layout containing line plots for PSD scores and related metrics.
    """
 


    ds = xr.open_dataset(filename)
    
    
    try : 
        method_name = ds.attrs['method']
    except:
        method_name = ''
    
    if box_lonlat is not None:
        lon_min = box_lonlat['lon_min']
        lon_max = box_lonlat['lon_max']
        lat_min = box_lonlat['lat_min']
        lat_max = box_lonlat['lat_max']
        #ds = ds.sel({'lon':slice(lon_min,lon_max),'lat':slice(lat_min,lat_max)}) 
        
        ds = regional_zoom(ds, [lon_min,lon_max], [lat_min,lat_max], namelon='lon', namelat='lat', change_lon=True)
    
    
    fig, axs = plt.subplots(nrows=1,ncols=1,
                        subplot_kw={'projection': ccrs.PlateCarree()},
                        figsize=(5.5,3.25))

    
    vmin = 100.
    vmax= 500.
    p0 = axs.pcolormesh(ds.lon, ds.lat, ds.effective_resolution, vmin=vmin, vmax=vmax, cmap='Spectral_r')
    axs.set_title('SSH Map Effective resolution')
    axs.add_feature(cfeature.LAND, color='w', zorder=12)
    axs.coastlines(resolution='10m', lw=0.5, zorder=13)
    # optional add grid lines
    p0.axes.gridlines(color='black', alpha=0., linestyle='--')
    # draw parallels/meridiens and write labels
    gl = p0.axes.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                          linewidth=0.1, color='black', alpha=0.5, linestyle='--')
    # adjust labels to taste
    gl.xlabels_top = False
    gl.ylabels_right = False
    gl.ylocator = mticker.FixedLocator([-90, -60, -30, 0, 30, 60, 90])
    gl.xformatter = LONGITUDE_FORMATTER
    gl.yformatter = LATITUDE_FORMATTER
    gl.xlabel_style = {'size': 10, 'color': 'black'}
    gl.ylabel_style = {'size': 10, 'color': 'black'}
    
    
    cax = fig.add_axes([0.92, 0.25, 0.04, 0.6])
    fig.colorbar(p0, cax=cax, orientation='vertical')
    cax.set_ylabel('Effective resolution [km]', fontweight='bold')
    
    fig.subplots_adjust(bottom=0.2, top=0.9, left=0.1, right=0.9,
                    wspace=0.02, hspace=0.01)
    
    plt.savefig("../figures/Maps_"+str(method_name)+"_effres_"+region+".png", bbox_inches='tight')


def plot_psd_scores(filename):
    """
    Generate and display plots related to Power Spectral Density (PSD) scores.

    Parameters
    ----------
    filename : str
        Path to the input NetCDF file containing required data.

    Returns
    -------
    holoviews layout
        Layout containing line plots for PSD scores and related metrics.
    """
    
    ds = xr.open_dataset(filename)
    
    ds = ds.assign_coords({'wavelenght':1./ds['wavenumber']})
    ds['psd_score'] = 1 - ds['psd_diff']/ds['psd_ref']
    ds['psd_ratio'] = ds['psd_study']/ds['psd_ref']
    
    # Compute noise level from PSD as in Dufau et al. (2018), for example
    def func_y0(x, noise):
        return 0*x + noise
    
    noise = (ds.where(1./ds.wavenumber < 23, drop=True)).psd_ref.curvefit(coords='wavenumber', func=func_y0).curvefit_coefficients[:, :, 0]
    ds['noise'] = noise.expand_dims(dim={'wavenumber':ds.wavenumber.size})
    
    fig1 = ((ds.psd_ref.hvplot.line(x='wavelenght', 
                                y='psd_ref', 
                                label='PSD_alongtrack', 
                                xlabel='wavelenght [km]', 
                                ylabel='PSD', 
                                logx=True, 
                                logy=True, 
                                flip_xaxis=True,
                                title='Power spectral density', 
                                xticks=[20, 50, 100, 200, 300, 400, 600, 800], 
                                grid=True))*(ds.noise.hvplot.line(x='wavelenght', 
                                                                  y='noise', 
                                                                  label='NOISE_alongtrack', 
                                                                  xlabel='wavelenght [km]', 
                                                                  ylabel='NOISE', 
                                                                  logx=True, 
                                                                  logy=True, 
                                                                  flip_xaxis=True,
                                                                  xticks=[20, 50, 100, 200, 300, 400, 600, 800], 
                                                                  grid=True))*(ds.psd_study.hvplot.line(x='wavelenght', 
                                                                                                        y='psd_study', 
                                                                                                        label='PSD_map', 
                                                                                                        logx=True, 
                                                                                                        logy=True, 
                                                                                                        flip_xaxis=True))*(ds.psd_diff.hvplot.line(x='wavelenght', 
                                                                                                                                                   y='psd_diff', 
                                                                                                                                                   label='PSD_err', 
                                                                                                                                                   logx=True, 
                                                                                                                                                   logy=True, 
                                                                                                                                                   flip_xaxis=True))).opts(width=500)
    
    
    fig2 = ((ds.psd_ratio.hvplot.line(x='wavelenght', 
                                      y='psd_ratio', 
                                      xlabel='wavelenght [km]', 
                                      ylabel='PSD_ratio',
                                      ylim=(0,1),
                                      label='PSD_map/PSD_ref', 
                                      logx=True, 
                                      logy=False, 
                                      flip_xaxis=True, 
                                      title='PSD ratio', 
                                      xticks=[20, 50, 100, 200, 300, 400, 600, 800], 
                                      grid=True))*((0.5*ds.coherence/ds.coherence).hvplot.line(x='wavelenght', 
                                                                                               y='coherence', 
                                                                                               c='r', 
                                                                                               line_width=0.5, 
                                                                                               logx=True, 
                                                                                               logy=False, 
                                                                                               flip_xaxis=True))).opts(width=500)
    
    
    fig3 = (ds.psd_score.hvplot.line(x='wavelenght', 
                                     y='psd_score', 
                                     xlabel='wavelenght [km]', 
                                     ylabel='PSD_score', 
                                     logx=True, 
                                     logy=False, 
                                     flip_xaxis=True, 
                                     title='PSD_score = 1. - PSD_err/PSD_ref', 
                                     xticks=[20, 50, 100, 200, 300, 400, 600, 800], 
                                     grid=True)*((0.5*ds.coherence/ds.coherence).hvplot.line(x='wavelenght', 
                                                                                             y='coherence', 
                                                                                             c='r', 
                                                                                             line_width=0.5, 
                                                                                             logx=True, 
                                                                                             logy=False, 
                                                                                             flip_xaxis=True))).opts(width=500)
    
    fig4 = (ds.coherence.hvplot.line(x='wavelenght', 
                                     y='coherence', 
                                     xlabel='wavelenght [km]', 
                                     ylabel='MSC', 
                                     logx=True, 
                                     logy=False, 
                                     flip_xaxis=True, 
                                     title='Magnitude Squared Coherence', 
                                     xticks=[20, 50, 100, 200, 300, 400, 600, 800], 
                                     grid=True)*((0.5*ds.coherence/ds.coherence).hvplot.line(x='wavelenght', 
                                                                                             y='coherence', 
                                                                                             c='r', 
                                                                                             line_width=0.5, 
                                                                                             logx=True, 
                                                                                             logy=False, 
                                                                                             flip_xaxis=True))).opts(width=500)
    
    
    
    
    
    return (fig1 + fig2 + fig3 + fig4).cols(2)




def plot_stat_score_map_uv(filename):
    """
    Generate and display quadmesh plots of variance and explained variance for zonal and meridional currents.

    Parameters
    ----------
    filename : str
        Path to the input NetCDF file containing required data.

    Returns
    -------
    holoviews layout
        Layout containing quadmesh plots for variance and explained variance of zonal and meridional currents.
    """
    
    ds_binning_allscale = xr.open_dataset(filename, group='all_scale')
    
    fig1 = ds_binning_allscale['variance_mapping_err_u'].hvplot.quadmesh(x='lon',
                                                              y='lat',
                                                              clim=(0, 0.1),
                                                              cmap='Reds',
                                                              rasterize=True,
                                                              title='Error variance zonal current [All scale]')
    
    fig2 = ds_binning_allscale['variance_mapping_err_v'].hvplot.quadmesh(x='lon',
                                                              y='lat',
                                                              clim=(0, 0.1),
                                                              cmap='Reds',
                                                              rasterize=True,
                                                              title='Error variance meridional current [All scale]')
    
    fig3 = (1 - ds_binning_allscale['variance_mapping_err_u']/ds_binning_allscale['variance_drifter_u']).hvplot.quadmesh(x='lon',
                                                              y='lat',
                                                              clim=(0, 1),
                                                              cmap='RdYlGn',
                                                              rasterize=True,
                                                              title='Explained variance zonal current [All scale]')
    
    fig4 = (1 - ds_binning_allscale['variance_mapping_err_v']/ds_binning_allscale['variance_drifter_v']).hvplot.quadmesh(x='lon',
                                                              y='lat',
                                                              clim=(0, 1),
                                                              cmap='RdYlGn',
                                                              rasterize=True,
                                                              title='Explained variance meridional current [All scale]')
    
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
    
    return (fig1 + fig2 + fig3 + fig4).cols(2)



def plot_psd_scores_currents(filename):
    """
    Generate and display plots related to Power Spectral Density (PSD) scores for current data.

    Parameters
    ----------
    filename : str
        Path to the input NetCDF file containing required data.

    Returns
    -------
    holoviews layout
        Layout containing quadmesh plots for various PSD-related metrics for current data.
    """
    
    ds_psd = xr.open_dataset(filename)
    
    
    fig1 = np.log10(ds_psd.psd_ref).hvplot.quadmesh(x='wavenumber', y='lat', clim=(-4, 0), cmap='Spectral_r', width=400, height=600, title='Rotary Spectra Drifters', ylim=(-60, 60))
    fig2 = np.log10(ds_psd.psd_study).hvplot.quadmesh(x='wavenumber', y='lat', clim=(-4, 0), cmap='Spectral_r', width=400, height=600, title='Rotary Spectra Maps', ylim=(-60, 60))
    fig3 = np.log10(ds_psd.psd_diff).hvplot.quadmesh(x='wavenumber', y='lat', cmap='Reds', clim=(-4, 0), width=400, height=600, title='Rotary Spectra Error', ylim=(-60, 60))
    fig4 = ds_psd.coherence.hvplot.quadmesh(x='wavenumber', y='lat', cmap='RdYlGn', clim=(0, 1), width=400, height=600, title='Coherence', ylim=(-60, 60))
    fig5 = (1. - ds_psd.psd_diff/ds_psd.psd_ref).hvplot.quadmesh(x='wavenumber', y='lat', cmap='RdYlGn', clim=(0, 1), width=400, height=600, title='PSDerr/PSDref', ylim=(-60, 60))
    
    return (fig1+fig2+fig3+fig4+fig5).cols(3) 


def plot_psd_scores_currents_1D(filename):
    """
    Generate and display 1D line plots of Power Spectral Density (PSD) scores for current data.

    Parameters
    ----------
    filename : str
        Path to the input NetCDF file containing required data.

    Returns
    -------
    holoviews layout
        Layout containing line plots for various PSD-related metrics for current data.
    """
    
    ds_psd = xr.open_dataset(filename)
    
    ds_psd['psd_err_psd_ref'] = (1. - ds_psd.psd_diff/ds_psd.psd_ref)
    fig1 = ds_psd.hvplot.line(x='wavenumber', y='psd_ref', logy=True, grid=True, label='PSD drifters', width=600)
    fig2 = ds_psd.hvplot.line(x='wavenumber', y='psd_study', logy=True, grid=True, label='PSD maps', width=600)
    
    fig3 = ds_psd.hvplot.line(x='wavenumber', y='coherence', grid=True, label='Coherence', width=600, ylim=(0,1))
    fig4 = ds_psd.hvplot.line(x='wavenumber', y='psd_err_psd_ref', grid=True, label='PSD_err/PSDref', width=600, ylim=(0,1))
    return (fig1*fig2 + fig3*fig4).cols(2)


def plot_polarization(filename):
    """
    Generate and display polarization plots of rotary spectrum for current data.

    Parameters
    ----------
    filename : str
        Path to the input NetCDF file containing required data.

    Returns
    -------
    holoviews layout
        Layout containing quadmesh plots for polarization of the rotary spectrum.
    """
    
    ds_psd = xr.open_dataset(filename)
    
    Splus_ref = ds_psd.psd_ref.where(ds_psd.wavenumber > 0, drop=True)
    Sminus_ref = ds_psd.psd_ref.where(ds_psd.wavenumber < 0, drop=True)
    Sminus_ref = np.flip(Sminus_ref, axis=1)
    Sminus_ref['wavenumber'] = np.abs(Sminus_ref['wavenumber'])
    Sminus_ref = Sminus_ref.where(Sminus_ref.wavenumber == Splus_ref.wavenumber, drop=True)
    
    r_ref = (Sminus_ref - Splus_ref)/(Sminus_ref + Splus_ref)
    
    Splus_study = ds_psd.psd_study.where(ds_psd.wavenumber > 0, drop=True)
    Sminus_study = ds_psd.psd_study.where(ds_psd.wavenumber < 0, drop=True)
    Sminus_study = np.flip(Sminus_study, axis=1)
    Sminus_study['wavenumber'] = np.abs(Sminus_study['wavenumber'])
    Sminus_study = Sminus_study.where(Sminus_study.wavenumber == Splus_study.wavenumber, drop=True)
    
    r_study = (Sminus_study - Splus_study)/(Sminus_study + Splus_study)
    
    fig1 = r_ref.hvplot.quadmesh(x='wavenumber', y='lat', clim=(-1, 1), cmap='Spectral_r', width=400, height=600, title='Polarization of the rotary spectrum')
    fig2 = r_study.hvplot.quadmesh(x='wavenumber', y='lat', clim=(-1, 1), cmap='Spectral_r', width=400, height=600, title='Polarization of the rotary spectrum')
    
    return fig1 + fig2





def plot_stat_score_map_uv_png(filename,region='glob',box_lonlat=None):    
    """
    Plot zonal and meridional current variance and explained variance maps as PNG files.

    Parameters
    ----------
    filename : str
        The name of the input NetCDF file.
    region : str, optional
        The region to plot ('glob' for global, or specify a region name).
    box_lonlat : dict, optional
        A dictionary with 'lon_min', 'lon_max', 'lat_min', and 'lat_max' keys defining a bounding box for the plot.

    Returns
    -------
    None
    """
    
    ds_binning_allscale = xr.open_dataset(filename, group='all_scale')
    
    
    if box_lonlat is not None:
        lon_min = box_lonlat['lon_min']
        lon_max = box_lonlat['lon_max']
        lat_min = box_lonlat['lat_min']
        lat_max = box_lonlat['lat_max'] 
        ds_binning_allscale = regional_zoom(ds_binning_allscale, [lon_min,lon_max], [lat_min,lat_max], namelon='lon', namelat='lat', change_lon=True)
    
    
    
    try : 
        method_name = ds_binning_allscale.attrs['method']
    except:
        method_name = ''
        
    
    fig, axs = plt.subplots(nrows=1,ncols=2,
                        subplot_kw={'projection': ccrs.PlateCarree()},
                        figsize=(11,7.5))

    axs=axs.flatten()
    
    vmin = 0.
    vmax= 0.1
    p0 = axs[0].pcolormesh(ds_binning_allscale.lon, ds_binning_allscale.lat, ds_binning_allscale.variance_mapping_err_u, vmin=vmin, vmax=vmax, cmap='Reds')
    axs[0].set_title('Zonal current [All scale]')
    axs[0].coastlines(resolution='10m', lw=0.5)
    # optional add grid lines
    p0.axes.gridlines(color='black', alpha=0., linestyle='--')
    # draw parallels/meridiens and write labels
    gl = p0.axes.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                          linewidth=0.1, color='black', alpha=0.5, linestyle='--')
    # adjust labels to taste
    gl.xlabels_top = False
    gl.ylabels_right = False
    gl.xlabels_bottom = False
    gl.ylocator = mticker.FixedLocator([-90, -60, -30, 0, 30, 60, 90])
    gl.xformatter = LONGITUDE_FORMATTER
    gl.yformatter = LATITUDE_FORMATTER
    gl.xlabel_style = {'size': 10, 'color': 'black'}
    gl.ylabel_style = {'size': 10, 'color': 'black'}
    
    p1 = axs[1].pcolormesh(ds_binning_allscale.lon, ds_binning_allscale.lat, ds_binning_allscale.variance_mapping_err_v, vmin=vmin, vmax=vmax, cmap='Reds')
    axs[1].set_title('Meridional current [All scale]')
    axs[1].coastlines(resolution='10m', lw=0.5)
    # optional add grid lines
    p1.axes.gridlines(color='black', alpha=0., linestyle='--')
    # draw parallels/meridiens and write labels
    gl = p1.axes.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                          linewidth=0.1, color='black', alpha=0.5, linestyle='--')
    # adjust labels to taste
    gl.xlabels_top = False
    gl.ylabels_right = False
    gl.ylabels_left = False
    gl.xlabels_bottom = False
    gl.ylocator = mticker.FixedLocator([-90, -60, -30, 0, 30, 60, 90])
    gl.xformatter = LONGITUDE_FORMATTER
    gl.yformatter = LATITUDE_FORMATTER
    gl.xlabel_style = {'size': 10, 'color': 'black'}
    gl.ylabel_style = {'size': 10, 'color': 'black'}
    
    cax = fig.add_axes([0.92, 0.42, 0.02, 0.25])
    cbar = fig.colorbar(p1, cax=cax, orientation='vertical')
    cax.set_ylabel('Error variance [m$^2$.s$^{-2}$]', fontweight='bold')
    
    plt.savefig("../figures/Maps_"+str(method_name)+"_errvar_"+region+"_uv.png", bbox_inches='tight')
    
    fig.subplots_adjust(bottom=0.2, top=0.9, left=0.1, right=0.9,
                    wspace=0.02, hspace=0.01) 
    
    fig, axs = plt.subplots(nrows=1,ncols=2,
                        subplot_kw={'projection': ccrs.PlateCarree()},
                        figsize=(11,7.5))

    axs=axs.flatten()
    
    vmin = 0.
    vmax= 1
    p2 = axs[0].pcolormesh(ds_binning_allscale.lon, ds_binning_allscale.lat, (1 - ds_binning_allscale['variance_mapping_err_u']/ds_binning_allscale['variance_drifter_u']), vmin=vmin, vmax=vmax, cmap='RdYlGn')
    axs[0].set_title('Zonal current [All scale]')
    axs[0].coastlines(resolution='10m', lw=0.5)
    # optional add grid lines
    p2.axes.gridlines(color='black', alpha=0., linestyle='--')
    # draw parallels/meridiens and write labels
    gl = p2.axes.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                          linewidth=0.1, color='black', alpha=0.5, linestyle='--')
    # adjust labels to taste
    gl.xlabels_top = False
    gl.ylabels_right = False
    gl.ylocator = mticker.FixedLocator([-90, -60, -30, 0, 30, 60, 90])
    gl.xformatter = LONGITUDE_FORMATTER
    gl.yformatter = LATITUDE_FORMATTER
    gl.xlabel_style = {'size': 10, 'color': 'black'}
    gl.ylabel_style = {'size': 10, 'color': 'black'}
    
    
    p3 = axs[1].pcolormesh(ds_binning_allscale.lon, ds_binning_allscale.lat, (1 - ds_binning_allscale['variance_mapping_err_v']/ds_binning_allscale['variance_drifter_v']), vmin=vmin, vmax=vmax, cmap='RdYlGn')
    axs[1].set_title('Meridional current [All scale]')
    axs[1].coastlines(resolution='10m', lw=0.5)
    # optional add grid lines
    p3.axes.gridlines(color='black', alpha=0., linestyle='--')
    # draw parallels/meridiens and write labels
    gl = p3.axes.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                          linewidth=0.1, color='black', alpha=0.5, linestyle='--')
    # adjust labels to taste
    gl.xlabels_top = False
    gl.ylabels_left = False
    gl.ylabels_right = False
    gl.ylocator = mticker.FixedLocator([-90, -60, -30, 0, 30, 60, 90])
    gl.xformatter = LONGITUDE_FORMATTER
    gl.yformatter = LATITUDE_FORMATTER
    gl.xlabel_style = {'size': 10, 'color': 'black'}
    gl.ylabel_style = {'size': 10, 'color': 'black'}
    
    cax = fig.add_axes([0.92, 0.42, 0.02, 0.25])
    fig.colorbar(p3, cax=cax, orientation='vertical')
    cax.set_ylabel('Explained variance', fontweight='bold')
    
    
    fig.subplots_adjust(bottom=0.2, top=0.9, left=0.1, right=0.9,
                    wspace=0.02, hspace=0.01) 
    
    plt.savefig("../figures/Maps_"+str(method_name)+"_explvar_"+region+"_uv.png", bbox_inches='tight')
    
    
def plot_psd_scores_currents_png(filename,region='glob'):    
    """
    Plot various spectral and coherence maps for zonal and meridional currents as PNG files.

    Parameters
    ----------
    filename : str
        The name of the input NetCDF file.
    region : str, optional
        The region to plot ('glob' for global, or specify a region name).

    Returns
    -------
    None
    """
    
    def coriolis_parameter(lat):
        """
        Compute the Coriolis parameter for the given latitude:
        ``f = 2*omega*sin(lat)``, where omega is the angular velocity
        of the Earth.

        Parameters
        ----------
        lat : array
        Latitude [degrees].
        """
        omega = 7.2921159e-05  # angular velocity of the Earth [rad/s]
        fc = 2 * omega * np.sin(lat * np.pi / 180.0)
        # avoid zero near equator, bound fc by min val as 1.e-8
        return np.maximum(abs(fc), 1.0e-8) * ((fc >= 0) * 2 - 1)*(24*3600)
    
    
    ds_psd = xr.open_dataset(filename)
    
    
    
    try : 
        method_name = ds_psd.attrs['method']
    except:
        method_name = ''
     
    fig, axs = plt.subplots(nrows=2,ncols=3,figsize=(11,20))

    axs=axs.flatten()
    
    vmin = -4.
    vmax= 0.
    p0 = axs[0].pcolormesh(ds_psd.wavenumber, ds_psd.lat, np.log10(ds_psd.psd_ref), vmin=vmin, vmax=vmax, cmap='Spectral_r')
    axs[0].set_xscale('symlog', linthresh=0.05)
    axs[0].set_title('Zonally averaged rotary spectra\n DRIFTERS')
    axs[0].set_ylim(-60, 60)
    axs[0].set_xlim(-2, 2)
    axs[0].set_xticklabels([-1, -0.1, '',0, '',0.1, 1])
    axs[0].set_yticklabels(['60°S', '40°S', '20°S' ,'0°', '20°N', '40°N', '60°N'])
    axs[0].set_xlabel('wavenumber [cpd]')
    axs[0].grid(alpha=0.5)
    ax2 = axs[0].twiny()
    ax2.pcolormesh(ds_psd.wavenumber, ds_psd.lat, np.log10(ds_psd.psd_ref), vmin=vmin, vmax=vmax, cmap='Spectral_r')
    ax2.plot(-coriolis_parameter(ds_psd.lat)/6, ds_psd.lat, color='k', alpha=0.5, ls='--')      # /6 because 6h drifter database
    ax2.set_xlim(-2, 2)
    ax2.set_xscale('symlog', linthresh=0.05)
    ax2.set_xticklabels(['1d', '10d', '','inf', '','10d', '1d'])
   
    
    
    p1 = axs[1].pcolormesh(ds_psd.wavenumber, ds_psd.lat, np.log10(ds_psd.psd_study), vmin=vmin, vmax=vmax, cmap='Spectral_r')
    axs[1].set_xscale('symlog', linthresh=0.05)
    axs[1].set_title('Zonally averaged rotary spectra \n MAP')
    axs[1].set_ylim(-60, 60)
    axs[1].set_xlim(-2, 2)
    axs[1].set_xticklabels([-1, -0.1, '',0, '',0.1, 1])
    axs[1].set_xlabel('wavenumber [cpd]')
    axs[1].grid(alpha=0.5)
    axs[1].set_yticks([])
    ax2 = axs[1].twiny()
    ax2.pcolormesh(ds_psd.wavenumber, ds_psd.lat, np.log10(ds_psd.psd_study), vmin=vmin, vmax=vmax, cmap='Spectral_r')
    ax2.plot(-coriolis_parameter(ds_psd.lat)/6, ds_psd.lat, color='k', alpha=0.5, ls='--')
    ax2.set_xlim(-2, 2)
    ax2.set_xscale('symlog', linthresh=0.05)
    ax2.set_xticklabels(['1d', '10d', '','inf', '','10d', '1d'])
    ax2.set_yticks([])
    
    
    p2 = axs[2].pcolormesh(ds_psd.wavenumber, ds_psd.lat, np.log10(ds_psd.psd_diff), vmin=vmin, vmax=vmax, cmap='Reds')
    axs[2].set_xscale('symlog', linthresh=0.05)
    axs[2].set_title('Zonally averaged rotary spectra \n ERR MAP-DRIFTERS')
    axs[2].set_ylim(-60, 60)
    axs[2].set_xlim(-2, 2)
    axs[2].set_xticklabels([-1, -0.1, '',0, '',0.1, 1])
    axs[2].set_xlabel('wavenumber [cpd]')
    axs[2].grid(alpha=0.5)
    axs[2].set_yticks([])
    ax2 = axs[2].twiny()
    ax2.pcolormesh(ds_psd.wavenumber, ds_psd.lat, np.log10(ds_psd.psd_diff), vmin=vmin, vmax=vmax, cmap='Reds')
    ax2.plot(-coriolis_parameter(ds_psd.lat)/6, ds_psd.lat, color='k', alpha=0.5, ls='--')
    ax2.set_xlim(-2, 2)
    ax2.set_xscale('symlog', linthresh=0.05)
    ax2.set_xticklabels(['1d', '10d', '','inf', '','10d', '1d'])
    ax2.set_yticks([])
    
    
    p3 = axs[3].pcolormesh(ds_psd.wavenumber, ds_psd.lat, ds_psd.coherence, vmin=0, vmax=1, cmap='RdYlGn')
    axs[3].set_xscale('symlog', linthresh=0.05)
    axs[3].set_title('Zonally averaged Coherence')
    axs[3].set_ylim(-60, 60)
    axs[3].set_xlim(-2, 2)
    axs[3].set_xticklabels([-1, -0.1, '',0, '',0.1, 1])
    axs[3].set_yticklabels(['60°S', '40°S', '20°S' ,'0°', '20°N', '40°N', '60°N'])
    axs[3].set_xlabel('wavenumber [cpd]')
    axs[3].grid(alpha=0.5)
    ax2 = axs[3].twiny()
    ax2.pcolormesh(ds_psd.wavenumber, ds_psd.lat, ds_psd.coherence, vmin=0, vmax=1, cmap='RdYlGn')
    ax2.contour(ds_psd.wavenumber, ds_psd.lat, ds_psd.coherence, levels=[0.5], colors='k', lws=2)
    ax2.plot(-coriolis_parameter(ds_psd.lat)/6, ds_psd.lat, color='k', alpha=0.5, ls='--')
    ax2.set_xlim(-2, 2)
    ax2.set_xscale('symlog', linthresh=0.05)
    ax2.set_xticklabels(['1d', '10d', '','inf', '','10d', '1d'])
    
    p4 = axs[4].pcolormesh(ds_psd.wavenumber, ds_psd.lat, (1. - ds_psd.psd_diff/ds_psd.psd_ref), vmin=0, vmax=1, cmap='RdYlGn')
    axs[4].set_xscale('symlog', linthresh=0.05)
    axs[4].set_title('Zonally averaged PSD$_{err}$/PSD$_{uv}$')
    axs[4].set_ylim(-60, 60)
    axs[4].set_xlim(-2, 2)
    axs[4].set_xticklabels([-1, -0.1, '',0, '',0.1, 1])
    axs[4].set_xlabel('wavenumber [cpd]')
    axs[4].grid(alpha=0.5)
    ax2 = axs[4].twiny()
    ax2.pcolormesh(ds_psd.wavenumber, ds_psd.lat, (1. - ds_psd.psd_diff/ds_psd.psd_ref), vmin=0, vmax=1, cmap='RdYlGn')
    c = ax2.contour(ds_psd.wavenumber, ds_psd.lat, (1. - ds_psd.psd_diff/ds_psd.psd_ref), levels=[0.5], colors='k', lws=2)
    ax2.plot(-coriolis_parameter(ds_psd.lat)/6, ds_psd.lat, color='k', alpha=0.5, ls='--')
    ax2.set_xlim(-2, 2)
    ax2.set_xscale('symlog', linthresh=0.05)
    ax2.set_xticklabels(['1d', '10d', '','inf', '','10d', '1d'])
    ax2.set_yticks([])
    
    
    cax = fig.add_axes([0.12, 0.56, 0.47, 0.01])
    fig.colorbar(p0, cax=cax, orientation='horizontal')
    cax.set_xlabel('log10(PSD) [m$^2$.s$^{-2}$/cpd]', fontweight='bold')
    
    cax = fig.add_axes([0.65, 0.56, 0.47*0.5, 0.01])
    fig.colorbar(p2, cax=cax, orientation='horizontal')
    cax.set_xlabel('log10(PSD$_{err}$) [m$^2$.s$^{-2}$/cpd]', fontweight='bold')
    
    cax = fig.add_axes([0.12, 0.16, 0.47, 0.01])
    cbar = fig.colorbar(p3, cax=cax, orientation='horizontal')
    cax.set_xlabel('Unresolved scales <          0.5         < Resolved scales     ', fontweight='bold')
    cbar.add_lines(c)
    
    fig.subplots_adjust(bottom=0.2, top=0.9, left=0.1, right=0.9,
                    wspace=0.02, hspace=0.3)
    
    fig.delaxes(axs[-1])
    
    
    plt.savefig("../figures/Maps_"+str(method_name)+"_effres_"+region+"_uv.png", bbox_inches='tight')
    
       
    
def plot_stat_score_map_png(filename,region='glob',box_lonlat=None, change_lon=True): 
    """
    Plot statistical score maps for zonal and meridional currents and save them as PNG files.

    Parameters
    ----------
    filename : str
        The name of the input NetCDF file.
    region : str, optional
        The region to plot ('glob' for global, or specify a region name).
    box_lonlat : dict, optional
        A dictionary with 'lon_min', 'lon_max', 'lat_min', and 'lat_max' keys
        defining a bounding box for the plot.

    Returns
    -------
    None

    Notes
    -----
    This function reads data from the input NetCDF file, creates statistical score maps for zonal and meridional currents,
    and saves the resulting plots as PNG files.

    If a bounding box is provided, the function restricts the plot to the specified region.

    The function saves two sets of plots: error variance maps and explained variance maps, both for two different scales.

    Error Variance Maps:
    - The first set of plots displays the error variance for zonal and meridional currents at two different scales (All scale and 65-500km).

    Explained Variance Maps:
    - The second set of plots displays the explained variance for zonal and meridional currents at two different scales (All scale and 65-500km).

    The PNG files are saved in the '../figures/' directory with informative filenames based on the method used.

    Examples
    --------
    >>> plot_stat_score_map_png("ocean_data.nc", region='global')
    >>> plot_stat_score_map_png("ocean_data.nc", region='Pacific', box_lonlat={'lon_min': 120, 'lon_max': 240, 'lat_min': -30, 'lat_max': 30})
    """

    ds_binning_allscale = xr.open_dataset(filename, group='all_scale')
    ds_binning_filtered = xr.open_dataset(filename, group='filtered')
     
    try : 
        method_name = ds_binning_allscale.attrs['method']
    except:
        method_name = ''
    
    if box_lonlat is not None:
        lon_min = box_lonlat['lon_min']
        lon_max = box_lonlat['lon_max']
        lat_min = box_lonlat['lat_min']
        lat_max = box_lonlat['lat_max']
        ds_binning_allscale = regional_zoom(ds_binning_allscale, [lon_min,lon_max], [lat_min,lat_max], namelon='lon', namelat='lat', change_lon=change_lon)
        ds_binning_filtered = regional_zoom(ds_binning_filtered, [lon_min,lon_max], [lat_min,lat_max], namelon='lon', namelat='lat', change_lon=change_lon)
   
    
    
    fig, axs = plt.subplots(nrows=1,ncols=2,
                        subplot_kw={'projection': ccrs.PlateCarree()},
                        figsize=(11,7.5))

    axs=axs.flatten()
     
    
    vmin = 0.
    vmax= 0.002
    p0 = axs[0].pcolormesh(ds_binning_allscale.lon, ds_binning_allscale.lat, ds_binning_allscale.variance_mapping_err, vmin=vmin, vmax=vmax, cmap='Reds')
    axs[0].set_title('SSH [All scale]')
    axs[0].coastlines(resolution='10m', lw=0.5)  
    # optional add grid lines
    p0.axes.gridlines(color='black', alpha=0., linestyle='--')
    # draw parallels/meridiens and write labels
    gl = p0.axes.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                          linewidth=0.1, color='black', alpha=0.5, linestyle='--')
    # adjust labels to taste
    gl.xlabels_top = False
    gl.ylabels_right = False
    gl.xlabels_bottom = False
    gl.ylocator = mticker.FixedLocator([-90, -60, -30, 0, 30, 60, 90])
    gl.xformatter = LONGITUDE_FORMATTER
    gl.yformatter = LATITUDE_FORMATTER
    gl.xlabel_style = {'size': 10, 'color': 'black'}
    gl.ylabel_style = {'size': 10, 'color': 'black'}
    
     
    
    p1 = axs[1].pcolormesh(ds_binning_filtered.lon, ds_binning_filtered.lat, ds_binning_filtered.variance_mapping_err, vmin=vmin, vmax=vmax, cmap='Reds')
    axs[1].set_title('SSH [65-500km]')
    axs[1].coastlines(resolution='10m', lw=0.5)
    # optional add grid lines
    p1.axes.gridlines(color='black', alpha=0., linestyle='--')
    # draw parallels/meridiens and write labels
    gl = p1.axes.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                          linewidth=0.1, color='black', alpha=0.5, linestyle='--')
    # adjust labels to taste
    gl.xlabels_top = False
    gl.ylabels_right = False
    gl.ylabels_left = False
    gl.xlabels_bottom = False
    gl.ylocator = mticker.FixedLocator([-90, -60, -30, 0, 30, 60, 90])
    gl.xformatter = LONGITUDE_FORMATTER
    gl.yformatter = LATITUDE_FORMATTER
    gl.xlabel_style = {'size': 10, 'color': 'black'}
    gl.ylabel_style = {'size': 10, 'color': 'black'}
    
    
    cax = fig.add_axes([0.92, 0.42, 0.02, 0.25])
    cbar = fig.colorbar(p1, cax=cax, orientation='vertical')
    cax.set_ylabel('Error variance [m$^2$]', fontweight='bold')
    
    plt.savefig("../figures/Maps_"+str(method_name)+"_errvar_"+region+".png", bbox_inches='tight')
    
    
    fig, axs = plt.subplots(nrows=1,ncols=2,
                        subplot_kw={'projection': ccrs.PlateCarree()},
                        figsize=(11,7.5))

    axs=axs.flatten()
    vmin = 0.
    vmax= 1
    p2 = axs[0].pcolormesh(ds_binning_allscale.lon, ds_binning_allscale.lat, (1 - ds_binning_allscale['variance_mapping_err']/ds_binning_allscale['variance_track']), vmin=vmin, vmax=vmax, cmap='RdYlGn')
    axs[0].set_title('SSH [All scale]')
    axs[0].coastlines(resolution='10m', lw=0.5)
    # optional add grid lines
    p2.axes.gridlines(color='black', alpha=0., linestyle='--')
    # draw parallels/meridiens and write labels
    gl = p2.axes.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                          linewidth=0.1, color='black', alpha=0.5, linestyle='--')
    # adjust labels to taste
    gl.xlabels_top = False
    gl.ylabels_right = False
    gl.ylocator = mticker.FixedLocator([-90, -60, -30, 0, 30, 60, 90])
    gl.xformatter = LONGITUDE_FORMATTER
    gl.yformatter = LATITUDE_FORMATTER
    gl.xlabel_style = {'size': 10, 'color': 'black'}
    gl.ylabel_style = {'size': 10, 'color': 'black'}
    
    
    p3 = axs[1].pcolormesh(ds_binning_allscale.lon, ds_binning_allscale.lat, (1 - ds_binning_filtered['variance_mapping_err']/ds_binning_filtered['variance_track']), vmin=vmin, vmax=vmax, cmap='RdYlGn')
    axs[1].set_title('SSH [65-500km]')
    axs[1].coastlines(resolution='10m', lw=0.5)
    # optional add grid lines
    p3.axes.gridlines(color='black', alpha=0., linestyle='--')
    # draw parallels/meridiens and write labels
    gl = p3.axes.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                          linewidth=0.1, color='black', alpha=0.5, linestyle='--')
    # adjust labels to taste
    gl.xlabels_top = False
    gl.ylabels_left = False
    gl.ylabels_right = False
    gl.ylocator = mticker.FixedLocator([-90, -60, -30, 0, 30, 60, 90])
    gl.xformatter = LONGITUDE_FORMATTER
    gl.yformatter = LATITUDE_FORMATTER
    gl.xlabel_style = {'size': 10, 'color': 'black'}
    gl.ylabel_style = {'size': 10, 'color': 'black'}
    
    cax = fig.add_axes([0.92, 0.42, 0.02, 0.25])
    fig.colorbar(p3, cax=cax, orientation='vertical')
    cax.set_ylabel('Explained variance', fontweight='bold')
    
    fig.subplots_adjust(bottom=0.2, top=0.9, left=0.1, right=0.9,
                    wspace=0.02, hspace=0.01)
        
    plt.savefig("../figures/Maps_"+str(method_name)+"_explvar_"+region+".png", bbox_inches='tight')
    
    
    
def compare_stat_score_map(study_filename, ref_filename):
    """
    Compare statistical score maps between a study dataset and a reference dataset and visualize the differences.

    Parameters
    ----------
    study_filename : str
        The name of the study NetCDF file.
    ref_filename : str
        The name of the reference NetCDF file.

    Returns
    -------
    holoviews.core.layout.Layout
        A layout containing Holoviews plots displaying the comparison results.

    Notes
    -----
    This function reads data from both the study and reference NetCDF files and generates various comparison plots
    between the two datasets.

    Comparison Plots:
    1. Error Variance Maps (Reference and Study) for All Scales.
    2. Error Variance Maps (Reference and Study) for 65-500km Scale.
    3. Percentage Change in Error Variance (Study vs. Reference) for All Scales.
    4. Percentage Change in Error Variance (Study vs. Reference) for 65-500km Scale.
    5. Explained Variance Maps (Reference) for All Scales.
    6. Explained Variance Maps (Reference) for 65-500km Scale.
    7. Gain(+)/Loss(-) in Explained Variance (Study vs. Reference) for All Scales.
    8. Gain(+)/Loss(-) in Explained Variance (Study vs. Reference) for 65-500km Scale.

    The generated plots are returned as a Holoviews layout for visualization.

    Examples
    --------
    >>> study_file = "study_data.nc"
    >>> reference_file = "reference_data.nc"
    >>> comparison_plots = compare_stat_score_map(study_file, reference_file)
    >>> comparison_plots.show()  # Display the comparison plots.
    """
    
    ds_ref_binning_allscale = xr.open_dataset(ref_filename, group='all_scale')
    ds_ref_binning_filtered = xr.open_dataset(ref_filename, group='filtered')
    
    explained_variance_ref_all_scale = 1 - ds_ref_binning_allscale['variance_mapping_err']/ds_ref_binning_allscale['variance_track']
    explained_variance_ref_filtered = 1 - ds_ref_binning_filtered['variance_mapping_err']/ds_ref_binning_filtered['variance_track']
    
    ds_study_binning_allscale = xr.open_dataset(study_filename, group='all_scale')
    ds_study_binning_filtered = xr.open_dataset(study_filename, group='filtered')
    
    explained_variance_study_all_scale = 1 - ds_study_binning_allscale['variance_mapping_err']/ds_study_binning_allscale['variance_track']
    explained_variance_study_filtered = 1 - ds_study_binning_filtered['variance_mapping_err']/ds_study_binning_filtered['variance_track']
    
    
    
    fig1 = ds_ref_binning_allscale['variance_mapping_err'].hvplot.quadmesh(x='lon',
                                                              y='lat',
                                                              clim=(0, 0.002),
                                                              cmap='Reds',
                                                              rasterize=True,
                                                              title='Reference Error variance [All scale]')
    
    fig2 = ds_ref_binning_filtered['variance_mapping_err'].hvplot.quadmesh(x='lon',
                                                              y='lat',
                                                              clim=(0, 0.002),
                                                              cmap='Reds',
                                                              rasterize=True,
                                                              title='Reference Error variance [65:500km]')
    
    fig3 = (100*(ds_study_binning_allscale['variance_mapping_err'] - ds_ref_binning_allscale['variance_mapping_err'])/ds_ref_binning_allscale['variance_mapping_err']).hvplot.quadmesh(x='lon',
                                                              y='lat',
                                                              clim=(-20, 20),
                                                              cmap='coolwarm',
                                                              rasterize=True,
                                                              title='Reference Error variance [All scale]')
    
    fig4 = (100*(ds_study_binning_filtered['variance_mapping_err'] - ds_ref_binning_filtered['variance_mapping_err'])/ds_ref_binning_filtered['variance_mapping_err']).hvplot.quadmesh(x='lon',
                                                              y='lat',
                                                              clim=(-20, 20),
                                                              cmap='coolwarm',
                                                              rasterize=True,
                                                              title='Reference Error variance [65:500km]')
    
    fig5 = explained_variance_ref_all_scale.hvplot.quadmesh(x='lon',
                                                              y='lat',
                                                              clim=(0, 1),
                                                              cmap='RdYlGn',
                                                              rasterize=True,
                                                              title='Reference Explained variance [All scale]')
    
    fig6 = explained_variance_ref_filtered.hvplot.quadmesh(x='lon',
                                                              y='lat',
                                                              clim=(0, 1),
                                                              cmap='RdYlGn',
                                                              rasterize=True,
                                                              title='Reference Explained variance [65:500km]')
    
    fig7 = (explained_variance_study_all_scale - explained_variance_ref_all_scale).hvplot.quadmesh(x='lon',
                                                              y='lat',
                                                              clim=(-0.2, 0.2),
                                                              cmap='coolwarm_r',
                                                              rasterize=True,
                                                              title='Gain(+)/Loss(-) Explained variance [All scale]')
    
    fig8 = (explained_variance_study_filtered - explained_variance_ref_filtered).hvplot.quadmesh(x='lon',
                                                              y='lat',
                                                              clim=(-0.2, 0.2),
                                                              cmap='coolwarm_r',
                                                              rasterize=True,
                                                              title='Gain(+)/Loss(-) Explained variance [65:500km]')
     
    
    return (fig1 + fig2 + fig3 + fig4 + fig5 + fig6 +fig7 +fig8).cols(2)



def compare_stat_score_map_png(study_filename, ref_filename):
    """
    Generate PNG plots comparing statistical score maps between study and reference datasets.

    Parameters
    ----------
    study_filename : str
        The filename of the study dataset in NetCDF format.
    ref_filename : str
        The filename of the reference dataset in NetCDF format.

    Returns
    -------
    None

    Notes
    -----
    This function generates a set of eight subplots comparing statistical score maps between a study dataset and a
    reference dataset. The comparison includes error variance, percentage change in error variance, and explained
    variance for both global and 65-500km scales. These subplots are arranged in a 4x2 grid, providing a visual
    comparison of the datasets. Colormaps are used to represent the data, and coastlines and grid lines are added to
    enhance visualization.

    Colorbars are added to indicate the color scale for each subplot. The function saves the resulting figure as a PNG
    file.

    Examples
    --------
    >>> study_file = "study_data.nc"
    >>> reference_file = "reference_data.nc"
    >>> compare_stat_score_map_png(study_file, reference_file)
    """
    
    
    ds_ref_binning_allscale = xr.open_dataset(ref_filename, group='all_scale')
    ds_ref_binning_filtered = xr.open_dataset(ref_filename, group='filtered')
    
    explained_variance_ref_all_scale = 1 - ds_ref_binning_allscale['variance_mapping_err']/ds_ref_binning_allscale['variance_track']
    explained_variance_ref_filtered = 1 - ds_ref_binning_filtered['variance_mapping_err']/ds_ref_binning_filtered['variance_track']
    
    ds_study_binning_allscale = xr.open_dataset(study_filename, group='all_scale')
    ds_study_binning_filtered = xr.open_dataset(study_filename, group='filtered')
    
    explained_variance_study_all_scale = 1 - ds_study_binning_allscale['variance_mapping_err']/ds_study_binning_allscale['variance_track']
    explained_variance_study_filtered = 1 - ds_study_binning_filtered['variance_mapping_err']/ds_study_binning_filtered['variance_track']
    
    
    fig, axs = plt.subplots(nrows=4,ncols=2,
                        subplot_kw={'projection': ccrs.PlateCarree()},
                        figsize=(11,15))

    axs=axs.flatten()
    
    vmin = 0.
    vmax= 0.002
    p0 = axs[0].pcolormesh(ds_ref_binning_allscale.lon, ds_ref_binning_allscale.lat, ds_ref_binning_allscale.variance_mapping_err, vmin=vmin, vmax=vmax, cmap='Reds')
    axs[0].set_title('SSH [All scale]')
    axs[0].coastlines(resolution='10m', lw=0.5)
    # optional add grid lines
    p0.axes.gridlines(color='black', alpha=0., linestyle='--')
    # draw parallels/meridiens and write labels
    gl = p0.axes.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                          linewidth=0.1, color='black', alpha=0.5, linestyle='--')
    # adjust labels to taste
    gl.xlabels_top = False
    gl.ylabels_right = False
    gl.xlabels_bottom = False
    gl.ylocator = mticker.FixedLocator([-90, -60, -30, 0, 30, 60, 90])
    gl.xformatter = LONGITUDE_FORMATTER
    gl.yformatter = LATITUDE_FORMATTER
    gl.xlabel_style = {'size': 10, 'color': 'black'}
    gl.ylabel_style = {'size': 10, 'color': 'black'}
    
    p1 = axs[1].pcolormesh(ds_ref_binning_filtered.lon, ds_ref_binning_filtered.lat, ds_ref_binning_filtered.variance_mapping_err, vmin=vmin, vmax=vmax, cmap='Reds')
    axs[1].set_title('SSH [65-500km]')
    axs[1].coastlines(resolution='10m', lw=0.5)
    # optional add grid lines
    p1.axes.gridlines(color='black', alpha=0., linestyle='--')
    # draw parallels/meridiens and write labels
    gl = p1.axes.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                          linewidth=0.1, color='black', alpha=0.5, linestyle='--')
    # adjust labels to taste
    gl.xlabels_top = False
    gl.ylabels_right = False
    gl.ylabels_left = False
    gl.xlabels_bottom = False
    gl.ylocator = mticker.FixedLocator([-90, -60, -30, 0, 30, 60, 90])
    gl.xformatter = LONGITUDE_FORMATTER
    gl.yformatter = LATITUDE_FORMATTER
    gl.xlabel_style = {'size': 10, 'color': 'black'}
    gl.ylabel_style = {'size': 10, 'color': 'black'}
    
    vmin = -20.
    vmax= 20
    p2 = axs[2].pcolormesh(ds_ref_binning_allscale.lon, ds_ref_binning_allscale.lat, 100*(ds_study_binning_allscale.variance_mapping_err - ds_ref_binning_allscale.variance_mapping_err)/ds_ref_binning_allscale.variance_mapping_err, vmin=vmin, vmax=vmax, cmap='coolwarm')
    axs[2].set_title('SSH [All scale]')
    axs[2].coastlines(resolution='10m', lw=0.5)
    # optional add grid lines
    p2.axes.gridlines(color='black', alpha=0., linestyle='--')
    # draw parallels/meridiens and write labels
    gl = p2.axes.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                          linewidth=0.1, color='black', alpha=0.5, linestyle='--')
    # adjust labels to taste
    gl.xlabels_top = False
    gl.ylabels_right = False
    gl.ylocator = mticker.FixedLocator([-90, -60, -30, 0, 30, 60, 90])
    gl.xformatter = LONGITUDE_FORMATTER
    gl.yformatter = LATITUDE_FORMATTER
    gl.xlabel_style = {'size': 10, 'color': 'black'}
    gl.ylabel_style = {'size': 10, 'color': 'black'}
    
    
    p3 = axs[3].pcolormesh(ds_study_binning_filtered.lon, 
                           ds_study_binning_filtered.lat, 
                           100*(ds_study_binning_filtered.variance_mapping_err - ds_ref_binning_filtered.variance_mapping_err)/ds_ref_binning_filtered.variance_mapping_err, 
                           vmin=vmin, 
                           vmax=vmax, 
                           cmap='coolwarm')
    axs[3].set_title('SSH [65-500km]')
    axs[3].coastlines(resolution='10m', lw=0.5)
    # optional add grid lines
    p3.axes.gridlines(color='black', alpha=0., linestyle='--')
    # draw parallels/meridiens and write labels
    gl = p3.axes.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                          linewidth=0.1, color='black', alpha=0.5, linestyle='--')
    # adjust labels to taste
    gl.xlabels_top = False
    gl.ylabels_left = False
    gl.ylabels_right = False
    gl.ylocator = mticker.FixedLocator([-90, -60, -30, 0, 30, 60, 90])
    gl.xformatter = LONGITUDE_FORMATTER
    gl.yformatter = LATITUDE_FORMATTER
    gl.xlabel_style = {'size': 10, 'color': 'black'}
    gl.ylabel_style = {'size': 10, 'color': 'black'}
    
    
    
    vmin = 0.
    vmax= 1
    p4 = axs[4].pcolormesh(ds_ref_binning_allscale.lon, ds_ref_binning_allscale.lat, explained_variance_ref_all_scale, vmin=vmin, vmax=vmax, cmap='RdYlGn')
    axs[4].set_title('SSH [All scale]')
    axs[4].coastlines(resolution='10m', lw=0.5)
    # optional add grid lines
    p4.axes.gridlines(color='black', alpha=0., linestyle='--')
    # draw parallels/meridiens and write labels
    gl = p4.axes.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                          linewidth=0.1, color='black', alpha=0.5, linestyle='--')
    # adjust labels to taste
    gl.xlabels_top = False
    gl.ylabels_right = False
    gl.xlabels_bottom = False
    gl.ylocator = mticker.FixedLocator([-90, -60, -30, 0, 30, 60, 90])
    gl.xformatter = LONGITUDE_FORMATTER
    gl.yformatter = LATITUDE_FORMATTER
    gl.xlabel_style = {'size': 10, 'color': 'black'}
    gl.ylabel_style = {'size': 10, 'color': 'black'}
    
    p5 = axs[5].pcolormesh(ds_ref_binning_filtered.lon, ds_ref_binning_filtered.lat, explained_variance_ref_filtered, vmin=vmin, vmax=vmax, cmap='RdYlGn')
    axs[5].set_title('SSH [65-500km]')
    axs[5].coastlines(resolution='10m', lw=0.5)
    # optional add grid lines
    p5.axes.gridlines(color='black', alpha=0., linestyle='--')
    # draw parallels/meridiens and write labels
    gl = p5.axes.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                          linewidth=0.1, color='black', alpha=0.5, linestyle='--')
    # adjust labels to taste
    gl.xlabels_top = False
    gl.ylabels_right = False
    gl.ylabels_left = False
    gl.xlabels_bottom = False
    gl.ylocator = mticker.FixedLocator([-90, -60, -30, 0, 30, 60, 90])
    gl.xformatter = LONGITUDE_FORMATTER
    gl.yformatter = LATITUDE_FORMATTER
    gl.xlabel_style = {'size': 10, 'color': 'black'}
    gl.ylabel_style = {'size': 10, 'color': 'black'}
    
    vmin = -0.2
    vmax= 0.2
    p6 = axs[6].pcolormesh(ds_ref_binning_allscale.lon, 
                           ds_ref_binning_allscale.lat, 
                           (explained_variance_study_all_scale - explained_variance_ref_all_scale), 
                           vmin=vmin, 
                           vmax=vmax, 
                           cmap='coolwarm_r')
    axs[6].set_title('SSH [All scale]')
    axs[6].coastlines(resolution='10m', lw=0.5)
    # optional add grid lines
    p6.axes.gridlines(color='black', alpha=0., linestyle='--')
    # draw parallels/meridiens and write labels
    gl = p6.axes.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                          linewidth=0.1, color='black', alpha=0.5, linestyle='--')
    # adjust labels to taste
    gl.xlabels_top = False
    gl.ylabels_right = False
    gl.ylocator = mticker.FixedLocator([-90, -60, -30, 0, 30, 60, 90])
    gl.xformatter = LONGITUDE_FORMATTER
    gl.yformatter = LATITUDE_FORMATTER
    gl.xlabel_style = {'size': 10, 'color': 'black'}
    gl.ylabel_style = {'size': 10, 'color': 'black'}
    
    
    p7 = axs[7].pcolormesh(ds_study_binning_filtered.lon, 
                           ds_study_binning_filtered.lat, 
                           (explained_variance_study_filtered - explained_variance_ref_filtered), 
                           vmin=vmin, 
                           vmax=vmax, 
                           cmap='coolwarm_r')
    axs[7].set_title('SSH [65-500km]')
    axs[7].coastlines(resolution='10m', lw=0.5)
    # optional add grid lines
    p7.axes.gridlines(color='black', alpha=0., linestyle='--')
    # draw parallels/meridiens and write labels
    gl = p7.axes.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                          linewidth=0.1, color='black', alpha=0.5, linestyle='--')
    # adjust labels to taste
    gl.xlabels_top = False
    gl.ylabels_left = False
    gl.ylabels_right = False
    gl.ylocator = mticker.FixedLocator([-90, -60, -30, 0, 30, 60, 90])
    gl.xformatter = LONGITUDE_FORMATTER
    gl.yformatter = LATITUDE_FORMATTER
    gl.xlabel_style = {'size': 10, 'color': 'black'}
    gl.ylabel_style = {'size': 10, 'color': 'black'}
    
    
    
    
    
    cax = fig.add_axes([0.92, 0.57, 0.02, 0.13])
    fig.colorbar(p3, cax=cax, orientation='vertical')
    cax.set_ylabel('Loss(-)/Gain(+)\n Error variance [%]', fontweight='bold')
    
    cax = fig.add_axes([0.92, 0.75, 0.02, 0.13])
    cbar = fig.colorbar(p1, cax=cax, orientation='vertical')
    cax.set_ylabel('Error variance [m$^2$]', fontweight='bold')
    
    cax = fig.add_axes([0.92, 0.22, 0.02, 0.13])
    fig.colorbar(p7, cax=cax, orientation='vertical')
    cax.set_ylabel('Loss(-)/Gain(+)\n Explained variance [%]', fontweight='bold')
    
    cax = fig.add_axes([0.92, 0.4, 0.02, 0.13])
    cbar = fig.colorbar(p5, cax=cax, orientation='vertical')
    cax.set_ylabel('Explained variance', fontweight='bold')
    
    fig.subplots_adjust(bottom=0.2, top=0.9, left=0.1, right=0.9,
                    wspace=0.02, hspace=0.01)
    
    
def compare_stat_score_map_uv_png(study_filename, ref_filename):
    """
    Generate PNG plots comparing statistical score maps for zonal and meridional currents between study and reference datasets.

    Parameters
    ----------
    study_filename : str
        The filename of the study dataset in NetCDF format.
    ref_filename : str
        The filename of the reference dataset in NetCDF format.

    Returns
    -------
    None

    Notes
    -----
    This function generates a set of eight subplots comparing statistical score maps for zonal and meridional currents
    between a study dataset and a reference dataset. The comparison includes error variance, percentage change in error
    variance, and explained variance for both zonal and meridional currents at all scales. These subplots are arranged
    in a 4x2 grid, providing a visual comparison of the datasets. Colormaps are used to represent the data, and coastlines
    and grid lines are added to enhance visualization.

    Colorbars are added to indicate the color scale for each subplot. The function saves the resulting figure as a PNG
    file.

    Examples
    --------
    >>> study_file = "study_data.nc"
    >>> reference_file = "reference_data.nc"
    >>> compare_stat_score_map_uv_png(study_file, reference_file)
    """

    ds_ref_binning_allscale = xr.open_dataset(ref_filename, group='all_scale')
    
    explained_variance_u_ref_all_scale = 1 - ds_ref_binning_allscale['variance_mapping_err_u']/ds_ref_binning_allscale['variance_drifter_u']
    explained_variance_v_ref_all_scale = 1 - ds_ref_binning_allscale['variance_mapping_err_v']/ds_ref_binning_allscale['variance_drifter_v']
    
    ds_study_binning_allscale = xr.open_dataset(study_filename, group='all_scale')
    
    explained_variance_u_study_all_scale = 1 - ds_study_binning_allscale['variance_mapping_err_u']/ds_study_binning_allscale['variance_drifter_u']
    explained_variance_v_study_all_scale = 1 - ds_study_binning_allscale['variance_mapping_err_v']/ds_study_binning_allscale['variance_drifter_v']
    
    
    fig, axs = plt.subplots(nrows=4,ncols=2,
                        subplot_kw={'projection': ccrs.PlateCarree()},
                        figsize=(11,15))

    axs=axs.flatten()
    
    vmin = 0.
    vmax= 0.1
    p0 = axs[0].pcolormesh(ds_ref_binning_allscale.lon, ds_ref_binning_allscale.lat, ds_ref_binning_allscale.variance_mapping_err_u, vmin=vmin, vmax=vmax, cmap='Reds')
    axs[0].set_title('Zonal current [All scale]')
    axs[0].coastlines(resolution='10m', lw=0.5)
    # optional add grid lines
    p0.axes.gridlines(color='black', alpha=0., linestyle='--')
    # draw parallels/meridiens and write labels
    gl = p0.axes.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                          linewidth=0.1, color='black', alpha=0.5, linestyle='--')
    # adjust labels to taste
    gl.xlabels_top = False
    gl.ylabels_right = False
    gl.xlabels_bottom = False
    gl.ylocator = mticker.FixedLocator([-90, -60, -30, 0, 30, 60, 90])
    gl.xformatter = LONGITUDE_FORMATTER
    gl.yformatter = LATITUDE_FORMATTER
    gl.xlabel_style = {'size': 10, 'color': 'black'}
    gl.ylabel_style = {'size': 10, 'color': 'black'}
    
    p1 = axs[1].pcolormesh(ds_ref_binning_allscale.lon, ds_ref_binning_allscale.lat, ds_ref_binning_allscale.variance_mapping_err_v, vmin=vmin, vmax=vmax, cmap='Reds')
    axs[1].set_title('Meridional current [All scale]')
    axs[1].coastlines(resolution='10m', lw=0.5)
    # optional add grid lines
    p1.axes.gridlines(color='black', alpha=0., linestyle='--')
    # draw parallels/meridiens and write labels
    gl = p1.axes.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                          linewidth=0.1, color='black', alpha=0.5, linestyle='--')
    # adjust labels to taste
    gl.xlabels_top = False
    gl.ylabels_right = False
    gl.ylabels_left = False
    gl.xlabels_bottom = False
    gl.ylocator = mticker.FixedLocator([-90, -60, -30, 0, 30, 60, 90])
    gl.xformatter = LONGITUDE_FORMATTER
    gl.yformatter = LATITUDE_FORMATTER
    gl.xlabel_style = {'size': 10, 'color': 'black'}
    gl.ylabel_style = {'size': 10, 'color': 'black'}
    
    vmin = -20.
    vmax= 20
    p2 = axs[2].pcolormesh(ds_ref_binning_allscale.lon, 
                           ds_ref_binning_allscale.lat, 
                           100*(ds_study_binning_allscale.variance_mapping_err_u - ds_ref_binning_allscale.variance_mapping_err_u)/ds_ref_binning_allscale.variance_mapping_err_u, 
                           vmin=vmin, 
                           vmax=vmax, 
                           cmap='coolwarm')
    axs[2].set_title('Zonal current [All scale]')
    axs[2].coastlines(resolution='10m', lw=0.5)
    # optional add grid lines
    p2.axes.gridlines(color='black', alpha=0., linestyle='--')
    # draw parallels/meridiens and write labels
    gl = p2.axes.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                          linewidth=0.1, color='black', alpha=0.5, linestyle='--')
    # adjust labels to taste
    gl.xlabels_top = False
    gl.ylabels_right = False
    gl.ylocator = mticker.FixedLocator([-90, -60, -30, 0, 30, 60, 90])
    gl.xformatter = LONGITUDE_FORMATTER
    gl.yformatter = LATITUDE_FORMATTER
    gl.xlabel_style = {'size': 10, 'color': 'black'}
    gl.ylabel_style = {'size': 10, 'color': 'black'}
    
    
    p3 = axs[3].pcolormesh(ds_ref_binning_allscale.lon, 
                           ds_ref_binning_allscale.lat, 
                           100*(ds_study_binning_allscale.variance_mapping_err_v - ds_ref_binning_allscale.variance_mapping_err_v)/ds_ref_binning_allscale.variance_mapping_err_v, 
                           vmin=vmin, 
                           vmax=vmax, 
                           cmap='coolwarm')
    axs[3].set_title('Meridional current [All scale]')
    axs[3].coastlines(resolution='10m', lw=0.5)
    # optional add grid lines
    p3.axes.gridlines(color='black', alpha=0., linestyle='--')
    # draw parallels/meridiens and write labels
    gl = p3.axes.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                          linewidth=0.1, color='black', alpha=0.5, linestyle='--')
    # adjust labels to taste
    gl.xlabels_top = False
    gl.ylabels_left = False
    gl.ylabels_right = False
    gl.ylocator = mticker.FixedLocator([-90, -60, -30, 0, 30, 60, 90])
    gl.xformatter = LONGITUDE_FORMATTER
    gl.yformatter = LATITUDE_FORMATTER
    gl.xlabel_style = {'size': 10, 'color': 'black'}
    gl.ylabel_style = {'size': 10, 'color': 'black'}
    
    
    
    vmin = 0.
    vmax= 1
    p4 = axs[4].pcolormesh(ds_ref_binning_allscale.lon, ds_ref_binning_allscale.lat, explained_variance_u_ref_all_scale, vmin=vmin, vmax=vmax, cmap='RdYlGn')
    axs[4].set_title('Zonal current [All scale]')
    axs[4].coastlines(resolution='10m', lw=0.5)
    # optional add grid lines
    p4.axes.gridlines(color='black', alpha=0., linestyle='--')
    # draw parallels/meridiens and write labels
    gl = p4.axes.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                          linewidth=0.1, color='black', alpha=0.5, linestyle='--')
    # adjust labels to taste
    gl.xlabels_top = False
    gl.ylabels_right = False
    gl.xlabels_bottom = False
    gl.ylocator = mticker.FixedLocator([-90, -60, -30, 0, 30, 60, 90])
    gl.xformatter = LONGITUDE_FORMATTER
    gl.yformatter = LATITUDE_FORMATTER
    gl.xlabel_style = {'size': 10, 'color': 'black'}
    gl.ylabel_style = {'size': 10, 'color': 'black'}
    
    p5 = axs[5].pcolormesh(ds_ref_binning_allscale.lon, ds_ref_binning_allscale.lat, explained_variance_v_ref_all_scale, vmin=vmin, vmax=vmax, cmap='RdYlGn')
    axs[5].set_title('Meridional current [All scale]')
    axs[5].coastlines(resolution='10m', lw=0.5)
    # optional add grid lines
    p5.axes.gridlines(color='black', alpha=0., linestyle='--')
    # draw parallels/meridiens and write labels
    gl = p5.axes.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                          linewidth=0.1, color='black', alpha=0.5, linestyle='--')
    # adjust labels to taste
    gl.xlabels_top = False
    gl.ylabels_right = False
    gl.ylabels_left = False
    gl.xlabels_bottom = False
    gl.ylocator = mticker.FixedLocator([-90, -60, -30, 0, 30, 60, 90])
    gl.xformatter = LONGITUDE_FORMATTER
    gl.yformatter = LATITUDE_FORMATTER
    gl.xlabel_style = {'size': 10, 'color': 'black'}
    gl.ylabel_style = {'size': 10, 'color': 'black'}
    
    vmin = -0.2
    vmax= 0.2
    p6 = axs[6].pcolormesh(ds_ref_binning_allscale.lon, 
                           ds_ref_binning_allscale.lat, 
                           (explained_variance_u_study_all_scale - explained_variance_u_ref_all_scale), 
                           vmin=vmin, 
                           vmax=vmax, 
                           cmap='coolwarm_r')
    axs[6].set_title('Zonal current [All scale]')
    axs[6].coastlines(resolution='10m', lw=0.5)
    # optional add grid lines
    p6.axes.gridlines(color='black', alpha=0., linestyle='--')
    # draw parallels/meridiens and write labels
    gl = p6.axes.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                          linewidth=0.1, color='black', alpha=0.5, linestyle='--')
    # adjust labels to taste
    gl.xlabels_top = False
    gl.ylabels_right = False
    gl.ylocator = mticker.FixedLocator([-90, -60, -30, 0, 30, 60, 90])
    gl.xformatter = LONGITUDE_FORMATTER
    gl.yformatter = LATITUDE_FORMATTER
    gl.xlabel_style = {'size': 10, 'color': 'black'}
    gl.ylabel_style = {'size': 10, 'color': 'black'}
    
    
    p7 = axs[7].pcolormesh(ds_ref_binning_allscale.lon, 
                           ds_ref_binning_allscale.lat, 
                           (explained_variance_v_study_all_scale - explained_variance_v_ref_all_scale), 
                           vmin=vmin, 
                           vmax=vmax, 
                           cmap='coolwarm_r')
    axs[7].set_title('Meridional current [All scale]')
    axs[7].coastlines(resolution='10m', lw=0.5)
    # optional add grid lines
    p7.axes.gridlines(color='black', alpha=0., linestyle='--')
    # draw parallels/meridiens and write labels
    gl = p7.axes.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                          linewidth=0.1, color='black', alpha=0.5, linestyle='--')
    # adjust labels to taste
    gl.xlabels_top = False
    gl.ylabels_left = False
    gl.ylabels_right = False
    gl.ylocator = mticker.FixedLocator([-90, -60, -30, 0, 30, 60, 90])
    gl.xformatter = LONGITUDE_FORMATTER
    gl.yformatter = LATITUDE_FORMATTER
    gl.xlabel_style = {'size': 10, 'color': 'black'}
    gl.ylabel_style = {'size': 10, 'color': 'black'}
    
    
    
    
    
    cax = fig.add_axes([0.92, 0.57, 0.02, 0.13])
    fig.colorbar(p3, cax=cax, orientation='vertical')
    cax.set_ylabel('Loss(-)/Gain(+)\n Error variance [%]', fontweight='bold')
    
    cax = fig.add_axes([0.92, 0.75, 0.02, 0.13])
    cbar = fig.colorbar(p1, cax=cax, orientation='vertical')
    cax.set_ylabel('Error variance [m$^2$]', fontweight='bold')
    
    cax = fig.add_axes([0.92, 0.22, 0.02, 0.13])
    fig.colorbar(p7, cax=cax, orientation='vertical')
    cax.set_ylabel('Loss(-)/Gain(+)\n Explained variance [%]', fontweight='bold')
    
    cax = fig.add_axes([0.92, 0.4, 0.02, 0.13])
    cbar = fig.colorbar(p5, cax=cax, orientation='vertical')
    cax.set_ylabel('Explained variance', fontweight='bold')
    
    fig.subplots_adjust(bottom=0.2, top=0.9, left=0.1, right=0.9,
                    wspace=0.02, hspace=0.01)



def compare_psd_score(study_filename, ref_filename):
    """
    Compare Power Spectral Density (PSD) scores between a study dataset and a reference dataset.

    Parameters
    ----------
    study_filename : str
        The filename of the study dataset in NetCDF format.
    ref_filename : str
        The filename of the reference dataset in NetCDF format.

    Returns
    -------
    panel : hvplot.core.panel.Panel
        A panel containing two subplots comparing PSD scores.

    Notes
    -----
    This function compares PSD scores between a study dataset and a reference dataset for the effective resolution.
    It generates two subplots: one showing the effective resolution of the reference dataset and the other showing the
    percentage change in effective resolution between the study and reference datasets. The subplots use colormaps to
    represent the data, and coastlines and grid lines are added for context.

    Examples
    --------
    >>> study_file = "study_data.nc"
    >>> reference_file = "reference_data.nc"
    >>> compare_psd_score(study_file, reference_file)
    """
    
    ds_ref = xr.open_dataset(ref_filename)
    ds_study = xr.open_dataset(study_filename)
    
    fig0 = ds_ref.effective_resolution.hvplot.quadmesh(x='lon', 
                                                   y='lat', 
                                                   cmap='Spectral_r', 
                                                   clim=(100, 500), 
                                                   title='Effective resolution [km]',
                                                   rasterize=True, 
                                                   projection=ccrs.PlateCarree(), 
                                                   project=True, 
                                                   geo=True, 
                                                   coastline=True)
    
    fig1 = (100*(ds_study.effective_resolution - ds_ref.effective_resolution)/ds_ref.effective_resolution).hvplot.quadmesh(x='lon', 
                                                   y='lat', 
                                                   cmap='coolwarm', 
                                                   clim=(-20, 20), 
                                                   title='Gain(-)/loss(+) Effective resolution [%]',
                                                   rasterize=True, 
                                                   projection=ccrs.PlateCarree(), 
                                                   project=True, 
                                                   geo=True, 
                                                   coastline=True)
    
    return (fig0 + fig1).cols(1)
    
    
def compare_psd_score_png(study_filename, ref_filename):
    """
    Generate PNG plots comparing Power Spectral Density (PSD) scores for effective resolution between a study dataset and a reference dataset.

    Parameters
    ----------
    study_filename : str
        The filename of the study dataset in NetCDF format.
    ref_filename : str
        The filename of the reference dataset in NetCDF format.

    Returns
    -------
    None

    Notes
    -----
    This function generates two subplots comparing PSD scores for effective resolution between a study dataset and a
    reference dataset. The subplots display the effective resolution of the reference dataset and the percentage change
    in effective resolution between the study and reference datasets. Colormaps are used to represent the data, and
    coastlines, land features, and grid lines are added to enhance visualization.

    Colorbars are added to indicate the color scale for each subplot. The resulting figure is saved as a PNG file.

    Examples
    --------
    >>> study_file = "study_data.nc"
    >>> reference_file = "reference_data.nc"
    >>> compare_psd_score_png(study_file, reference_file)
    """
    
    ds_ref = xr.open_dataset(ref_filename)
    ds_study = xr.open_dataset(study_filename)
    
    fig, axs = plt.subplots(nrows=3,ncols=1,
                        subplot_kw={'projection': ccrs.PlateCarree()},
                        figsize=(6,10))

    axs=axs.flatten()
    
    vmin = 100.
    vmax= 500.
    p0 = axs[0].pcolormesh(ds_ref.lon, ds_ref.lat, ds_ref.effective_resolution, vmin=vmin, vmax=vmax, cmap='Spectral_r')
    axs[0].set_title('SSH ')
    axs[0].coastlines(resolution='10m', lw=0.5, zorder=13)
    axs[0].add_feature(cfeature.LAND, color='w', zorder=12)
    # optional add grid lines
    p0.axes.gridlines(color='black', alpha=0., linestyle='--')
    # draw parallels/meridiens and write labels
    gl = p0.axes.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                          linewidth=0.1, color='black', alpha=0.5, linestyle='--')
    # adjust labels to taste
    gl.xlabels_top = False
    gl.ylabels_right = False
    gl.xlabels_bottom = False
    gl.ylocator = mticker.FixedLocator([-90, -60, -30, 0, 30, 60, 90])
    gl.xformatter = LONGITUDE_FORMATTER
    gl.yformatter = LATITUDE_FORMATTER
    gl.xlabel_style = {'size': 10, 'color': 'black'}
    gl.ylabel_style = {'size': 10, 'color': 'black'}
    
    vmin = -40.
    vmax= 40.
    p1 = axs[1].pcolormesh(ds_ref.lon, ds_ref.lat, ds_study.effective_resolution - ds_ref.effective_resolution, vmin=vmin, vmax=vmax, cmap='coolwarm')
    axs[1].set_title('SSH')
    axs[1].coastlines(resolution='10m', lw=0.5, zorder=13)
    axs[1].add_feature(cfeature.LAND, color='w', zorder=12)
    # optional add grid lines
    p1.axes.gridlines(color='black', alpha=0., linestyle='--')
    # draw parallels/meridiens and write labels
    gl = p1.axes.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                          linewidth=0.1, color='black', alpha=0.5, linestyle='--')
    # adjust labels to taste
    gl.xlabels_top = False
    gl.ylabels_right = False
    gl.ylabels_left = True
    gl.xlabels_bottom = False
    gl.ylocator = mticker.FixedLocator([-90, -60, -30, 0, 30, 60, 90])
    gl.xformatter = LONGITUDE_FORMATTER
    gl.yformatter = LATITUDE_FORMATTER
    gl.xlabel_style = {'size': 10, 'color': 'black'}
    gl.ylabel_style = {'size': 10, 'color': 'black'}
    
    vmin = -10.
    vmax= 10
    p2 = axs[2].pcolormesh(ds_ref.lon, ds_ref.lat, 100*(ds_study.effective_resolution - ds_ref.effective_resolution)/ds_ref.effective_resolution, vmin=vmin, vmax=vmax, cmap='coolwarm')
    axs[2].set_title('SSH')
    axs[2].coastlines(resolution='10m', lw=0.5, zorder=13)
    axs[2].add_feature(cfeature.LAND, color='w', zorder=12)
    # optional add grid lines
    p2.axes.gridlines(color='black', alpha=0., linestyle='--')
    # draw parallels/meridiens and write labels
    gl = p2.axes.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                          linewidth=0.1, color='black', alpha=0.5, linestyle='--')
    # adjust labels to taste
    gl.xlabels_top = False
    gl.ylabels_right = False
    gl.ylocator = mticker.FixedLocator([-90, -60, -30, 0, 30, 60, 90])
    gl.xformatter = LONGITUDE_FORMATTER
    gl.yformatter = LATITUDE_FORMATTER
    gl.xlabel_style = {'size': 10, 'color': 'black'}
    gl.ylabel_style = {'size': 10, 'color': 'black'}
    
    cax = fig.add_axes([0.92, 0.12, 0.02, 0.2])
    fig.colorbar(p2, cax=cax, orientation='vertical')
    cax.set_ylabel('Gain(-)/Loss(+)\n Effective resolution [%]', fontweight='bold')
    
    cax = fig.add_axes([0.92, 0.4, 0.02, 0.2])
    cbar = fig.colorbar(p1, cax=cax, orientation='vertical')
    cax.set_ylabel('Gain(-)/Loss(+)\n Effective resolution [km]', fontweight='bold')
    
    cax = fig.add_axes([0.92, 0.67, 0.02, 0.2])
    fig.colorbar(p0, cax=cax, orientation='vertical')
    cax.set_ylabel('Efective resolition [km]', fontweight='bold')
    
    
    
    

def movie(ds, name_var, method='DUACS', region='Global', dir_output='../results/',dim_name=['time','latitude','longitude'],
          framerate=24,Display=True,clim=None,cmap='Spectral', newmovie=False):

    # For memory leak when saving multiple png files...
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.gridspec as gridspec
    import gc
    import os
    import subprocess
    from IPython.display import Video
    
    
    moviename = f'movie_{method}_{region}_{name_var}.mp4'
    
    if not os.path.isfile(os.path.join(dir_output, moviename)) or  newmovie: 
     
        # Create merged dataset 
        name_dim_rmse = (dim_name[0],)
        coords = (dim_name[0],dim_name[1],dim_name[2])

        ds = ds.chunk({dim_name[0]:1})

        if name_var == 'uv':
            ds = ds.assign(uv=1/2*np.sqrt(ds["ugos"]**2+ds["vgos"]**2))
            ds.uv.attrs = ds.vgos.attrs
            ds.uv.attrs['long_name'] = 'Absolute geostrophic velocity magnitud'
            ds.uv.attrs['standard_name'] = 'surface_geostrophic_sea_water_velocity_magnitud'

        if clim is None:
            clim = (ds[name_var].to_dataset().apply(np.nanmin)[name_var].values,ds[name_var].to_dataset().apply(np.nanmax)[name_var].values) 

        # Plotting function
        def _save_single_frame(ds, tt,clim=clim,cmap=cmap):

            if tt==0:
                return

            fig = plt.figure(figsize=(8,6))
            gs = gridspec.GridSpec(1,1,width_ratios=(1,))

            date = str(ds[dim_name[0]][tt].values)[:13]

            ids = ds.isel(time=tt) # dim_name[0] ?

            ax1 = fig.add_subplot(gs[0, 0])
            
            extend_style = 'both'
            if clim[0] == 0:
                extend_style = 'max'
            if clim[1] == 0:
                extend_style = 'min'
            ids[name_var].plot(ax=ax1,cmap=cmap,vmin=clim[0],vmax=clim[1],add_colorbar=True, extend=extend_style)
            ax1.set_title(method+': '+date)  

            fig.savefig(f'{dir_output}/frame_{method}_{region}_{name_var}_{str(tt).zfill(5)}.png',dpi=100)

            plt.close(fig)
            del fig
            gc.collect(2)


        # Compute and save frames 
        for tt in range(ds[dim_name[0]].size):
            _save_single_frame(ds, tt)

        # Create movie
        sourcefolder = dir_output
        frame_pattern = f'frame_{method}_{region}_{name_var}_*.png'
        ffmpeg_options="-c:v libx264 -preset veryslow -crf 15 -pix_fmt yuv420p -loglevel quiet"

        command = 'ffmpeg -f image2 -r %i -pattern_type glob -i %s -y %s -r %i %s' % (
                framerate,
                os.path.join(sourcefolder, frame_pattern),
                ffmpeg_options,
                framerate,
                os.path.join(dir_output, moviename),
            )
        #print(command)

        _ = subprocess.run(command.split(' '),stdout=subprocess.DEVNULL)
         

        # Delete frames
        os.system(f'rm {os.path.join(sourcefolder, frame_pattern)}')
        
    else:  
        print("Movie already exists. Playing the existing movie (if Display==True). To rewrite movie set newmovie argument to True.") 

    # Display movie
    if Display: 
        return Video(os.path.join(dir_output, moviename),embed=True)
