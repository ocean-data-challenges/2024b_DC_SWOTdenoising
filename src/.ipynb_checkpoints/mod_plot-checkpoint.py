import hvplot.xarray
import xarray as xr
import cartopy.crs as ccrs
import pandas as pd


def plot_stat_score_map(filename):
    
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
    


def plot_effective_resolution(filename):
    
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


def plot_psd_scores(filename):
    
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
    
    
    