
import numpy as np 
import xarray as xr 
import matplotlib.pyplot as plt 

def plot_snapshots(ds_SWOT, variable='ssh', name_denoised='ssha_denoised_unet', date_plot=np.datetime64('2023-09-12'), 
                   region_info=None, method = 'Unet_baseline', orbit = '1d', dir_output='../figures/', savefig = True, fig=None):

    if variable == 'ssh':
        cmap = 'Spectral'
        vmin=-0.6
        vmax =0.6
        noisy_var = 'ssha_new_editing'
        
    if variable == 'grad':
        cmap = 'inferno'
        vmin= 0.
        vmax =0.0005
        noisy_var = 'grad_new_editing'
        
    if variable == 'lapl':
        cmap = 'viridis'
        vmin = -0.005
        vmax = 0.005
        noisy_var = 'lapl_new_editing'
        
    if fig == None:     
        fig = plt.figure(figsize=(14,6)) 
    plt.suptitle(date_plot,fontsize=18)

    ax = plt.subplot(121)
    plt.title('Noisy '+variable,fontsize=16)
    plt.scatter(ds_SWOT.longitude,ds_SWOT.latitude,c=ds_SWOT[noisy_var],s=2,cmap=cmap, vmin=vmin, vmax = vmax) 
    plt.colorbar()

    ax.set_xlim(region_info['lon_min'], region_info['lon_max'])
    ax.set_ylim(region_info['lat_min'], region_info['lat_max'])
    ax.set_xticks(range(region_info['lon_min'], region_info['lon_max'], 2))
    ax.set_xticklabels(region_info['lon_ticks'],fontsize=12)
    ax.set_yticks(range(region_info['lat_min'], region_info['lat_max'], 2))
    ax.set_yticklabels(region_info['lat_ticks'],fontsize=12)
    ax.set_xlabel('Longitude',fontsize=16)
    ax.set_ylabel('Latitude',fontsize=16)

    ax = plt.subplot(122)
    plt.title(method+' denoised '+variable,fontsize=16)
    plt.scatter(ds_SWOT.longitude,ds_SWOT.latitude,c=ds_SWOT[name_denoised],s=2,cmap=cmap, vmin=vmin, vmax = vmax) 
    plt.colorbar() 

    ax.set_xlim(region_info['lon_min'], region_info['lon_max'])
    ax.set_ylim(region_info['lat_min'], region_info['lat_max'])
    ax.set_xticks(range(region_info['lon_min'], region_info['lon_max'], 2))
    ax.set_xticklabels(region_info['lon_ticks'],fontsize=12)
    ax.set_yticks(range(region_info['lat_min'], region_info['lat_max'], 2))
    ax.set_yticklabels(region_info['lat_ticks'],fontsize=12)
    ax.set_xlabel('Longitude',fontsize=16)
    ax.set_ylabel('Latitude',fontsize=16)

    if savefig :
        plt.savefig(dir_output+'plot_'+method+'_'+region_info['name']+'_'+orbit+'_'+str(date_plot)+'_'+variable+'.png') 
        plt.show()
    else :     
        return fig
    
    
    

def create_movie(ds_passes, var_type = 'ssh',name_denoised='ssha_denoised_unet', method='SWOT', region_info=None, 
                 time0 = np.datetime64('2023-06-03'), ndays = 38, dir_output='../figures/', 
                 dim_name=['time','latitude','longitude'],orbit = '1d', framerate=3, Display=True, newmovie=False):
 


    # For memory leak when saving multiple png files...
    import matplotlib.pyplot as plt
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.gridspec as gridspec
    import gc
    import os
    import subprocess
    from IPython.display import Video


    moviename = 'movie_'+method+'_'+region_info['name']+'_'+var_type+'_'+orbit+'.mp4'

    if not os.path.isfile(os.path.join(dir_output, moviename)) or  newmovie: 

        # Create merged dataset 
        name_dim_rmse = (dim_name[0],)
        coords = (dim_name[0],dim_name[1],dim_name[2])




        # Plotting function
        def _save_single_frame(ids, tt, time0, name_denoised):

            if tt==0:
                return

            fig = plt.figure(figsize=(14,6))

            gs = gridspec.GridSpec(1,1,width_ratios=(1,))


            ax1 = fig.add_subplot(gs[0, 0])
            
            
            fig = plot_snapshots(ids, variable = var_type, name_denoised=name_denoised, 
                                 date_plot=time0, region_info=region_info, method=method, savefig=False, fig=fig)
 
 

            fig.savefig(dir_output+'/frame_'+method+'_'+region_info['name']+'_'+var_type+'_'+str(tt).zfill(5)+'.png',dpi=100)



            plt.close(fig)
            del fig
            gc.collect(2)


        # Compute and save frames 

        tt = 0 
        for idays in range(ndays):

            time1 = time0 + np.timedelta64(1,'D')

            print(time0)

            ds_passes_now=ds_passes.where(ds_passes.time>time0,drop=True)
            ds_passes_now=ds_passes_now.where(ds_passes_now.time<time1,drop=True) 
             

            if ds_passes_now.num_lines.size != 0: 
                _save_single_frame(ds_passes_now, tt, time0, name_denoised)

            time0 += np.timedelta64(1,'D')
            tt += 1

        # Create movie
        sourcefolder = dir_output
        frame_pattern = 'frame_'+method+'_'+region_info['name']+'_'+var_type+'_*.png'
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
    if Display :
        return Video(os.path.join(dir_output, moviename),embed=True)

    
    

def create_movie_old(ds_passes, name_var = 'ssha', method='SWOT', region='GS', lonlat_minmax=[309,314,30,42], 
                 time0 = np.datetime64('2023-06-03'), ndays = 38, dir_output='../figures/', 
                 dim_name=['time','latitude','longitude'], framerate=3, Display=True, clim=[-0.6,0.6], 
                 cmap='Spectral', newmovie=False):
 


    # For memory leak when saving multiple png files...
    import matplotlib.pyplot as plt
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




        # Plotting function
        def _save_single_frame(ids, tt, time0,clim=clim,cmap=cmap):

            if tt==0:
                return

            fig = plt.figure(figsize=(5,5))


            gs = gridspec.GridSpec(1,1,width_ratios=(1,))


            ax1 = fig.add_subplot(gs[0, 0])

            extend_style = 'both'
            if clim[0] == 0:
                extend_style = 'max'
            if clim[1] == 0:
                extend_style = 'min'


            plt.scatter(ids.longitude[:,:],ids.latitude[:,:],c=ids[name_var][:,:],cmap=cmap, vmin=clim[0],vmax=clim[1])  
            plt.axis(lonlat_minmax)
            plt.colorbar(extend=extend_style)

            ax1.set_title(time0)   

            fig.savefig(f'{dir_output}/frame_{method}_{region}_{name_var}_{str(tt).zfill(5)}.png',dpi=100)



            plt.close(fig)
            del fig
            gc.collect(2)


        # Compute and save frames 

        tt = 0 
        for idays in range(ndays):

            time1 = time0 + np.timedelta64(1,'D')

            print(time0)

            ds_passes_now=ds_passes.where(ds_passes.time>time0,drop=True)
            ds_passes_now=ds_passes_now.where(ds_passes_now.time<time1,drop=True) 
             

            if ds_passes_now.num_lines.size != 0: 
                _save_single_frame(ds_passes_now, tt, time0)

            time0 += np.timedelta64(1,'D')
            tt += 1

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
    if Display :
        return Video(os.path.join(dir_output, moviename),embed=True)
