
import numpy as np 
import xarray as xr 

def create_movie(ds_passes, name_var = 'ssha', method='SWOT', region='GS', lonlat_minmax=[309,314,30,42], 
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


            plt.scatter(ids.longitude[:,:34],ids.latitude[:,:34],c=ids[name_var][:,:34],cmap=cmap, vmin=clim[0],vmax=clim[1]) 
            plt.scatter(ids.longitude[:,35:],ids.latitude[:,35:],c=ids[name_var][:,35:],cmap=cmap, vmin=clim[0],vmax=clim[1])
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
