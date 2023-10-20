import logging

import cartopy.crs as ccrs
import hvplot.xarray
import numpy as np
import pyinterp
import scipy.signal
import xarray as xr
from netCDF4 import Dataset 
import sys
import operator   
import numpy.ma as ma
import dask.array as dsar
from dask import delayed
from functools import reduce
import scipy.signal as sps
import scipy.linalg as spl
import scipy.signal as signal
sys.path.append('..')
import src.mod_powerspec as ps 
import src.mod_xscale as xfft

#from src.mod_interp import *

import warnings
warnings.filterwarnings("ignore")


def compute_median_dx(dataset):
    """
    Compute the median distance between consecutive points in a dataset.

    Parameters
    ----------
    dataset : xarray.Dataset
        Dataset containing longitude and latitude data.

    Returns
    -------
    float
        The median distance in meters.
    """
        
    return 0.001*np.median(pyinterp.geodetic.coordinate_distances(dataset['longitude'][:-1].values,
                                                              dataset['latitude'][:-1].values,
                                                              dataset['longitude'][1:].values,
                                                              dataset['latitude'][1:].values
                                                             ))


def compute_median_lon_lat(vlon, vlat, sub_segment_point, npt):
    """
    Compute the median longitude and latitude for a sub-segment of points.

    Parameters
    ----------
    vlon : numpy.ndarray
        Array of longitudes.
    vlat : numpy.ndarray
        Array of latitudes.
    sub_segment_point : int
        Starting index of the sub-segment.
    npt : int
        Number of points in the sub-segment.

    Returns
    -------
    float
        The median longitude.
    float
        The median latitude.
    """
    
    # Near Greenwhich case
    if ((vlon[sub_segment_point + npt - 1] < 50.)
        and (vlon[sub_segment_point] > 320.)) \
            or ((vlon[sub_segment_point + npt - 1] > 320.)
                and (vlon[sub_segment_point] < 50.)):
        
        # make lon negative for lon > 180.
        tmp_lon = np.where(vlon[sub_segment_point:sub_segment_point + npt] > 180,
                           vlon[sub_segment_point:sub_segment_point + npt] - 360,
                           vlon[sub_segment_point:sub_segment_point + npt])
        # compute median longitude
        median_lon_sub_segment = np.median(tmp_lon)%360.
    
    # Far from greenwich
    else:
        median_lon_sub_segment = np.median(vlon[sub_segment_point:sub_segment_point + npt])

    median_lat_sub_segment = np.median(vlat[sub_segment_point:sub_segment_point + npt])
    
    return median_lon_sub_segment, median_lat_sub_segment


def compute_segment(ds_interp, npt, ref_var_name='sla_unfiltered', study_var_name='msla_interpolated', max_time_gap=np.timedelta64(2, 's'), segment_overlapping=0.25):
    """
    Compute segments of data based on specified criteria.

    Parameters
    ----------
    ds_interp : xarray.Dataset
        Interpolated dataset containing relevant variables.
    npt : int
        Number of points in each segment.
    ref_var_name : str, optional
        Name of the reference variable. Default is 'sla_unfiltered'.
    study_var_name : str, optional
        Name of the study variable. Default is 'msla_interpolated'.
    max_time_gap : numpy.timedelta64, optional
        Maximum time gap to consider for segmenting. Default is 2 seconds.
    segment_overlapping : float, optional
        Overlapping factor between segments. Default is 0.25.

    Returns
    -------
    numpy.ndarray
        Array of median longitudes for each segment.
    numpy.ndarray
        Array of median latitudes for each segment.
    numpy.ndarray
        Array of segmented reference variable data.
    numpy.ndarray
        Array of segmented study variable data.
    """
    
    #ds_interp = ds_interp.sortby('time').dropna('time')
    
    lon_along_track = ds_interp['longitude'].values
    lat_along_track = ds_interp['latitude'].values
    sla = ds_interp[ref_var_name].values
    msla = ds_interp[study_var_name].values
        
    indi = np.where((np.diff(ds_interp['time']) > max_time_gap))[0]
    if len(indi) > 0:
        track_segment_lenght = np.insert(np.diff(indi), [0], indi[0])
        selected_track_segment = np.where(track_segment_lenght >= npt)[0]
    else:
        selected_track_segment = np.asarray([])
    
    # indi = np.where((np.diff(ds_interp['time']) > max_time_gap))[0]
    # track_segment_lenght = np.insert(np.diff(indi), [0], indi[0])
    # selected_track_segment = np.where(track_segment_lenght >= npt)[0]

    list_lat_segment = []
    list_lon_segment = []
    list_sla_segment = []
    list_msla_segment = []

    if selected_track_segment.size > 0:

        for track in selected_track_segment:

            if track-1 >= 0:
                index_start_selected_track = indi[track-1]
                index_end_selected_track = indi[track]
            else:
                index_start_selected_track = 0
                index_end_selected_track = indi[track]

            start_point = index_start_selected_track
            end_point = index_end_selected_track

            for sub_segment_point in range(start_point, end_point - npt, int(npt*segment_overlapping)):
            
                mean_lon_sub_segment, mean_lat_sub_segment = compute_median_lon_lat(lon_along_track, 
                                                                                    lat_along_track,
                                                                                    sub_segment_point,
                                                                                    npt)

                sla_segment = sla[sub_segment_point:sub_segment_point + npt]
                msla_segment = msla[sub_segment_point:sub_segment_point + npt]
                
                list_sla_segment.append(sla_segment)
                list_lon_segment.append(mean_lon_sub_segment)
                list_lat_segment.append(mean_lat_sub_segment)
                list_msla_segment.append(msla_segment)
                
    return np.asarray(list_lon_segment), np.asarray(list_lat_segment), np.asarray(list_sla_segment), np.asarray(list_msla_segment)


def compute_segment_v2(ds_interp, npt, delta_x, ref_var_name='sla_unfiltered', study_var_name='msla_interpolated', max_time_gap=np.timedelta64(2, 's'), segment_overlapping=0.25):
    """
    Compute various spectral features for segments of data.

    Parameters
    ----------
    ds_interp : xarray.Dataset
        Interpolated dataset containing relevant variables.
    npt : int
        Number of points in each segment.
    delta_x : float
        Spatial resolution for spectral computation.
    ref_var_name : str, optional
        Name of the reference variable. Default is 'sla_unfiltered'.
    study_var_name : str, optional
        Name of the study variable. Default is 'msla_interpolated'.
    max_time_gap : numpy.timedelta64, optional
        Maximum time gap to consider for segmenting. Default is 2 seconds.
    segment_overlapping : float, optional
        Overlapping factor between segments. Default is 0.25.

    Returns
    -------
    numpy.ndarray
        Array of wavenumbers.
    numpy.ndarray
        Array of latitudes.
    numpy.ndarray
        Array of longitudes.
    numpy.ndarray
        Array of segment counts per grid cell.
    numpy.ndarray
        Array of reference variable power spectral densities.
    numpy.ndarray
        Array of study variable power spectral densities.
    numpy.ndarray
        Array of power spectral densities of the difference between study and reference.
    numpy.ndarray
        Array of coherence values.
    numpy.ndarray
        Array of cross-spectral densities.
    """
 
    #ds_interp = ds_interp.sortby('time').dropna('time')
    
    lon_along_track = ds_interp['longitude'].values
    lat_along_track = ds_interp['latitude'].values
    sla = ds_interp[ref_var_name].values
    msla = ds_interp[study_var_name].values
        
    indi = np.where((np.diff(ds_interp['time']) > max_time_gap))[0]
    if len(indi) > 0:
        track_segment_lenght = np.insert(np.diff(indi), [0], indi[0])
        selected_track_segment = np.where(track_segment_lenght >= npt)[0]
    else:
        selected_track_segment = np.asarray([])
    
    # indi = np.where((np.diff(ds_interp['time']) > max_time_gap))[0]
    # track_segment_lenght = np.insert(np.diff(indi), [0], indi[0])
    # selected_track_segment = np.where(track_segment_lenght >= npt)[0]

    fs = 1.0 / delta_x
    
    list_lat_segment = []
    list_lon_segment = []
    list_sla_segment = []
    list_msla_segment = []
    
    list_psd_ref = []
    list_nb_segment = []
    list_frequency = []
    list_psd_study = []
    list_psd_diff_study_ref = []
    list_coherence = []
    list_cross_spectrum = []

    if selected_track_segment.size > 0:

        for track in selected_track_segment:

            if track-1 >= 0:
                index_start_selected_track = indi[track-1]
                index_end_selected_track = indi[track]
            else:
                index_start_selected_track = 0
                index_end_selected_track = indi[track]

            start_point = index_start_selected_track
            end_point = index_end_selected_track

            for sub_segment_point in range(start_point, end_point - npt, int(npt*segment_overlapping)):
            
                mean_lon_sub_segment, mean_lat_sub_segment = compute_median_lon_lat(lon_along_track, 
                                                                                    lat_along_track,
                                                                                    sub_segment_point,
                                                                                    npt)

                sla_segment = sla[sub_segment_point:sub_segment_point + npt]
                msla_segment = msla[sub_segment_point:sub_segment_point + npt]
                
                wavenumber, psd_ref = scipy.signal.welch(sla_segment,
                                                         fs=fs,
                                                         nperseg=npt,
                                                         scaling='density',
                                                         noverlap=0)
                
                list_psd_ref.append(psd_ref)
                list_frequency.append(wavenumber)
        
                diff_study_ref = msla_segment - sla_segment

                # Power spectrum density of the error between to field
                _, psd_diff_study_ref = scipy.signal.welch(diff_study_ref,
                                                                        fs=fs,
                                                                        nperseg=npt,
                                                                        scaling='density',
                                                                        noverlap=0)

                # Power spectrum density study field
                _, psd_study = scipy.signal.welch(msla_segment,
                                                               fs=fs,
                                                               nperseg=npt,
                                                               scaling='density',
                                                               noverlap=0)

                # Magnitude square coherence between the ref and study field
                _, coherence = scipy.signal.coherence(msla_segment,
                                              sla_segment,
                                              fs=fs,
                                              nperseg=npt,
                                              noverlap=0)

                # Cross spectrum
                # _, cross_spectrum = scipy.signal.csd(msla_segment,
                #                                                   sla_segment,
                #                                                   fs=fs,
                #                                                   nperseg=npt,
                #                                                   noverlap=0)

                list_psd_study.append(psd_study)
                list_psd_diff_study_ref.append(psd_diff_study_ref)
                list_coherence.append(coherence)
                # list_cross_spectrum.append(cross_spectrum)
                
                
                
                list_sla_segment.append(sla_segment)
                list_lon_segment.append(mean_lon_sub_segment)
                list_lat_segment.append(mean_lat_sub_segment)
                list_msla_segment.append(msla_segment)
                
    return list_lon_segment, list_lat_segment, list_sla_segment, list_msla_segment


def spectral_computation(lon_segment, lat_segment, ref_segments, study_segments, delta_x, npt):

    """
    Perform spectral computation and analysis for grid cells.

    Parameters
    ----------
    lon_segment : numpy.ndarray
        Array of segment longitudes.
    lat_segment : numpy.ndarray
        Array of segment latitudes.
    ref_segments : numpy.ndarray
        Array of reference variable segments.
    study_segments : numpy.ndarray
        Array of study variable segments.
    delta_x : float
        Spatial resolution for spectral computation.
    npt : int
        Number of points in each segment.

    Returns
    -------
    numpy.ndarray
        Array of wavenumbers.
    numpy.ndarray
        Array of latitudes.
    numpy.ndarray
        Array of longitudes.
    numpy.ndarray
        Array of segment counts per grid cell.
    numpy.ndarray
        Array of reference variable power spectral densities.
    numpy.ndarray
        Array of study variable power spectral densities.
    numpy.ndarray
        Array of power spectral densities of the difference between study and reference.
    numpy.ndarray
        Array of coherence values.
    numpy.ndarray
        Array of cross-spectral densities.
    """ 


    delta_lat = 10.
    delta_lon = 10.
    nb_min_segment = 2
    vlon = np.arange(0., 360., 1.)
    vlat = np.arange(-80., 91., 1)

    list_mean_psd_ref = []
    list_nb_segment = []
    list_mean_frequency = []
    list_mean_psd_study = []
    list_mean_psd_diff_study_ref = []
    list_mean_coherence = []
    list_mean_cross_spectrum = []

    fs = 1.0 / delta_x

    # Loop over output lon/lat boxes and selection of the segment within the box plus/minus delta_lon/lat
    for ilat in vlat:

        lat_min = ilat - 0.5*delta_lat
        lat_max = ilat + 0.5*delta_lat

        selected_lat_index = np.where(np.logical_and(lat_segment >= lat_min, lat_segment <= lat_max))[0]
        ref_segments_tmp = ref_segments[selected_lat_index]
        if study_segments is not None:
            study_segments_tmp = study_segments[selected_lat_index]
        else:
            study_segments_tmp = None

        for ilon in vlon:

            lon_min = ilon - 0.5*delta_lon
            lon_max = ilon + 0.5*delta_lon

            if (lon_min < 0.) and (lon_max > 0.):
                selected_segment = np.where(np.logical_or(lon_segment[selected_lat_index] % 360. >= lon_min + 360.,
                                                          lon_segment[selected_lat_index] % 360. <= lon_max))[0]
            elif (lon_min > 0.) and (lon_max > 360.):
                selected_segment = np.where(np.logical_or(lon_segment[selected_lat_index] % 360. >= lon_min,
                                                          lon_segment[selected_lat_index] % 360. <= lon_max - 360.))[0]
            else:
                selected_segment = np.where(np.logical_and(lon_segment[selected_lat_index] % 360. >= lon_min,
                                                           lon_segment[selected_lat_index] % 360. <= lon_max))[0]

            if len(selected_segment) > nb_min_segment:
                selected_ref_segments = np.ma.masked_where(ref_segments_tmp[selected_segment].flatten() > 1.E10,
                                                           ref_segments_tmp[selected_segment].flatten())

                # Power spectrum density reference field
                wavenumber, psd_ref = scipy.signal.welch(selected_ref_segments,
                                                         fs=fs,
                                                         nperseg=npt,
                                                         scaling='density',
                                                         noverlap=0)
                
                wavenumber_to_keep = wavenumber

                list_mean_frequency.append(wavenumber)
                list_mean_psd_ref.append(psd_ref)
                list_nb_segment.append(selected_segment.size)


                selected_study_segments = np.ma.masked_where(
                        study_segments_tmp[selected_segment].flatten() > 1.E10,
                        study_segments_tmp[selected_segment].flatten())

                # Compute diff study minus ref
                diff_study_ref = selected_study_segments - selected_ref_segments

                # Power spectrum density of the error between to field
                wavenumber, psd_diff_study_ref = scipy.signal.welch(diff_study_ref,
                                                                        fs=fs,
                                                                        nperseg=npt,
                                                                        scaling='density',
                                                                        noverlap=0)

                # Power spectrum density study field
                wavenumber, psd_study = scipy.signal.welch(selected_study_segments,
                                                               fs=fs,
                                                               nperseg=npt,
                                                               scaling='density',
                                                               noverlap=0)

                # Magnitude square coherence between the ref and study field
                wavenumber, coherence = scipy.signal.coherence(selected_study_segments,
                                                                   selected_ref_segments,
                                                                   fs=fs,
                                                                   nperseg=npt,
                                                                   noverlap=0)

                # Cross spectrum
                wavenumber, cross_spectrum = scipy.signal.csd(selected_study_segments,
                                                                  selected_ref_segments,
                                                                  fs=fs,
                                                                  nperseg=npt,
                                                                  noverlap=0)

                list_mean_psd_study.append(psd_study)
                list_mean_psd_diff_study_ref.append(psd_diff_study_ref)
                list_mean_coherence.append(coherence)
                list_mean_cross_spectrum.append(cross_spectrum)
                

            else:

                list_mean_frequency.append(np.zeros(npt))
                list_mean_psd_ref.append(np.zeros(npt))
                list_nb_segment.append(0.)
                list_mean_psd_study.append(np.zeros(npt))
                list_mean_psd_diff_study_ref.append(np.zeros(npt))
                list_mean_coherence.append(np.zeros(npt))
                list_mean_cross_spectrum.append(np.zeros(npt))
                
#                 list_mean_frequency.append(np.zeros((int(npt / 2) + 1)))
#                 list_mean_psd_ref.append(np.zeros((int(npt / 2) + 1)))
#                 list_nb_segment.append(0.)

#                 list_mean_psd_study.append(np.zeros((int(npt / 2) + 1)))
#                 list_mean_psd_diff_study_ref.append(np.zeros((int(npt / 2) + 1)))
#                 list_mean_coherence.append(np.zeros((int(npt / 2) + 1)))
#                 list_mean_cross_spectrum.append(np.zeros((int(npt / 2) + 1)))

    # wavenumber = np.asarray(list_mean_frequency)
    # wavenumber = np.ma.masked_where(wavenumber == 0., wavenumber)
    # try:
    #     wavenumber = np.ma.masked_invalid(wavenumber)
    #     wavenumber =  np.ma.mean(wavenumber, axis=0).filled(0.)
    # except:
    #     wavenumber =  np.ma.mean(wavenumber, axis=0)
        
    # wavenumber = np.ma.mean(np.ma.masked_invalid(np.ma.masked_where(np.asarray(list_mean_frequency) == 0,
    #                                                           np.asarray(list_mean_frequency))), axis=0).filled(0.)
    
    nb_segment = np.asarray(list_nb_segment).reshape((vlat.size, vlon.size))
    psd_ref = np.transpose(np.asarray(list_mean_psd_ref)).reshape((wavenumber_to_keep.size, vlat.size, vlon.size))
    psd_study = np.transpose(np.asarray(list_mean_psd_study)).reshape((wavenumber_to_keep.size, vlat.size, vlon.size))
    psd_diff = np.transpose(np.asarray(list_mean_psd_diff_study_ref)).reshape((wavenumber_to_keep.size, vlat.size, vlon.size))
    coherence = np.transpose(np.asarray(list_mean_coherence)).reshape((wavenumber_to_keep.size, vlat.size, vlon.size))
    cross_spectrum = np.transpose(np.asarray(list_mean_cross_spectrum)).reshape((wavenumber_to_keep.size, vlat.size, vlon.size))
    
    return wavenumber_to_keep, vlat, vlon, nb_segment, psd_ref, psd_study, psd_diff, coherence, cross_spectrum



def spectral_computation_v2(lon_segment, lat_segment, ref_segments, study_segments, delta_x, npt):
    """
    Perform spectral computation and analysis for grid cells.

    Parameters
    ----------
    lon_segment : numpy.ndarray
        Array of segment longitudes.
    lat_segment : numpy.ndarray
        Array of segment latitudes.
    ref_segments : numpy.ndarray
        Array of reference variable segments.
    study_segments : numpy.ndarray
        Array of study variable segments.
    delta_x : float
        Spatial resolution for spectral computation.
    npt : int
        Number of points in each segment.

    Returns
    -------
    numpy.ndarray
        Array of wavenumbers.
    numpy.ndarray
        Array of latitudes.
    numpy.ndarray
        Array of longitudes.
    numpy.ndarray
        Array of segment counts per grid cell.
    numpy.ndarray
        Array of reference variable power spectral densities.
    numpy.ndarray
        Array of study variable power spectral densities.
    numpy.ndarray
        Array of power spectral densities of the difference between study and reference.
    numpy.ndarray
        Array of coherence values.
    numpy.ndarray
        Array of cross-spectral densities.
    """

    delta_lat = 10.
    delta_lon = 10.
    nb_min_segment = 2
    vlon = np.arange(0., 360., 1.)
    vlat = np.arange(-80., 91., 1)

    list_mean_psd_ref = []
    list_nb_segment = []
    list_mean_frequency = []
    list_mean_psd_study = []
    list_mean_psd_diff_study_ref = []
    list_mean_coherence = []
    list_mean_cross_spectrum = []

    fs = 1.0 / delta_x

    # Loop over output lon/lat boxes and selection of the segment within the box plus/minus delta_lon/lat
    for ilat in vlat:

        lat_min = ilat - 0.5*delta_lat
        lat_max = ilat + 0.5*delta_lat

        selected_lat_index = np.where(np.logical_and(lat_segment >= lat_min, lat_segment <= lat_max))[0]
        ref_segments_tmp = ref_segments[selected_lat_index]
        if study_segments is not None:
            study_segments_tmp = study_segments[selected_lat_index]
        else:
            study_segments_tmp = None

        for ilon in vlon:

            lon_min = ilon - 0.5*delta_lon
            lon_max = ilon + 0.5*delta_lon

            if (lon_min < 0.) and (lon_max > 0.):
                selected_segment = np.where(np.logical_or(lon_segment[selected_lat_index] % 360. >= lon_min + 360.,
                                                          lon_segment[selected_lat_index] % 360. <= lon_max))[0]
            elif (lon_min > 0.) and (lon_max > 360.):
                selected_segment = np.where(np.logical_or(lon_segment[selected_lat_index] % 360. >= lon_min,
                                                          lon_segment[selected_lat_index] % 360. <= lon_max - 360.))[0]
            else:
                selected_segment = np.where(np.logical_and(lon_segment[selected_lat_index] % 360. >= lon_min,
                                                           lon_segment[selected_lat_index] % 360. <= lon_max))[0]

            if len(selected_segment) > nb_min_segment:
                selected_ref_segments = np.ma.masked_where(ref_segments_tmp[selected_segment].flatten() > 1.E10,
                                                           ref_segments_tmp[selected_segment].flatten())

                # Power spectrum density reference field
                wavenumber, psd_ref = scipy.signal.welch(selected_ref_segments,
                                                         fs=fs,
                                                         nperseg=npt,
                                                         scaling='density',
                                                         noverlap=0)
                
                wavenumber_to_keep = wavenumber

                list_mean_frequency.append(wavenumber)
                list_mean_psd_ref.append(psd_ref)
                list_nb_segment.append(selected_segment.size)


                selected_study_segments = np.ma.masked_where(
                        study_segments_tmp[selected_segment].flatten() > 1.E10,
                        study_segments_tmp[selected_segment].flatten())

                # Compute diff study minus ref
                diff_study_ref = selected_study_segments - selected_ref_segments

                # Power spectrum density of the error between to field
                wavenumber, psd_diff_study_ref = scipy.signal.welch(diff_study_ref,
                                                                        fs=fs,
                                                                        nperseg=npt,
                                                                        scaling='density',
                                                                        noverlap=0)

                # Power spectrum density study field
                wavenumber, psd_study = scipy.signal.welch(selected_study_segments,
                                                               fs=fs,
                                                               nperseg=npt,
                                                               scaling='density',
                                                               noverlap=0)

                # Magnitude square coherence between the ref and study field
                wavenumber, coherence = scipy.signal.coherence(selected_study_segments,
                                                                   selected_ref_segments,
                                                                   fs=fs,
                                                                   nperseg=npt,
                                                                   noverlap=0)

                # Cross spectrum
                wavenumber, cross_spectrum = scipy.signal.csd(selected_study_segments,
                                                                  selected_ref_segments,
                                                                  fs=fs,
                                                                  nperseg=npt,
                                                                  noverlap=0)

                list_mean_psd_study.append(psd_study)
                list_mean_psd_diff_study_ref.append(psd_diff_study_ref)
                list_mean_coherence.append(coherence)
                list_mean_cross_spectrum.append(cross_spectrum)
                

            else:

                # list_mean_frequency.append(np.zeros(npt))
                # list_mean_psd_ref.append(np.zeros(npt))
                # list_nb_segment.append(0.)
                # list_mean_psd_study.append(np.zeros(npt))
                # list_mean_psd_diff_study_ref.append(np.zeros(npt))
                # list_mean_coherence.append(np.zeros(npt))
                # list_mean_cross_spectrum.append(np.zeros(npt))
                
                list_mean_frequency.append(np.zeros((int(npt / 2) + 1)))
                list_mean_psd_ref.append(np.zeros((int(npt / 2) + 1)))
                list_nb_segment.append(0.)
                list_mean_psd_study.append(np.zeros((int(npt / 2) + 1)))
                list_mean_psd_diff_study_ref.append(np.zeros((int(npt / 2) + 1)))
                list_mean_coherence.append(np.zeros((int(npt / 2) + 1)))
                list_mean_cross_spectrum.append(np.zeros((int(npt / 2) + 1)))

    # wavenumber = np.asarray(list_mean_frequency)
    # wavenumber = np.ma.masked_where(wavenumber == 0., wavenumber)
    # try:
    #     wavenumber = np.ma.masked_invalid(wavenumber)
    #     wavenumber =  np.ma.mean(wavenumber, axis=0).filled(0.)
    # except:
    #     wavenumber =  np.ma.mean(wavenumber, axis=0)
        
    # wavenumber = np.ma.mean(np.ma.masked_invalid(np.ma.masked_where(np.asarray(list_mean_frequency) == 0,
    #                                                           np.asarray(list_mean_frequency))), axis=0).filled(0.)
    
    nb_segment = np.asarray(list_nb_segment).reshape((vlat.size, vlon.size))
    psd_ref = np.transpose(np.asarray(list_mean_psd_ref)).reshape((wavenumber_to_keep.size, vlat.size, vlon.size))
    psd_study = np.transpose(np.asarray(list_mean_psd_study)).reshape((wavenumber_to_keep.size, vlat.size, vlon.size))
    psd_diff = np.transpose(np.asarray(list_mean_psd_diff_study_ref)).reshape((wavenumber_to_keep.size, vlat.size, vlon.size))
    coherence = np.transpose(np.asarray(list_mean_coherence)).reshape((wavenumber_to_keep.size, vlat.size, vlon.size))
    cross_spectrum = np.transpose(np.asarray(list_mean_cross_spectrum)).reshape((wavenumber_to_keep.size, vlat.size, vlon.size))
    
    return wavenumber_to_keep, vlat, vlon, nb_segment, psd_ref, psd_study, psd_diff, coherence, cross_spectrum


def compute_crossing(array, wavenumber, threshold=0.5):
    """
    Compute the crossing of a threshold by an array of values.

    Parameters
    ----------
    array : numpy.ndarray
        Array of values.
    wavenumber : numpy.ndarray
        Array of wavenumbers.
    threshold : float, optional
        Threshold for crossing, defaults to 0.5.

    Returns
    -------
    float
        Effective resolution.
    bool
        Flag indicating multiple crossings.
    """


    flag_multiple_crossing = False
    zero_crossings = np.where(np.diff(np.sign(array - threshold)))[0]
    
    if len(zero_crossings) > 1:
        #print('Multiple crossing', len(zero_crossings))
        flag_multiple_crossing = True
        # MB add for large scale bais
        zero_crossings[0] = zero_crossings[-1]
         
    if len(zero_crossings) > 0:
        if zero_crossings[0] + 1 < array.size:

            array1 = array[zero_crossings[0]] - threshold
            array2 = array[zero_crossings[0] + 1] - threshold
            dist1 = np.log(wavenumber[zero_crossings[0]])
            dist2 = np.log(wavenumber[zero_crossings[0] + 1])
            log_wavenumber_crossing = dist1 - array1 * (dist1 - dist2) / (array1 - array2)
            resolution_scale = 1. / np.exp(log_wavenumber_crossing)

        else:
            resolution_scale = 0.

    else:
        resolution_scale = 0.

    return resolution_scale, flag_multiple_crossing


def compute_resolution(lon, lat, wavenumber, psd_diff, psd_ref):
    """
    Compute the resolution.

    Parameters
    ----------
    lon : numpy.ndarray
        Array of longitudes.
    lat : numpy.ndarray
        Array of latitudes.
    wavenumber : numpy.ndarray
        Array of wavenumbers.
    psd_diff : numpy.ndarray
        Power spectral density of the difference between study and reference.
    psd_ref : numpy.ndarray
        Power spectral density of the reference.

    Returns
    -------
    numpy.ndarray
        Array of resolution values.
    """
    
    ratio = psd_diff/psd_ref

    resolution = np.empty((lat.size, lon.size))
    for jj in range(lat.size):
        for ii in range(lon.size):
            if not np.ma.is_masked(ratio[:, jj, ii]):
                resolution[jj, ii], flag = compute_crossing(ratio[:, jj, ii], wavenumber)
    
    return resolution


def write_psd_output(output_netcdf_file, wavenumber, vlat, vlon, nb_segment, psd_ref, psd_study, psd_diff_ref_study, coherence, cross_spectrum, one_sided=True, method_name=' '): 
    """
    Write power spectral density and related variables to a NetCDF file.

    Parameters
    ----------
    output_netcdf_file : str
        Output NetCDF file path.
    wavenumber : numpy.ndarray
        Array of wavenumbers.
    vlat : numpy.ndarray
        Array of latitudes.
    vlon : numpy.ndarray
        Array of longitudes.
    nb_segment : numpy.ndarray
        Array of segment counts per grid cell.
    psd_ref : numpy.ndarray
        Power spectral density of the reference field.
    psd_study : numpy.ndarray
        Power spectral density of the study field.
    psd_diff_ref_study : numpy.ndarray
        Power spectral density of the difference between study and reference fields.
    coherence : numpy.ndarray
        Magnitude squared coherence between reference and study fields.
    cross_spectrum : numpy.ndarray
        Complex cross-spectrum between reference and study fields.
    one_sided : bool, optional
        Whether to use one-sided wavenumbers, defaults to True.
    method_name : str, optional
        Method name for the NetCDF file, defaults to ' '.
    """
    
    nc_out = Dataset(output_netcdf_file, 'w', format='NETCDF4')
    
    
    if one_sided:
        positive_index = np.where(wavenumber >=0)[0]
    else:
        positive_index = np.arange(wavenumber.size)
    
    nc_out.createDimension('wavenumber', wavenumber[positive_index].size)
    nc_out.createDimension('lat', vlat.size)
    nc_out.createDimension('lon', vlon.size)

    wavenumber_out = nc_out.createVariable('wavenumber', 'f8', 'wavenumber', zlib=True)
    wavenumber_out.units = "1/km"
    wavenumber_out.axis = 'T'
    wavenumber_out[:] = wavenumber[positive_index]

    nb_segment_out = nc_out.createVariable('nb_segment', 'f8', ('lat', 'lon'), zlib=True)
    nb_segment_out.long_name = "number of segment used in spectral computation"
    nb_segment_out[:, :] = np.ma.masked_where(nb_segment == 0., nb_segment).filled(np.nan)

    lat_out = nc_out.createVariable('lat', 'f4', 'lat', zlib=True)
    lat_out[:] = vlat
    lat_out.units = 'degree_north'
    lon_out = nc_out.createVariable('lon', 'f4', 'lon', zlib=True)
    lon_out[:] = vlon
    lon_out.units = 'degree_east'

    psd_ref_out = nc_out.createVariable('psd_ref', 'f8', ('wavenumber', 'lat', 'lon'), zlib=True)
    psd_ref_out.units = 'm2/km'
    psd_ref_out.coordinates = "wavenumber lat lon"
    psd_ref_out.long_name = "power spectrum density reference field"
    psd_ref_out[:, :, :] = np.ma.masked_invalid(np.ma.masked_where(psd_ref == 0, psd_ref)).filled(np.nan)[positive_index, :, :]
    
    psd_study_out = nc_out.createVariable('psd_study', 'f8', ('wavenumber', 'lat', 'lon'), zlib=True)
    psd_study_out.units = 'm2/km'
    psd_study_out.coordinates = "wavenumber lat lon"
    psd_study_out.long_name = "power spectrum density study field"
    psd_study_out[:, :, :] = np.ma.masked_invalid(np.ma.masked_where(psd_study == 0, psd_study)).filled(np.nan)[positive_index, :, :]
    
    psd_diff = nc_out.createVariable('psd_diff', 'f8', ('wavenumber', 'lat', 'lon'), zlib=True)
    psd_diff.units = 'm2/km'
    psd_diff.coordinates = "wavenumber lat lon"
    psd_diff.long_name = "power spectrum density of difference study minus reference field"
    psd_diff[:, :, :] = np.ma.masked_invalid(np.ma.masked_where(psd_diff_ref_study == 0, psd_diff_ref_study))[positive_index, :, :]
    
    resolution = compute_resolution(vlon, vlat, wavenumber[positive_index], psd_diff_ref_study[positive_index, :, :], psd_ref[positive_index, :, :])
    resolution_out = nc_out.createVariable('effective_resolution', 'f8', ('lat', 'lon'), zlib=True)
    resolution_out.units = 'km'
    resolution_out.coordinates = "lat lon"
    resolution_out.long_name = "Effective resolution computed as wavelenght where psd_err/psd_ref=0.5"
    resolution_out[:, :] = np.ma.masked_invalid(np.ma.masked_where(resolution == 0, resolution)).filled(np.nan)

    coherence_out = nc_out.createVariable('coherence', 'f8', ('wavenumber', 'lat', 'lon'), zlib=True)
    coherence_out[:, :, :] = np.ma.masked_invalid(np.ma.masked_where(coherence == 0, coherence)).filled(np.nan)[positive_index, :, :]
    coherence_out.coordinates = "wavenumber lat lon"
    coherence_out.long_name = "magnitude squared coherence between reference and study fields"

    cross_spectrum_real_out = nc_out.createVariable('cross_spectrum_real', 'f8', ('wavenumber', 'lat', 'lon'),
                                                        zlib=True)
    cross_spectrum_real_out[:, :, :] = np.ma.masked_invalid(np.ma.masked_where(np.real(cross_spectrum) == 0., np.real(cross_spectrum))).filled(np.nan)[positive_index, :, :]
    cross_spectrum_real_out.coordinates = "wavenumber lat lon"
    cross_spectrum_real_out.long_name = "real part of cross_spectrum between reference and study fields"
    cross_spectrum_imag_out = nc_out.createVariable('cross_spectrum_imag', 'f8', ('wavenumber', 'lat', 'lon'),
                                                        zlib=True)
    cross_spectrum_imag_out[:, :, :] = np.ma.masked_invalid(np.ma.masked_where(np.imag(cross_spectrum) == 0., np.imag(cross_spectrum))).filled(np.nan)[positive_index, :, :]
    cross_spectrum_imag_out.coordinates = "wavenumber lat lon"
    cross_spectrum_imag_out.long_name = "imaginary part of cross_spectrum between reference and study fields"
 
    nc_out.method=method_name
    
    nc_out.close()
    
    
def compute_psd_scores(ds_interp, output_filename, lenght_scale=1500. ):
    """
    Compute power spectral density (PSD) scores and save them to a NetCDF file.

    Parameters
    ----------
    ds_interp : xarray.Dataset
        Interpolated dataset.
    output_filename : str
        Output NetCDF file path.
    length_scale : float, optional
        Length scale, defaults to 1500.0.

    Returns
    -------
    None
    """    
    # logging.info("Interpolate SLA maps onto alongtrack")
    # ds_interp = run_interpolation(ds_maps, ds_alongtrack)
    # ds_interp = ds_interp.dropna('time')
    
    logging.info('Segment computation...')
    delta_x = compute_median_dx(ds_interp) # in km
    npt = int(lenght_scale / delta_x)
    lon_segment, lat_segment, sla_segment, msla_segment = compute_segment(ds_interp, npt)
    
    logging.info('Spectral analysis...')
    wavenumber, vlat, vlon, nb_segment, psd_ref, psd_study, psd_diff, coherence, cross_spectrum = spectral_computation(lon_segment,
                                                                                                                       lat_segment,
                                                                                                                       sla_segment + 1j*np.zeros(np.shape(sla_segment)),
                                                                                                                       msla_segment + 1j*np.zeros(np.shape(msla_segment)), 
                                                                                                                       delta_x,
                                                                                                                       npt
                                                                                                                      )
    logging.info('Saving ouput...')
    write_psd_output(output_filename, wavenumber, vlat, vlon, nb_segment, psd_ref, psd_study, psd_diff, coherence, cross_spectrum)
    logging.info("PSD file saved as: %s", output_filename)
    
    

def compute_psd_scores_v2(ds_interp, output_filename, lenght_scale=1500., method_name=' '):
    """
    Compute power spectral density (PSD) scores using a different method and save them to a NetCDF file.

    Parameters
    ----------
    ds_interp : xarray.Dataset
        Interpolated dataset.
    output_filename : str
        Output NetCDF file path.
    length_scale : float, optional
        Length scale, defaults to 1500.0.
    method_name : str, optional
        Method name for the NetCDF file, defaults to ' '.

    Returns
    -------
    None
    """    
    # logging.info("Interpolate SLA maps onto alongtrack")
    # ds_interp = run_interpolation(ds_maps, ds_alongtrack)
    # ds_interp = ds_interp.dropna('time')
    
    logging.info('Segment computation...')
    delta_x = compute_median_dx(ds_interp) # in km
    npt = int(lenght_scale / delta_x)
    lon_segment, lat_segment, sla_segment, msla_segment = compute_segment(ds_interp, npt)
    
    logging.info('Spectral analysis...')
    wavenumber, vlat, vlon, nb_segment, psd_ref, psd_study, psd_diff, coherence, cross_spectrum = spectral_computation_v2(lon_segment,
                                                                                                                       lat_segment,
                                                                                                                       sla_segment,
                                                                                                                       msla_segment, 
                                                                                                                       delta_x,
                                                                                                                       npt
                                                                                                                      )
    logging.info('Saving ouput...') 
    write_psd_output(output_filename, wavenumber, vlat, vlon, nb_segment, psd_ref, psd_study, psd_diff, coherence, cross_spectrum, method_name=method_name)
    logging.info("PSD file saved as: %s", output_filename)
    
 


def spectral_computation_drifters_by_latbin(lon_segment, lat_segment, ref_segments, study_segments, delta_x, npt):
    """
    Compute spectral properties of drifters by latitude bin.

    Parameters
    ----------
    lon_segment : numpy.ndarray
        Array of segmented longitudes.
    lat_segment : numpy.ndarray
        Array of segmented latitudes.
    ref_segments : numpy.ndarray
        Array of segmented reference segments.
    study_segments : numpy.ndarray
        Array of segmented study segments.
    delta_x : float
        Delta x value.
    npt : int
        Number of points.

    Returns
    -------
    wavenumber_to_keep : numpy.ndarray
        Array of wavenumbers.
    vlat : numpy.ndarray
        Array of latitudes.
    nb_segment : numpy.ndarray
        Array of segment counts.
    psd_ref : numpy.ndarray
        Power spectral density of the reference.
    psd_study : numpy.ndarray
        Power spectral density of the study.
    psd_diff : numpy.ndarray
        Power spectral density of the difference.
    coherence : numpy.ndarray
        Magnitude squared coherence.
    cross_spectrum : numpy.ndarray
        Complex cross-spectrum.
    """

    delta_lat = 1.
    nb_min_segment = 2
    vlat = np.arange(-80., 91., 1)

    list_mean_psd_ref = []
    list_nb_segment = []
    list_mean_frequency = []
    list_mean_psd_study = []
    list_mean_psd_diff_study_ref = []
    list_mean_coherence = []
    list_mean_cross_spectrum = []

    fs = 1.0 / delta_x

    # Loop over output lon/lat boxes and selection of the segment within the box plus/minus delta_lon/lat
    for ilat in vlat:

        lat_min = ilat - 0.5*delta_lat
        lat_max = ilat + 0.5*delta_lat

        selected_lat_index = np.where(np.logical_and(lat_segment >= lat_min, lat_segment <= lat_max))[0]
        
        if len(selected_lat_index) > nb_min_segment:
            ref_segments_tmp = ref_segments[selected_lat_index]
            if study_segments is not None:
                study_segments_tmp = study_segments[selected_lat_index]
            else:
                study_segments_tmp = None
        
            selected_ref_segments = np.ma.masked_where(ref_segments_tmp.flatten() > 1.E10,
                                                           ref_segments_tmp.flatten())

            # Power spectrum density reference field
            wavenumber, psd_ref = scipy.signal.welch(selected_ref_segments,
                                                         fs=fs,
                                                         nperseg=npt,
                                                         scaling='density',
                                                         noverlap=0)
                
            wavenumber_to_keep = wavenumber

            list_mean_frequency.append(wavenumber)
            list_mean_psd_ref.append(psd_ref)
            list_nb_segment.append(selected_ref_segments.size)


            selected_study_segments = np.ma.masked_where(
                    study_segments_tmp.flatten() > 1.E10,
                    study_segments_tmp.flatten())

            # Compute diff study minus ref
            diff_study_ref = selected_study_segments - selected_ref_segments

            # Power spectrum density of the error between to field
            wavenumber, psd_diff_study_ref = scipy.signal.welch(diff_study_ref,
                                                                        fs=fs,
                                                                        nperseg=npt,
                                                                        scaling='density',
                                                                        noverlap=0)

            # Power spectrum density study field
            wavenumber, psd_study = scipy.signal.welch(selected_study_segments,
                                                               fs=fs,
                                                               nperseg=npt,
                                                               scaling='density',
                                                               noverlap=0)

            # Magnitude square coherence between the ref and study field
            wavenumber, coherence = scipy.signal.coherence(selected_study_segments,
                                                                   selected_ref_segments,
                                                                   fs=fs,
                                                                   nperseg=npt,
                                                                   noverlap=0)

            # Cross spectrum
            wavenumber, cross_spectrum = scipy.signal.csd(selected_study_segments,
                                                                  selected_ref_segments,
                                                                  fs=fs,
                                                                  nperseg=npt,
                                                                  noverlap=0)

            list_mean_psd_study.append(psd_study)
            list_mean_psd_diff_study_ref.append(psd_diff_study_ref)
            list_mean_coherence.append(coherence)
            list_mean_cross_spectrum.append(cross_spectrum)
                

        else:

            list_mean_frequency.append(np.zeros(npt))
            list_mean_psd_ref.append(np.zeros(npt))
            list_nb_segment.append(0.)
            list_mean_psd_study.append(np.zeros(npt))
            list_mean_psd_diff_study_ref.append(np.zeros(npt))
            list_mean_coherence.append(np.zeros(npt))
            list_mean_cross_spectrum.append(np.zeros(npt))
                
    nb_segment = np.asarray(list_nb_segment).reshape((vlat.size))
    psd_ref = np.transpose(np.asarray(list_mean_psd_ref)).reshape((wavenumber_to_keep.size, vlat.size))
    psd_study = np.transpose(np.asarray(list_mean_psd_study)).reshape((wavenumber_to_keep.size, vlat.size))
    psd_diff = np.transpose(np.asarray(list_mean_psd_diff_study_ref)).reshape((wavenumber_to_keep.size, vlat.size))
    coherence = np.transpose(np.asarray(list_mean_coherence)).reshape((wavenumber_to_keep.size, vlat.size))
    cross_spectrum = np.transpose(np.asarray(list_mean_cross_spectrum)).reshape((wavenumber_to_keep.size, vlat.size))
    
    return wavenumber_to_keep, vlat, nb_segment, psd_ref, psd_study, psd_diff, coherence, cross_spectrum



def compute_psd_scores_current(ds_interp, output_filename, lenght_scale=np.timedelta64(40, 'D'), method_name=' '):
    """
    Compute power spectral density (PSD) scores for current data and save them to a NetCDF file.

    Parameters
    ----------
    ds_interp : xarray.Dataset
        Interpolated dataset.
    output_filename : str
        Output NetCDF file path.
    length_scale : np.timedelta64, optional
        Length scale, defaults to 40 days.
    method_name : str, optional
        Method name for the NetCDF file, defaults to ' '.

    Returns
    -------
    None
    """
    
    logging.info('Segment computation...')
    delta_t = np.median(np.diff(np.unique(ds_interp['time'])))
    delta_t_days = (delta_t / np.timedelta64(1, 'D')).astype(np.float64)
    npt = int(lenght_scale / delta_t)
    
    ds_interp['complex_current_drifter'] = (("time"), ds_interp['EWCT'].data + 1j* ds_interp['NSCT'].data)
    ds_interp['complex_current_maps'] = (("time"), ds_interp['ugos_interpolated'].data + 1j* ds_interp['vgos_interpolated'].data)
    
    list_of_lon = []
    list_of_lat = []
    list_of_udrifter = []
    list_of_umaps = []
    
    for sensor_id in np.unique(ds_interp['sensor_id']):
        ds_sel = ds_interp.where(ds_interp['sensor_id'] == sensor_id, drop=True)
        lon_segment, lat_segment, udfriter_segment, umap_segment = compute_segment(ds_sel, npt, 
                                                                              ref_var_name='complex_current_drifter', 
                                                                              study_var_name='complex_current_maps',
                                                                              max_time_gap=np.timedelta64(10, 'h'),
                                                                              segment_overlapping=0.5)
        del ds_sel
        if (udfriter_segment.size > 0) and (umap_segment.size > 0):
            list_of_lon.append(lon_segment)
            list_of_lat.append(lat_segment)
            list_of_udrifter.append(udfriter_segment)
            list_of_umaps.append(umap_segment)
        del lon_segment, lat_segment, udfriter_segment, umap_segment
            
    lon_segment = np.concatenate(list_of_lon)
    lat_segment = np.concatenate(list_of_lat)
    udfriter_segment = np.concatenate(list_of_udrifter)
    umap_segment = np.concatenate(list_of_umaps)
    
    logging.info('Spectral analysis...')
    wavenumber, vlat, nb_segment, psd_ref, psd_study, psd_diff, coherence, cross_spectrum = spectral_computation_drifters_by_latbin(lon_segment,
                                                                                                                       lat_segment,
                                                                                                                       udfriter_segment,
                                                                                                                       umap_segment, 
                                                                                                                       delta_t_days,
                                                                                                                       npt
                                                                                                                      )
    
    logging.info('Write output...')
    write_psd_current_output(output_filename, wavenumber, vlat, nb_segment, psd_ref, psd_study, psd_diff, coherence, cross_spectrum,method_name=method_name)
    logging.info("PSD file saved as: %s", output_filename)


def write_psd_current_output(output_netcdf_file, wavenumber, vlat, nb_segment, psd_ref, psd_study, psd_diff_ref_study, coherence, cross_spectrum,method_name=''):
    """
    Write power spectral density (PSD) current data to a NetCDF file.

    Parameters
    ----------
    output_netcdf_file : str
        Output NetCDF file path.
    wavenumber : numpy.ndarray
        Array of wavenumbers.
    vlat : numpy.ndarray
        Array of latitudes.
    nb_segment : numpy.ndarray
        Array of segment counts.
    psd_ref : numpy.ndarray
        Power spectral density of the reference.
    psd_study : numpy.ndarray
        Power spectral density of the study.
    psd_diff_ref_study : numpy.ndarray
        Power spectral density of the difference between study and reference.
    coherence : numpy.ndarray
        Magnitude squared coherence.
    cross_spectrum : numpy.ndarray
        Complex cross-spectrum.
    method_name : str, optional
        Method name for the NetCDF file, defaults to an empty string.

    Returns
    -------
    None
    """
 
    nc_out = Dataset(output_netcdf_file, 'w', format='NETCDF4')
    
    # Reformate data
    sorted_index = np.argsort(wavenumber)
    wavenumber = wavenumber[sorted_index]
    for jj in range(psd_ref[0, :].size):
        psd_ref[:, jj] = psd_ref[sorted_index, jj]
        psd_study[:, jj] = psd_study[sorted_index, jj]
        psd_diff_ref_study[:, jj] = psd_diff_ref_study[sorted_index, jj]
        coherence[:, jj] = coherence[sorted_index, jj]
        cross_spectrum[:, jj] = cross_spectrum[sorted_index, jj]

    nc_out.createDimension('wavenumber', wavenumber.size)
    nc_out.createDimension('lat', vlat.size)

    wavenumber_out = nc_out.createVariable('wavenumber', 'f8', 'wavenumber', zlib=True)
    wavenumber_out.units = "1/day"
    wavenumber_out.axis = 'T'
    wavenumber_out[:] = wavenumber

    nb_segment_out = nc_out.createVariable('nb_segment', 'f8', ('lat'), zlib=True)
    nb_segment_out.long_name = "number of segment used in spectral computation"
    nb_segment_out[:] = np.ma.masked_where(nb_segment == 0., nb_segment).filled(np.nan)

    lat_out = nc_out.createVariable('lat', 'f4', 'lat', zlib=True)
    lat_out[:] = vlat
    lat_out.units = 'degree_north'

    psd_ref_out = nc_out.createVariable('psd_ref', 'f8', ('lat', 'wavenumber'), zlib=True)
    psd_ref_out.units = 'm2/km'
    psd_ref_out.coordinates = "wavenumber lat lon"
    psd_ref_out.long_name = "power spectrum density reference field"
    psd_ref_out[:, :] = np.ma.masked_invalid(np.ma.masked_where(psd_ref.T == 0, psd_ref.T)).filled(np.nan)[:, :]
    
    psd_study_out = nc_out.createVariable('psd_study', 'f8', ('lat', 'wavenumber'), zlib=True)
    psd_study_out.units = 'm2/km'
    psd_study_out.coordinates = "wavenumber lat lon"
    psd_study_out.long_name = "power spectrum density study field"
    psd_study_out[:, :] = np.ma.masked_invalid(np.ma.masked_where(psd_study.T == 0, psd_study.T)).filled(np.nan)[:, :]
    
    psd_diff = nc_out.createVariable('psd_diff', 'f8', ('lat', 'wavenumber'), zlib=True)
    psd_diff.units = 'm2/km'
    psd_diff.coordinates = "wavenumber lat lon"
    psd_diff.long_name = "power spectrum density of difference study minus reference field"
    psd_diff[:, :] = np.ma.masked_invalid(np.ma.masked_where(psd_diff_ref_study.T == 0, psd_diff_ref_study.T))[:, :]
    
    # resolution = compute_resolution(vlon, vlat, wavenumber[positive_index], psd_diff_ref_study[positive_index, :, :], psd_ref[positive_index, :, :])
    # resolution_out = nc_out.createVariable('effective_resolution', 'f8', ('lat'), zlib=True)
    # resolution_out.units = 'km'
    # resolution_out.coordinates = "lat lon"
    # resolution_out.long_name = "Effective resolution computed as wavelenght where psd_err/psd_ref=0.5"
    # resolution_out[:, :] = np.ma.masked_invalid(np.ma.masked_where(resolution == 0, resolution)).filled(np.nan)

    coherence_out = nc_out.createVariable('coherence', 'f8', ('lat', 'wavenumber'), zlib=True)
    coherence_out[:, :] = np.ma.masked_invalid(np.ma.masked_where(coherence.T == 0, coherence.T)).filled(np.nan)[:, :]
    coherence_out.coordinates = "wavenumber lat lon"
    coherence_out.long_name = "magnitude squared coherence between reference and study fields"

    cross_spectrum_real_out = nc_out.createVariable('cross_spectrum_real', 'f8', ('lat', 'wavenumber'),
                                                        zlib=True)
    cross_spectrum_real_out[:, :] = np.ma.masked_invalid(np.ma.masked_where(np.real(cross_spectrum.T) == 0., np.real(cross_spectrum.T))).filled(np.nan)[:, :]
    cross_spectrum_real_out.coordinates = "wavenumber lat lon"
    cross_spectrum_real_out.long_name = "real part of cross_spectrum between reference and study fields"
    cross_spectrum_imag_out = nc_out.createVariable('cross_spectrum_imag', 'f8', ('lat', 'wavenumber'),
                                                        zlib=True)
    cross_spectrum_imag_out[:, :] = np.ma.masked_invalid(np.ma.masked_where(np.imag(cross_spectrum.T) == 0., np.imag(cross_spectrum.T))).filled(np.nan)[:, :]
    cross_spectrum_imag_out.coordinates = "wavenumber lat lon"
    cross_spectrum_imag_out.long_name = "imaginary part of cross_spectrum between reference and study fields"

    nc_out.method=method_name
    
    nc_out.close()

    
    
    
########################################            
def get_dx_dy(data,navlon,navlat):
    """
    Obtain dx and dy from navigation data.

    Parameters
    ----------
    data : numpy.ndarray
        Data.
    navlon : numpy.ndarray
        Longitude.
    navlat : numpy.ndarray
        Latitude.

    Returns
    -------
    dx : float
        Grid spacing in the x-direction (kilometers).
    dy : float
        Grid spacing in the y-direction (kilometers).
    """
    
    # Obtain dx and dy
    x,y,dat = ps.interpolate(data,np.array(navlon),np.array(navlat))
    x1,y1 = x[0,:],y[:,0]
    dx=np.int(np.ceil(x1[1]-x1[0]))
    dy=np.int(np.ceil(y1[1]-y1[0]))
    dx = dx/1E3
    dy = dy/1E3
    return dx,dy        

#########################################
def apply_window(da, dims, window_type='hanning'):
    """
    Create windows in dimensions dims.

    Parameters
    ----------
    da : xarray.DataArray
        Data array.
    dims : list of str
        Dimensions.
    window_type : str, optional
        Windowing type (default is 'hanning').

    Returns
    -------
    da : xarray.DataArray
        Data array with applied windowing.
    """

    if window_type not in ['hanning']:
        raise NotImplementedError("Only hanning window is supported for now.")

    numpy_win_func = getattr(np, window_type)

    if da.chunks:
        def dask_win_func(n):
            return dsar.from_delayed(
                delayed(numpy_win_func, pure=True)(n),
                (n,), float)
        win_func = dask_win_func
    else:
        win_func = numpy_win_func

    windows = [xr.DataArray(win_func(len(da[d])),
               dims=da[d].dims, coords=da[d].coords) for d in dims]

    return da * reduce(operator.mul, windows[::-1])

#########################################
def get_wavnum_kradial(kx,ky):
    """
    Compute a wavenumber vector and radial wavenumber.

    Parameters
    ----------
    kx : numpy.ndarray
        Wavenumber values in the x-direction.
    ky : numpy.ndarray
        Wavenumber values in the y-direction.

    Returns
    -------
    wavnum : numpy.ndarray
        Wavenumber values.
    kradial : numpy.ndarray
        Radial wavenumber values.
    """
    
    k, l = np.meshgrid(kx,ky)
    kradial = np.sqrt(k**2 + l**2)
    kmax = np.sqrt((k.max())**2 + (l.max())**2)/np.sqrt(2)
    
    dkx = np.abs(kx[2]-kx[1])
    dky = np.abs(ky[2]-ky[1])
    dkradial = min(dkx,dky)
    
    # radial wavenumber
    wavnum = (dkradial.data)*np.arange(1,int(kmax/dkradial))
    return wavnum,kradial

#########################################
def get_spec_1D(kradial,wavnum,spec_2D):
    """
    Compute the azimuthal average of the 2D spectrum.

    Parameters
    ----------
    kradial : numpy.ndarray
        Radial wavenumber values.
    wavnum : numpy.ndarray
        Wavenumber values.
    spec_2D : numpy.ndarray
        2D spectrum.

    Returns
    -------
    spec_1D : numpy.ndarray
        Azimuthal average of the 2D spectrum.
    """
    
    spec_1D = np.zeros(len(wavnum))
    for i in range(wavnum.size):
        kfilt =  (kradial>=wavnum[i] - wavnum[0]) & (kradial<=wavnum[i])
        N = kfilt.sum()
        spec_1D[i] = (spec_2D[kfilt].sum())*wavnum[i]/N  #
    return spec_1D

#########################################
def get_f_kx_ky(hat):
    """
    Get frequency, kx, and ky from a hat object.

    Parameters
    ----------
    hat : hat object
        The hat object.

    Returns
    -------
    f : numpy.ndarray
        Frequency values.
    kx : numpy.ndarray
        Wavenumber values in the x-direction.
    ky : numpy.ndarray
        Wavenumber values in the y-direction.
    """
    
    f = hat.f_time_counter
    kx = hat.f_x
    ky = hat.f_y
    return f,kx,ky

#########################################
def get_f_kx_ky_mit(hat):
    """
    Get frequency, kx, and ky from a MITgcm hat object.

    Parameters
    ----------
    hat : MITgcm hat object
        The MITgcm hat object.

    Returns
    -------
    f : numpy.ndarray
        Frequency values.
    kx : numpy.ndarray
        Wavenumber values in the x-direction.
    ky : numpy.ndarray
        Wavenumber values in the y-direction.
    """
    
    f = hat.f_time
    kx = hat.f_i
    ky = hat.f_j
    return f,kx,ky

#########################################
def get_f_kx_ky_flo(hat):
    """
    Get frequency, kx, and ky from a FLOWMOS hat object.

    Parameters
    ----------
    hat : FLOWMOS hat object
        The FLOWMOS hat object.

    Returns
    -------
    f : numpy.ndarray
        Frequency values.
    kx : numpy.ndarray
        Wavenumber values in the x-direction.
    ky : numpy.ndarray
        Wavenumber values in the y-direction.
    """    
    
    f = hat.f_time
    kx = hat.f_x
    ky = hat.f_y
    return f,kx,ky

#########################################
def get_f_kx_ky_jet(hat):
    """
    Get frequency, kx, and ky from a JET hat object.

    Parameters
    ----------
    hat : JET hat object
        The JET hat object.

    Returns
    -------
    f : numpy.ndarray
        Frequency values.
    kx : numpy.ndarray
        Wavenumber values in the x-direction.
    ky : numpy.ndarray
        Wavenumber values in the y-direction.
    """
    
    f = hat.f_time_counter
    kx = hat.f_x_rho
    ky = hat.f_y_rho
    return f,kx,ky


#########################################
def get_f_k_in_2D(kradial,wavnum,spec2D):
    """
    Compute 2D spectra in frequency-wavenumber space.

    Parameters
    ----------
    kradial : numpy.ndarray
        Radial wavenumber values.
    wavnum : numpy.ndarray
        Wavenumber values.
    spec2D : numpy.ndarray
        2D spectrum.

    Returns
    -------
    spec_1D : numpy.ndarray
        2D spectra in frequency-wavenumber space.
    """    
    
    _spec_1D = []
    for i in range(len(spec2D)):
        #if i%100==0: print(i)
        psd_2D = spec2D[i]
        spec1D = get_spec_1D(kradial,wavnum,psd_2D)
        _spec_1D.append(spec1D)
    spec_1D = np.array(_spec_1D)
    return spec_1D

#########################
def get_flux(wavnum2D,wavnum1D,spec_2D):
    """
    Compute kinetic energy flux.

    Parameters
    ----------
    wavnum2D : numpy.ndarray
        2D wavenumber values.
    wavnum1D : numpy.ndarray
        1D wavenumber values.
    spec_2D : numpy.ndarray
        2D spectrum.

    Returns
    -------
    flux : numpy.ndarray
        Kinetic energy flux.
    """
    
    flux = np.zeros(len(wavnum1D))
    for i in range(wavnum1D.size):
        kfilt =  (wavnum1D[i] <= wavnum2D ) 
        flux[i] = (spec_2D[kfilt]).sum()
    return flux

#########################
def get_flux_in_1D(kradial,wavnum,spec2D):
    """
    Compute 1D kinetic energy flux.

    Parameters
    ----------
    kradial : numpy.ndarray
        Radial wavenumber values.
    wavnum : numpy.ndarray
        Wavenumber values.
    spec2D : numpy.ndarray
        2D spectrum.

    Returns
    -------
    flux_1D : numpy.ndarray
        1D kinetic energy flux.
    """
    
    _flux_1D = []
    for i in range(len(spec2D)):
        #if i%100==0: print(i)
        psd_2D = spec2D[i]
        flux1D = get_flux(kradial,wavnum,psd_2D)
        _flux_1D.append(flux1D)
    flux_1D = np.array(_flux_1D)
    return flux_1D



#########################################
def detrendn(da, axes=None):
    """
    Detrend by subtracting out the least-square plane or least-square cubic fit
    depending on the number of axis.
    
    Parameters
    ----------
    da : `dask.array`
        The data to be detrended
    Returns
    -------
    da : `numpy.array`
        The detrended input data
    """
#     if da.ndim > 2:
#         raise ValueError('The data should only have two dimensions')
#     print(da.shape)
    N = [da.shape[n] for n in axes]
    M = []
    for n in range(da.ndim):
        if n not in axes:
            M.append(da.shape[n])

    if len(N) == 2:
        G = np.ones((N[0]*N[1],3))
        for i in range(N[0]):
            G[N[1]*i:N[1]*i+N[1], 1] = i+1
            G[N[1]*i:N[1]*i+N[1], 2] = np.arange(1, N[1]+1)
        if type(da) == xr.DataArray:
            d_obs = np.reshape(da.copy().values, (N[0]*N[1],1))
        else:
            d_obs = np.reshape(da.copy(), (N[0]*N[1],1))
    elif len(N) == 3:
        if type(da) == xr.DataArray:
            if da.ndim > 3:
                raise NotImplementedError("Cubic detrend is not implemented "
                                         "for 4-dimensional `xarray.DataArray`."
                                         " We suggest converting it to "
                                         "`dask.array`.")
            else:
                d_obs = np.reshape(da.copy().values, (N[0]*N[1]*N[2],1))
        else:
            d_obs = np.reshape(da.copy(), (N[0]*N[1]*N[2],1))

        G = np.ones((N[0]*N[1]*N[2],4))
        G[:,3] = np.tile(np.arange(1,N[2]+1), N[0]*N[1])
        ys = np.zeros(N[1]*N[2])
        for i in range(N[1]):
            ys[N[2]*i:N[2]*i+N[2]] = i+1
        G[:,2] = np.tile(ys, N[0])
        for i in range(N[0]):
            G[len(ys)*i:len(ys)*i+len(ys),1] = i+1
    else:
        raise NotImplementedError("Detrending over more than 4 axes is "
                                 "not implemented.")

    m_est = np.dot(np.dot(spl.inv(np.dot(G.T, G)), G.T), d_obs)
    d_est = np.dot(G, m_est)

    lin_trend = np.reshape(d_est, da.shape)

    return da - lin_trend


def velocity_derivatives(u, v, xdim, ydim, dx):
    """
    Compute the derivatives of velocity fields.

    Parameters
    ----------
    u : array-like
        Velocity component in the x-direction.
    v : array-like
        Velocity component in the y-direction.
    xdim : str
        Name of the x-dimension.
    ydim : str
        Name of the y-dimension.
    dx : float
        Grid spacing.

    Returns
    -------
    ds_derivatives : xarray.Dataset
        Dataset containing the derivatives:
        - 'u_x': Derivative of u in the x-direction.
        - 'u_y': Derivative of u in the y-direction.
        - 'v_x': Derivative of v in the x-direction.
        - 'v_y': Derivative of v in the y-direction.
    """
    
    uhat = xfft.fft(u, dim=(xdim, ydim), dx=dx, sym=True)
    vhat = xfft.fft(v, dim=(xdim, ydim), dx=dx, sym=True)
    k = uhat['f_%s' % xdim]
    l = vhat['f_%s' % ydim]
    u_x_hat = uhat * 2 * np.pi * 1j * k 
    u_y_hat = uhat * 2 * np.pi * 1j * l 
    v_x_hat = vhat * 2 * np.pi * 1j * k
    v_y_hat = vhat * 2 * np.pi * 1j * l
    ds_derivatives = xr.Dataset({'u_x': xfft.ifft(u_x_hat),
                                 'u_y': xfft.ifft(u_y_hat),
                                 'v_x': xfft.ifft(v_x_hat),
                                 'v_y': xfft.ifft(v_y_hat)
                                })
    return ds_derivatives



def compute_wk(lon,lat,ssh):
    """
    Compute the wavenumber-frequency spectrum of sea surface height (SSH) data.

    Parameters
    ----------
    lon : array-like
        Longitudes.
    lat : array-like
        Latitudes.
    ssh : array-like
        Sea surface height data.

    Returns
    -------
    wavenumber : array
        Isotropic wavenumber values.
    frequency : array
        Frequency values.
    SSH_wavenum_freq_spectrum : array
        Two-dimensional wavenumber-frequency spectrum of SSH.
    """
    
    if False:#xr.ufuncs.isnan(ssh).sum().values: 
        print('interpolating nan values')
        import pyinterp
        import pyinterp.backends.xarray 
        import pyinterp.fill 
        
        for itime in range(est1['time'].size):
            ssh0 = ssh[itime,:,:]

            ssh0['lon'].attrs["units"] = "degree_E"
            ssh0['lat'].attrs["units"] = "degree_N"

            interpolator = pyinterp.backends.xarray.Grid2D(ssh0[:,:],
                                                           increasing_axes=True)

            mx, my = numpy.meshgrid(ssh0['lon'],
                                        ssh0['lat'],
                                        indexing='ij')
            has_converged, filled = pyinterp.fill.gauss_seidel(interpolator)
            ssh0 = np.transpose(filled)

            ssh[itime,:,:] = ssh0
 
    
    # - get dx and dy
    #print('get dx and dy')
    dx,dy = get_dx_dy(ssh[0],lon,lat)

    #... Detrend data in all dimension ...
    #print('Detrend data in all dimension')
    ssh_detrended = detrendn(ssh,axes=[0,1,2])

    #... Apply hanning windowing ...') 
    #print('Apply hanning windowing')
    ssh_hanning = apply_window(ssh_detrended, ssh.dims, window_type='hanning')

    #... Apply hanning windowing ...') 
    #print('FFT ')
    ssh_hat = xfft.fft(ssh_hanning, dim=('time', 'lon', 'lat'), dx={'lon': dx, 'lat': dx}, sym=True)
     

    #... Apply hanning windowing ...') 
    #print('PSD ')
    ssh_psd = xfft.psd(ssh_hat)

    #... Get frequency and wavenumber ... 
    #print('Get frequency and wavenumber')
    frequency = ssh_hat.f_time
    kx = ssh_hat.f_lon
    ky = ssh_hat.f_lat 

    #... Get istropic wavenumber ... 
    #print('Get istropic wavenumber')
    wavenumber,kradial = get_wavnum_kradial(kx,ky)

    #... Get numpy array ... 
    #print('Get numpy array')
    ssh_psd_np = ssh_psd.values

    #... Get 2D frequency-wavenumber field ... 
    #print('Get f k in 2D')
    SSH_wavenum_freq_spectrum = get_f_k_in_2D(kradial,wavenumber,ssh_psd_np)
    
    return wavenumber, frequency, SSH_wavenum_freq_spectrum

# def compute_psd_scores_current(ds_interp, output_filename, lenght_scale=np.timedelta64(20, 'D')):
    
#     logging.info('Segment computation...')
#     # delta_t = np.median(np.diff(ds_interp['time']))
#     delta_t = np.median(np.diff(np.unique(ds_interp['time'])))
#     delta_t_days = (delta_t / np.timedelta64(1, 'D')).astype(np.float64)
#     npt = int(lenght_scale / delta_t)
    
#     ds_interp['complex_current_drifter'] = (("time"), ds_interp['EWCT'].data + 1j* ds_interp['NSCT'].data)
#     ds_interp['complex_current_maps'] = (("time"), ds_interp['ugos_interpolated'].data + 1j* ds_interp['vgos_interpolated'].data)
    
#     list_of_lon = []
#     list_of_lat = []
#     list_of_udrifter = []
#     list_of_umaps = []
    
#     for sensor_id in np.unique(ds_interp['sensor_id']):
#         ds_sel = ds_interp.where(ds_interp['sensor_id'] == sensor_id, drop=True)
#         lon_segment, lat_segment, udfriter_segment, umap_segment = compute_segment(ds_sel, npt, 
#                                                                               ref_var_name='complex_current_drifter', 
#                                                                               study_var_name='complex_current_maps',
#                                                                               max_time_gap=np.timedelta64(10, 'h'),
#                                                                               segment_overlapping=0.1)
#         del ds_sel
#         if (udfriter_segment.size > 0) and (umap_segment.size > 0):
#             list_of_lon.append(lon_segment)
#             list_of_lat.append(lat_segment)
#             list_of_udrifter.append(udfriter_segment)
#             list_of_umaps.append(umap_segment)
#         del lon_segment, lat_segment, udfriter_segment, umap_segment
            
#     lon_segment = np.concatenate(list_of_lon)
#     lat_segment = np.concatenate(list_of_lat)
#     udfriter_segment = np.concatenate(list_of_udrifter)
#     umap_segment = np.concatenate(list_of_umaps)
    
#     logging.info('Spectral analysis...')
#     wavenumber, vlat, vlon, nb_segment, psd_ref, psd_study, psd_diff, coherence, cross_spectrum = spectral_computation(lon_segment,
#                                                                                                                        lat_segment,
#                                                                                                                        udfriter_segment,
#                                                                                                                        umap_segment, 
#                                                                                                                        delta_t_days,
#                                                                                                                        npt
#                                                                                                                       )
#     logging.info('Saving ouput...')
#     write_psd_output(output_filename, wavenumber, vlat, vlon, nb_segment, psd_ref, psd_study, psd_diff, coherence, cross_spectrum, one_sided=False)
    
    
# def plot_psd_scores_currents(filename):
    
#     ds_psd = xr.open_dataset(filename)
    
#     ds_clockwise = ds_psd.where(ds_psd['wavenumber'] <= 0., drop=True)
#     ds_counter_clockwise = ds_psd.where(ds_psd['wavenumber'] >= 0., drop=True)

#     ds_clockwise['wavenumber'] = np.abs(ds_clockwise['wavenumber'] )
#     ds_clockwise['wavelenght'] = np.abs(1./ds_clockwise['wavenumber'])
#     ds_clockwise['wavelenght'].attrs['units'] = 'days'
#     ds_clockwise = ds_clockwise.assign_coords(wavelenght=ds_clockwise['wavelenght'])

#     ds_counter_clockwise['wavenumber'] = np.abs(ds_counter_clockwise['wavenumber'] )
#     ds_counter_clockwise['wavelenght'] = np.abs(1./ds_counter_clockwise['wavenumber'])
#     ds_counter_clockwise['wavelenght'].attrs['units'] = 'days'
#     ds_counter_clockwise = ds_counter_clockwise.assign_coords(wavelenght=ds_counter_clockwise['wavelenght'])
    
    
    
#     fig1 = (ds_clockwise.psd_ref.hvplot.line(x='wavelenght', logx=True, logy=True, label='PSD ref (clockwise)', color='k', flip_xaxis=True, line_width=3)*\
#     ds_clockwise.psd_study.hvplot.line(x='wavelenght', logx=True, logy=True, label='PSD study (clockwise)', color='r', flip_xaxis=True)*\
#     ds_clockwise.psd_diff.hvplot.line(x='wavelenght', logx=True, logy=True, label='PSD err (clockwise)', color='grey', flip_xaxis=True)*\
#     ds_counter_clockwise.psd_ref.hvplot.line(x='wavelenght', logx=True, logy=True, label='PSD ref (counterclockwise)', line_dash='dashed', color='k', flip_xaxis=True, line_width=3)*\
#     ds_counter_clockwise.psd_study.hvplot.line(x='wavelenght', logx=True, logy=True, label='PSD study (counterclockwise)', line_dash='dashed', color='r', flip_xaxis=True)*\
#     ds_counter_clockwise.psd_diff.hvplot.line(x='wavelenght', logx=True, logy=True, label='PSD err (counterclockwise)', line_dash='dashed', color='grey', flip_xaxis=True)).opts(legend_position='bottom_left')
    
#     fig2 = (ds_clockwise.coherence.hvplot.line(x='wavelenght', logx=True, logy=False, label='Coherence (clockwise) ', color='k', flip_xaxis=True, line_width=2, ylim=(0, 1))*\
#     (ds_clockwise.psd_diff/ds_clockwise.psd_ref).hvplot.line(x='wavelenght', logx=True, logy=False, label='PSD err/ PSD ref (clockwise)', color='r', flip_xaxis=True, line_width=2, ylim=(0, 1))*\
#     ds_counter_clockwise.coherence.hvplot.line(x='wavelenght', logx=True, logy=False, label='Coherence (counterclockwise) ', color='k', flip_xaxis=True, line_width=2, line_dash='dashed', ylim=(0, 1))*\
#     (ds_counter_clockwise.psd_diff/ds_clockwise.psd_ref).hvplot.line(x='wavelenght', logx=True, logy=False, label='PSD err/ PSD ref (counterclockwise)', color='r', flip_xaxis=True, line_width=2, line_dash='dashed', ylim=(0, 1))).opts(legend_position='bottom_left')
    
#     return (ds_psd.effective_resolution.hvplot.quadmesh(x='lon', y='lat', clim=(100, 500), cmap='jet') +\
# ds_psd.nb_segment.hvplot.quadmesh(x='lon', y='lat', cmap='jet') +\
# fig1+fig2).cols(2) 

