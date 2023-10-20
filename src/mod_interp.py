import datetime
import logging
from datetime import timedelta

import numpy as np 
import pandas
import pyinterp
import pyinterp.backends.xarray
import xarray as xr


class TimeSeries:
    """
    Manage a time series composed of a grid stack.

    Parameters
    ----------
    ds : xarray.Dataset
        Input dataset containing the time series data.

    Attributes
    ----------
    ds : xarray.Dataset
        Input dataset containing the time series data.
    series : pandas.Series
        Time series data loaded from the dataset.
    dt : datetime.timedelta
        Time step duration between consecutive data points in the series.

    Methods
    -------
    _is_sorted(array)
        Check if an array is sorted.
    _load_ts()
        Load the time series data into memory.
    _load_dataset(self, varname, start, end)
        Loading the time series into memory for the defined period.
    """

    def __init__(self, ds):
        """
        Initialize a TimeSeries object.

        Parameters
        ----------
        ds : xarray.Dataset
            Input dataset containing the time series data.
        """
        
        self.ds = ds
        self.series, self.dt = self._load_ts()

    @staticmethod
    def _is_sorted(array):
        """
        Check if an array is sorted.

        Parameters
        ----------
        array : numpy.ndarray
            Input array to check.

        Returns
        -------
        bool
            True if the array is sorted, False otherwise.
        """
        
        indices = np.argsort(array)
        return np.all(indices == np.arange(len(indices)))

    def _load_ts(self):
        """
        Load the time series data into memory.

        Returns
        -------
        pandas.Series
            Loaded time series data.
        datetime.timedelta
            Time step duration between consecutive data points in the series.
        """
        time = self.ds.time
        assert self._is_sorted(time)

        series = pandas.Series(time)
        frequency = set(
            np.diff(series.values.astype("datetime64[s]")).astype("int64"))
        if len(frequency) != 1:
            raise RuntimeError(
                "Time series does not have a constant step between two "
                f"grids: {frequency} seconds")
        #return series, datetime.timedelta(seconds=float(frequency.pop()))
        return series, timedelta(seconds=float(frequency.pop()))

    def _load_dataset(self, varname, start, end):
        """
        Loading the time series into memory for the defined period.
        
        Parameters
        ----------
        varname: str
            Name of the variable to be loaded into memory.
        start: datetime.datetime
               Date of the first map to be loaded.
        end: datetime.datetime
               Date of the last map to be loaded.
        
        Returns
        -------  
        pyinterp.backends.xarray.Grid3D: 
                The interpolator handling the interpolation of the grid series.
        """
        
        if start < self.series.min():
            start = self.series.min()
        if end > self.series.max():
            end = self.series.max()
        
        #if start < self.series.min() or end > self.series.max():
        #    raise IndexError(
        #        f"period [{start}, {end}] out of range [{self.series.min()}, "
        #        f"{self.series.max()}]")
            
        first = start - self.dt
        last = end + self.dt

        selected = self.series[(self.series >= first) & (self.series < last)]
        logging.info("fetch data from %s to %s",selected.min(), selected.max())

        data_array = self.ds[varname].isel(time=selected.index)
        return pyinterp.backends.xarray.Grid3D(data_array)
    
    
def periods(df, time_series, var_name="sla_unfiltered", frequency='W'):
    """
    Return the list of periods covering the time series loaded in memory.

    Parameters
    ----------
    df : pandas.DataFrame
        Input DataFrame containing time series data.
    time_series : TimeSeries
        Time series data and properties.
    var_name : str, optional
        Name of the variable to consider, by default "sla_unfiltered".
    frequency : str, optional
        Frequency for period grouping, by default 'W' (weekly).

    Yields
    ------
    tuple
        A tuple containing the start and end timestamps of each period.
    """
    period_start = df.groupby(
        df.index.to_period(frequency))[var_name].count().index

    for start, end in zip(period_start, period_start[1:]):
        start = start.to_timestamp()
        if start < time_series.series[0]:
            start = time_series.series[0]
        end = end.to_timestamp()
        yield start, end
    yield end, df.index[-1] + time_series.dt
    
    

def interpolate(df, time_series, start, end, var='sla'):
    """
    Interpolate the time series over the defined period.

    Parameters
    ----------
    df : pandas.DataFrame
        Input DataFrame containing time series data.
    time_series : TimeSeries
        Time series data and properties.
    start : pandas.Timestamp
        Start timestamp of the interpolation period.
    end : pandas.Timestamp
        End timestamp of the interpolation period.
    """
    interpolator = time_series._load_dataset(var, start, end)
    mask = (df.index >= start) & (df.index < end)
    selected = df.loc[mask, ["longitude", "latitude"]]
    df.loc[mask, ["msla_interpolated"]] = interpolator.trivariate(
        dict(longitude=selected["longitude"].values,
             latitude=selected["latitude"].values,
             time=selected.index.values),
        interpolator="inverse_distance_weighting",
        num_threads=0)
    
    if var == 'ssh':
        df.msla_interpolated = df.msla_interpolated - df.mdt
    
    
def run_interpolation(ds_maps, ds_alongtrack, frequency='M', var='sla'):
    """
    Interpolate time series data over specified periods.

    Parameters
    ----------
    ds_maps : xarray.Dataset
        Input dataset containing maps data.
    ds_alongtrack : xarray.Dataset
        Input dataset containing along-track data.
    frequency : str, optional
        Frequency for period grouping, by default 'M' (monthly).

    Returns
    -------
    xarray.Dataset
        Interpolated dataset.
    """
    
    time_series = TimeSeries(ds_maps)
    
    df = ds_alongtrack.to_dataframe()

    for start, end in periods(df, time_series, frequency=frequency):
        interpolate(df, time_series, start, end, var)
        
    ds = df.to_xarray()
        
    return ds


def interpolate_current(df, time_series, start, end):
    """
    Interpolate the current time series over the defined period.

    Parameters
    ----------
    df : pandas.DataFrame
        Input DataFrame containing time series data.
    time_series : TimeSeries
        Time series data and properties.
    start : pandas.Timestamp
        Start timestamp of the interpolation period.
    end : pandas.Timestamp
        End timestamp of the interpolation period.
    """
    interpolator = time_series._load_dataset("ugos", start, end)
    mask = (df.index >= start) & (df.index < end)
    selected = df.loc[mask, ["longitude", "latitude"]]
    df.loc[mask, ["ugos_interpolated"]] = interpolator.trivariate(
        dict(longitude=selected["longitude"].values,
             latitude=selected["latitude"].values,
             time=selected.index.values),
        interpolator="inverse_distance_weighting",
        num_threads=0)
    
    interpolator = time_series._load_dataset("vgos", start, end)
    mask = (df.index >= start) & (df.index < end)
    selected = df.loc[mask, ["longitude", "latitude"]]
    df.loc[mask, ["vgos_interpolated"]] = interpolator.trivariate(
        dict(longitude=selected["longitude"].values,
             latitude=selected["latitude"].values,
             time=selected.index.values),
        interpolator="inverse_distance_weighting",
        num_threads=0)


def reformat_drifter_dataset(ds):
    """
    Reformat a drifter dataset, extracting relevant variables.

    Parameters
    ----------
    ds : xarray.Dataset
        Input drifter dataset.

    Returns
    -------
    xarray.Dataset
        Reformatted drifter dataset.
    """
    
    ds = ds.isel(DEPTH=1)
    lat = ds['LATITUDE'].values
    lon = ds['LONGITUDE'].values
    drop_vars = ['TIME_QC', 'POSITION_QC', 'DEPH_QC', 'EWCT_QC', 'NSCT_QC', 'EWCT_WS_QC', 'NSCT_WS_QC', 'WS_TYPE_OF_PROCESSING', 'TEMP', 'TEMP_QC', 'LONGITUDE', 'LATITUDE']
    ds = ds.drop_vars(drop_vars)
    ds = ds.rename({'TIME':'time'})
    ds['longitude'] = (("time"), lon)
    ds['latitude'] = (("time"), lat)
    ds['sensor_id'] = (("time"), ds.platform_code*np.ones(ds['time'].size))
    return ds
 

def run_interpolation_drifters(ds_maps, ds_drifter, time_min, time_max, frequency='M'):
    """
    Interpolate drifters data over specified periods.

    Parameters
    ----------
    ds_maps : xarray.Dataset
        Input dataset containing maps data.
    ds_drifter : xarray.Dataset
        Input dataset containing drifter data.
    time_min : numpy.datetime64
        Minimum time for interpolation.
    time_max : numpy.datetime64
        Maximum time for interpolation.
    frequency : str, optional
        Frequency for period grouping, by default 'M' (monthly).

    Returns
    -------
    xarray.Dataset
        Interpolated drifter dataset.
    """
 
    time_series = TimeSeries(ds_maps)
    
    # Convert to dataframe and interpolate
    df = ds_drifter.to_dataframe()
    for start, end in periods(df, time_series, var_name='NSCT', frequency=frequency):
        interpolate_current(df, time_series, start, end)
        
    ds = df.to_xarray()
    
    return ds


def interp2d(ds,name_vars,lon_out,lat_out):
    """
    Interpolate 2D data on a new grid defined by lon_out and lat_out.

    Parameters
    ----------
    ds : xarray.Dataset
        Dataset containing the data to be interpolated.
    name_vars : dict
        A dictionary specifying the variable names and dimension names.
        Example: {'lon': 'longitude', 'lat': 'latitude', 'var': 'data_variable'}
    lon_out : numpy.ndarray
        2D array of longitudes for the output grid.
    lat_out : numpy.ndarray
        2D array of latitudes for the output grid.

    Returns
    -------
    var_out : numpy.ndarray
        2D array of the interpolated data on the new grid.
    """
    
    from scipy import interpolate
    ds = ds.assign_coords(
                 {name_vars['lon']:ds[name_vars['lon']], # or {name_vars['lon']:(ds[name_vars['lon']] % 360), be careful of that!
                  name_vars['lat']:ds[name_vars['lat']]})
            
    if ds[name_vars['var']].shape[0]!=ds[name_vars['lat']].shape[0]:
        ds[name_vars['var']] = ds[name_vars['var']].transpose()
        
    if len(ds[name_vars['lon']].shape)==2:
        dlon = (ds[name_vars['lon']][:,1:].values - ds[name_vars['lon']][:,:-1].values).max()
        dlat = (ds[name_vars['lat']][1:,:].values - ds[name_vars['lat']][:-1,:].values).max()

        
        ds = ds.where((ds[name_vars['lon']]<=lon_out.max()+dlon) &\
                      (ds[name_vars['lon']]>=lon_out.min()-dlon) &\
                      (ds[name_vars['lat']]<=lat_out.max()+dlat) &\
                      (ds[name_vars['lat']]>=lat_out.min()-dlat),drop=True)
            
        lon_sel = ds[name_vars['lon']].values
        lat_sel = ds[name_vars['lat']].values
            
    else: 
        dlon = (ds[name_vars['lon']][1:].values - ds[name_vars['lon']][:-1].values).max()
        dlat = (ds[name_vars['lat']][1:].values - ds[name_vars['lat']][:-1].values).max()
         
        
        ds = ds.where((ds[name_vars['lon']]<=lon_out.max()+dlon) &\
                      (ds[name_vars['lon']]>=lon_out.min()-dlon) &\
                      (ds[name_vars['lat']]<=lat_out.max()+dlat) &\
                      (ds[name_vars['lat']]>=lat_out.min()-dlat),drop=True)
            
        lon_sel,lat_sel = np.meshgrid(
            ds[name_vars['lon']].values,
            ds[name_vars['lat']].values)
    
    var_sel = ds[name_vars['var']].values

    # Interpolate to state grid   
    var_out = interpolate.griddata((lon_sel.ravel(),lat_sel.ravel()),
                   var_sel.ravel(),
                   (lon_out.ravel(),lat_out.ravel())).reshape((lat_out.shape))
    
    return var_out