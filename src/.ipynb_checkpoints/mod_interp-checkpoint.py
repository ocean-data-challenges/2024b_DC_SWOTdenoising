import numpy
import xarray as xr
import logging

import pyinterp
import pyinterp.backends.xarray

import pandas
import datetime
from datetime import timedelta



class TimeSeries:
    """Manage a time series composed of a grid stack."""

    def __init__(self, ds):
        self.ds = ds
        self.series, self.dt = self._load_ts()

    @staticmethod
    def _is_sorted(array):
        indices = numpy.argsort(array)
        return numpy.all(indices == numpy.arange(len(indices)))

    def _load_ts(self):
        """Loading the time series into memory."""
        time = self.ds.time
        assert self._is_sorted(time)

        series = pandas.Series(time)
        frequency = set(
            numpy.diff(series.values.astype("datetime64[s]")).astype("int64"))
        if len(frequency) != 1:
            raise RuntimeError(
                "Time series does not have a constant step between two "
                f"grids: {frequency} seconds")
        #return series, datetime.timedelta(seconds=float(frequency.pop()))
        return series, timedelta(seconds=float(frequency.pop()))

    def load_dataset(self, varname, start, end):
        """Loading the time series into memory for the defined period.
        Args:
            varname (str): Name of the variable to be loaded into memory.
            start (datetime.datetime): Date of the first map to be loaded.
            end (datetime.datetime): Date of the last map to be loaded.
        Returns:
            pyinterp.backends.xarray.Grid3D: The interpolator handling the
            interpolation of the grid series.
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
    """Return the list of periods covering the time series loaded in memory."""
    period_start = df.groupby(
        df.index.to_period(frequency))[var_name].count().index

    for start, end in zip(period_start, period_start[1:]):
        start = start.to_timestamp()
        if start < time_series.series[0]:
            start = time_series.series[0]
        end = end.to_timestamp()
        yield start, end
    yield end, df.index[-1] + time_series.dt
    
    

def interpolate(df, time_series, start, end):
    """Interpolate the time series over the defined period."""
    interpolator = time_series.load_dataset("sla", start, end)
    mask = (df.index >= start) & (df.index < end)
    selected = df.loc[mask, ["longitude", "latitude"]]
    df.loc[mask, ["msla_interpolated"]] = interpolator.trivariate(
        dict(longitude=selected["longitude"].values,
             latitude=selected["latitude"].values,
             time=selected.index.values),
        interpolator="inverse_distance_weighting",
        num_threads=0)
    
    
def run_interpolation(ds_maps, ds_alongtrack, frequency='M'):
    
    time_series = TimeSeries(ds_maps)
    
    df = ds_alongtrack.to_dataframe()

    for start, end in periods(df, time_series, frequency=frequency):
        interpolate(df, time_series, start, end)
        
    ds = df.to_xarray()
        
    return ds


def interpolate_current(df, time_series, start, end):
    """Interpolate the time series over the defined period."""
    interpolator = time_series.load_dataset("ugos", start, end)
    mask = (df.index >= start) & (df.index < end)
    selected = df.loc[mask, ["longitude", "latitude"]]
    df.loc[mask, ["ugos_interpolated"]] = interpolator.trivariate(
        dict(longitude=selected["longitude"].values,
             latitude=selected["latitude"].values,
             time=selected.index.values),
        interpolator="inverse_distance_weighting",
        num_threads=0)
    
    interpolator = time_series.load_dataset("vgos", start, end)
    mask = (df.index >= start) & (df.index < end)
    selected = df.loc[mask, ["longitude", "latitude"]]
    df.loc[mask, ["vgos_interpolated"]] = interpolator.trivariate(
        dict(longitude=selected["longitude"].values,
             latitude=selected["latitude"].values,
             time=selected.index.values),
        interpolator="inverse_distance_weighting",
        num_threads=0)


def reformat_drifter_dataset(ds):
    
    ds = ds.isel(DEPTH=1)
    lat = ds['LATITUDE'].values
    lon = ds['LONGITUDE'].values
    drop_vars = ['TIME_QC', 'POSITION_QC', 'DEPH_QC', 'EWCT_QC', 'NSCT_QC', 'EWCT_WS_QC', 'NSCT_WS_QC', 'WS_TYPE_OF_PROCESSING', 'TEMP', 'TEMP_QC', 'LONGITUDE', 'LATITUDE']
    ds = ds.drop_vars(drop_vars)
    ds = ds.rename({'TIME':'time'})
    ds['longitude'] = (("time"), lon)
    ds['latitude'] = (("time"), lat)
    ds['sensor_id'] = (("time"), ds.platform_code*numpy.ones(ds['time'].size))
    return ds


# def run_interpolation_drifters(ds_maps, list_filename_drifters, time_min, time_max, frequency='M'):
    
#     time_series = TimeSeries(ds_maps)
    
#     # Load drifters data file and reformat dataset
#     list_of_ds = []
#     for drifters_filename in list_filename_drifters:
#         ds_drifter = xr.open_dataset(drifters_filename)
#         ds_drifter = reformat_drifter_dataset(ds_drifter)
#         try:        
#             ds_drifter = ds_drifter.where((ds_drifter.time >= numpy.datetime64(time_min)) & (ds_drifter.time <=  numpy.datetime64(time_max)), drop=True)
#             list_of_ds.append(ds_drifter)
#         except IndexError:
#             # there is gap in dataset
#             pass
#         del ds_drifter

#     ds_out = xr.concat(list_of_ds, dim='time')
    
#     # Convert to dataframe and interpolate
#     df = ds_out.to_dataframe()
#     for start, end in periods(df, time_series, var_name='NSCT', frequency=frequency):
#         interpolate_current(df, time_series, start, end)
        
#     ds = df.to_xarray()
    
#     return ds

def run_interpolation_drifters(ds_maps, ds_drifter, time_min, time_max, frequency='M'):
    
    time_series = TimeSeries(ds_maps)
    
    # Convert to dataframe and interpolate
    df = ds_drifter.to_dataframe()
    for start, end in periods(df, time_series, var_name='NSCT', frequency=frequency):
        interpolate_current(df, time_series, start, end)
        
    ds = df.to_xarray()
    
    return ds