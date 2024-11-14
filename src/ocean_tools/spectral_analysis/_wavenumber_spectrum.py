import dask.array as da

from scipy.fftpack import fftn, fftfreq
import numpy as np
from scipy import signal
from scipy import ndimage
from scipy.signal import detrend, tukey, get_window

from ._detrend import detrend_2d

def zonal_wavenumber_spectrum(
    data,
    time,
    dx,
    window_x,
    window_time,
    x_axis: int = 0,
    y_axis: int = 1,
    time_axis: int = 2,
    x_unit_factor: float = 1,
    time_unit: str = "h"):
    """ 
    This function is used to create spectra. First, it detrends and tapers the data
    
    Parameters
    ----------
    data
        Data over which we wish to plot the 
    time
        Time of the data
    dx
        Distances between each point along the x axis
    window_x
        Name of the tapering window along the x axis. It relies on scipy.signal.get_window
        so it accepts the same arguments: a simple string "tukey" or a tuple ("tukey", 0.5) if
        additionnal parameters are needed.
    window_time
        Name of the tapering window along the time axis
    x_axis
        Index of the longitude axis
    y_axis
        Index of the latitude axis
    time_axis
        Index of the time axis
    
    Returns
    -------
    zonal_psd
        PSD of the input signal averaged over the latitudes
    frequencies_time
        FRequencies over the time axis. Only the positive frequencies have been kept
    frequencies_x
        Frequencies over the x axis
    """
    # Mask invalid values
    data_masked = da.ma.filled(da.ma.masked_invalid(data), fill_value=0)

    # Detrending
    data_detrended = detrend_2d(data_masked, x_axis, time_axis)

    # Tapering
    tapering_window = np.outer(
        get_window(window_x, data.shape[x_axis], fftbins=True),
        get_window(window_time, data.shape[time_axis], fftbins=True)
    )
    
    tapering_window = np.expand_dims(
        tapering_window,
        [x for x in np.arange(data_detrended.ndim) if x not in (x_axis, time_axis)]
    )

    data_tapered = data_detrended * tapering_window

    # Compute frequencies in cycle per unit
    dx= np.median(dx) * x_unit_factor
    frequencies_x = fftfreq(data_tapered.shape[x_axis], d=dx) 

    dt = np.median(np.diff(time)) / np.timedelta64(1, time_unit)
    frequencies_time = fftfreq(data_tapered.shape[y_axis], d=dt)

    # FFT spectrum and PSD
    tusp = fftn(data_tapered, axes=(x_axis, time_axis))
    psd = tusp * np.conjugate(tusp)

    # PSD normalization
    ld = np.sum(tapering_window ** 2) / dx / dt
    psd /= ld

    # Mean on latitudes (hence the "zonal" qualification)
    zonal_psd = np.nanmean(psd.real, axis=y_axis)

    # Keep only the positive time frequencies
    zonal_psd = zonal_psd[frequencies_time >= 0] * 2

    return frequencies_time, frequencies_x, zonal_psd