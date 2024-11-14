from scipy import signal
import dask.array as da

def _wrapped_sos_filtered(variable, sos, time_axis):
    """Wraps the sos filter function to be compatible with map_blocks."""
    return signal.sosfiltfilt(sos, variable, axis=time_axis)


def filter_butterworth(z, time, cutoff_period = 48, order = 3, time_axis=0):
    """Butterworth filtering over the time dimension.
    
    This function applies a temporal filtering based on the wanted cutoff frequency,
    and a lowpass Butterworth filter on order N

    Parameters
    ----------

    z
        Variable to filter
    time
        Time steps associated with the variable
    cutoff_period
        Cutoff period of the butterworth filter in hours. Default is 48h
    order
        Order of the Butterworth filter. Default is 3
    time_axis
        Time axis index in the variable array

    Returns
    -------
    z_filtered: np.ndarray
        Filtered variable with butterworth filter over the time dimension
    """

    original_chunks = list(z.chunks)
    if len(original_chunks[time_axis]) > 0:
        print(
            "Warning: variable is chunked along the filtered dimension: fusing chunks can degrade performance")
        chunks = list(original_chunks)
        chunks[time_axis] = -1
        z = z.rechunk(chunks)

    # Compute the sampling period in seconds
    sampling_period = (time[1] - time[0]).astype("timedelta64[s]").item().seconds

    # Creation of the temporal filter
    sampling_frequency = 1 / sampling_period
    cutoff_frequency = 1 / (cutoff_period * 3600)

    sos = signal.butter(order, cutoff_frequency, btype='lowpass', output='sos', fs=sampling_frequency)

    # Parallized processing
    return da.map_blocks(
        _wrapped_sos_filtered,
        z,
        sos,
        time_axis,
        dtype=float).rechunk(original_chunks)

