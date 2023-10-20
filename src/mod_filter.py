import numpy as np
import pyinterp
from numba import njit


@njit(cache=True, fastmath=True)
def distance(lon0, lat0, lon1, lat1):
    """
    Compute distance between two geographical locations.

    :param float lon0: Longitude of first location.
    :param float lat0: Latitude of first location.
    :param float lon1: Longitude of second location.
    :param float lat1: Latitude of second location.
    :return: distance (in m)
    :rtype: array
    """
    D2R = np.pi / 180.0
    sin_dlat = np.sin((lat1 - lat0) * 0.5 * D2R)
    sin_dlon = np.sin((lon1 - lon0) * 0.5 * D2R)
    cos_lat1 = np.cos(lat0 * D2R)
    cos_lat2 = np.cos(lat1 * D2R)
    a_val = sin_dlon**2 * cos_lat1 * cos_lat2 + sin_dlat**2
    return 6370997.0 * 2 * np.arctan2(a_val**0.5, (1 - a_val) ** 0.5)


@njit(cache=True)
def median_filter(half_window, x, z):
    """
    Apply a median filter on z field

    :param float half_window: half window where apply median
    :param array x: must be growing for each track but could be irregular
    :param array z: field to apply median
    """
    nb = z.shape[0]
    z_new = np.empty(z.shape, dtype=z.dtype)
    i_previous, i_next = 0, 0
    for i in range(nb):
        while x[i] - x[i_previous] > half_window:
            i_previous += 1
        while i_next < nb and x[i_next] - x[i] <= half_window:
            i_next += 1
        z_new[i] = np.median(z[i_previous:i_next])
    return z_new


@njit(cache=True)
def lanczos_filter(wave_length, x, z, order=1):
    """
    Apply a lanczos filter on z field

    :param float wave_length: half window where apply lanczos in x units
    :param array x: must be growing for each track but could be irregular
    :param array z: field to apply lanczos
    """
    nb = z.shape[0]
    last = nb - 1
    z_new = np.empty(z.shape, dtype=z.dtype)
    for i in range(nb):
        z_sum = z[i]
        w_sum = 1
        if i != 0:
            # from the computed value to the left bounds of window
            i_previous = i - 1
            dx = (x[i] - x[i_previous]) / wave_length
            while dx < order and i_previous >= 0:
                w = order * np.sin(np.pi * dx) * np.sin(np.pi * dx / order) / (np.pi * dx) ** 2
                z_sum += z[i_previous] * w
                w_sum += w
                i_previous -= 1
                dx = (x[i] - x[i_previous]) / wave_length
        if i != last:
            # from the computed value to the left bounds of window
            i_next = i + 1
            dx = (x[i_next] - x[i]) / wave_length
            while dx < order and i_next != last:
                w = order * np.sin(np.pi * dx) * np.sin(np.pi * dx / order) / (np.pi * dx) ** 2
                z_sum += z[i_next] * w
                w_sum += w
                i_next += 1
                dx = (x[i_next] - x[i]) / wave_length
        z_new[i] = z_sum / w_sum
    return z_new


def compute_median_dx(dataset):
    """
    Compute the median spacing between along-track measurements.

    Parameters
    ----------
    dataset : xarray.Dataset
        Input dataset containing longitude and latitude coordinates.

    Returns
    -------
    float
        The median spacing between along-track measurements in kilometers.
    """
        
    return 0.001*np.median(pyinterp.geodetic.coordinate_distances(dataset['longitude'][:-1].values,
                                                              dataset['latitude'][:-1].values,
                                                              dataset['longitude'][1:].values,
                                                              dataset['latitude'][1:].values
                                                             ))


def apply_bandpass_filter(ds, lambda_min=65., lambda_max=500.):
    """
    Apply a bandpass filter to a dataset.

    Parameters
    ----------
    ds : xarray.Dataset
        Input dataset containing relevant variables.
    lambda_min : float, optional
        Minimum wavelength for the filter in kilometers, by default 65.
    lambda_max : float, optional
        Maximum wavelength for the filter in kilometers, by default 500.

    Returns
    -------
    xarray.Dataset
        The filtered dataset with additional variables 'msla_filtered', 'sla_filtered', and 'mapping_err_filtered'.
    """
    
    # Compute median spacinf between along-track measurement (in km)
    dx = compute_median_dx(ds)
    
    # Filter scales in map of SLA
    filter_lambda_max = lanczos_filter(lambda_max/dx*np.timedelta64(1,'s'), ds['time'].values, ds['msla_interpolated'].values)
    filter_lambda_min = lanczos_filter(lambda_min/dx*np.timedelta64(1,'s'), ds['time'].values, ds['msla_interpolated'].values)
    ds['msla_filtered'] = ('time', filter_lambda_min - filter_lambda_max)
    
    # Filter scales in alongtrack
    sla = ds['sla_unfiltered'] - ds['lwe']  # we remove lwe to be independent from DUACS processing !!!!
    filter_lambda_max = lanczos_filter(lambda_max/dx*np.timedelta64(1,'s'), ds['time'].values, sla.values)
    filter_lambda_min = lanczos_filter(lambda_min/dx*np.timedelta64(1,'s'), ds['time'].values, sla.values)
    ds['sla_filtered'] = ('time', filter_lambda_min - filter_lambda_max)
    
    # compute mapping error on filtered fields
    ds['mapping_err_filtered'] = ds['msla_filtered'] - ds['sla_filtered']
    
    return ds