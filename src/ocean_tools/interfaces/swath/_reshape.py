import xarray as xr
import numpy as np
from typing import Optional

import logging
logger = logging.getLogger(__name__)


def split_cycle_pass(dataset: xr.Dataset, duplicate_coords: bool = False) -> xr.Dataset:
    """
    Takes the default swath dataset arranged with (num_lines, num_pixels)
    dimensions and rearrange it per track and cycle (num_lines, num_pixels, pass, cycle)
    
    The rearrangement means that all cycles have the same number of measures,
    and that the localization of the measures are identical between cycles. This
    method only checks the number of measures and assumes the 'longitude' and 
    'latitude' variables are the same between cycles.

    Moreover, the pass should have the same length in one cycle, else the stack
    will fail

    Parameters 
    ----------
    dataset
        Xarray dataset of a swath. It must contain the 'ncycle' variable in
        order to stack the cycle along a new dimension
    duplicate_coords
        In case we wish to have the coordinates duplicated over one cycle

    Returns
    -------
    :
        xarray dataset with extra 'pass' and 'cycle' dimensions built from the
        splitting of num_lines
    """
    duplicate_coords = "minimal" if not duplicate_coords else "all"

    # While all coordinates are identical for the cycles, the tracks have
    # different lon/lat so the output dataset should be structured as follows:
    # ssh (cycle, pass, num_lines, num_pixels)
    # longitude (pass, num_lines, num_pixels)
    # In case we want to duplicate the coordinates for the cycles, we can output
    # ssh (cycle, pass, num_lines, num_pixels)
    # longitude (cycle, pass, num_lines, num_pixels)
    dataset = split_unique_variable(
        dataset,
        "cycle_number",
        "cycle",
        coords=duplicate_coords,
        remove_last=True,
        egalize_counts=False)
    dataset = split_unique_variable(
        dataset,
        "pass_number",
        "pass",
        coords='all',
        remove_last=True,
        egalize_counts=True)
    return dataset.transpose("cycle", "pass", ...)


def split_unique_variable(
    dataset:xr.Dataset,
    variable: str,
    dimension: str,
    compat: str ='override',
    coords: str ='minimal',
    remove_last: bool = False,
    egalize_counts: bool = False) -> xr.Dataset:
    """
    Takes a swath dataset (num_lines, num_pixels) and split the num_lines
    dimension with a unique variable (usually cycle_number or pass_number)
    
    The last unique batch of variable is discarded if it has not the same size
    of the previous one.
        
    Parameters 
    ----------
    dataset
        Xarray dataset of a swath. It must contain the 'variable'
    variable
        Name of the variable whose unique values will become a new dimension
    dimension
        Name of the dimension associated with the 'variable' unique values
    compat
        xarray concatenation compatibility. It is set to 'override' by default
        to speed up the concatenation
    coords
        xarray concatenation coords. It is set to 'minimal' to speed up the
        concatenation
    remove_last
        If set to True, the last grouped element can be dropped if there are not
        enough element (ex. the last cycle can be incomplete)
    egalize_counts
        If set to True, the tracks can have 1 line of difference before being
        stacked in a dataset. The tracks with one more line will be cropped

    Returns
    -------
    :
        xarray dataset with an extra 'cycle' dimension to enable time-averages
    """
    # In case the input dataset is already formed with a cycle dimension, do
    # nothing. This can happens if the storage backend already gives a stacked
    # dataset
    if dimension in dataset.dims.keys():
        logger.info(f"Dataset has already unique {dimension}, nothing to do")
        return dataset

    # Ensure num_lines is the first axis
    dataset = dataset.transpose("num_lines", ...)

    try:
        # Check dimension numbers are sorted. This is mandatory for performance issued.
        # Moreover, the hypothesis is not too strong because data are stored by
        # chronological order
        if np.any(np.diff(dataset[variable]) < 0):
            raise Exception(
                f"{variable} are not sorted. Input dataset is malformed")
    except KeyError: 
            raise Exception(
                f"The '{variable}' variable is needed for computation of time "
                f"averages. Please add this variable to the dataset.")

    # Check the indexes 
    unique_numbers, indexes, counts = np.unique(
        dataset[variable],
        return_counts=True,
        return_index=True,
        axis=0,
    )

    # It is possible that tracks have one line 
    if egalize_counts:
        mask = (counts - np.min(counts)) == 1
        if mask.sum() > 0:
            logger.warning(f"Cropping last element to egalizing '{variable}' counts")
        counts[mask] -= 1

    # Last cycle is treated differently from the others because it could be
    # still in the process of acquisition. If this is the case, it is ignored
    # and the analysis can continue.    
    if remove_last and counts[-1] != np.max(counts):
        logger.warning(f"Last '{variable}' is incomplete and has been left out")
        unique_numbers = unique_numbers[:-1]              
        indexes = indexes[:-1]
        counts = counts[:-1]
        

    # Split of num_lines should all have the same size except the last one.
    if not np.all(counts[:-1] == counts[0]):
        idx = (counts[:-1] != counts[0]).nonzero()[0][0]
        raise ValueError(
            f"Dataset malformed, cannot split and stack '{variable}':\n"
            f" {variable}: {unique_numbers[0]} ->{counts[0]} elements \n"
            f" {variable}: {unique_numbers[idx]} -> {counts[idx]} elements")
    
        
    # We make sure to use slices because indexing with boolean or a range has
    # poor performance. Slice triggers the fastest xarray pattern
    datasets=[
        dataset.isel(num_lines=slice(ind, ind+count))
        for ind, count in zip(indexes, counts)
    ]

    dataset = xr.concat(datasets, dimension, compat=compat, coords=coords)
    dataset = dataset.assign_coords({dimension : np.unique(unique_numbers)})
    dataset = dataset.drop_vars([variable])

    return dataset
