import glob
from typing import Any, Dict, List, Optional

import xarray as xr
import numpy as np

from ocean_tools.utilities.area import ROI


def open_llc(grid_store, *vars_store, roi: ROI.GLOBAL):
    """Open data that is split between multiple ZARR stores.
    
    LLC data has been split in multiple stores. It is not optimal
    but we can work with that.

    Some variables are only computed on the equilibrium period (spinup omitted),
    which cut the time dimension (SSU.time.size < Eta.time.size). The dataset will
    be returned with the outer join (nan where values are missing).
    """
    ds = xr.open_zarr(grid_store)
    for var_store in vars_store:
        # Override will skip comparing the coordinates to speed up merging
        ds = xr.merge([ds, xr.open_zarr(var_store)],
                      compat="override",
                      join="outer")

    return crop_over_region(ds, roi)


def crop_over_region(ds: xr.Dataset, roi: ROI.GLOBAL) -> xr.Dataset:
    """Crop an LLC dataset over a given region.

    Crop must take into account the two set of spatial axes (i, j)
    and (i_g, j_g).

    Parameters
    ----------
    ds
        Input LLC4320 dataset
    roi
        Region of interest
    """
    if roi != ROI.GLOBAL:
        ds_where = ds.where(
            np.logical_and(
                np.logical_and(ds.XC >= roi.area[0], ds.XC <= roi.area[2]),
                np.logical_and(ds.YC >= roi.area[1], ds.YC <= roi.area[3]))
        , drop=True)

        i_min, i_max = ds_where.i.data[0], ds_where.i.data[-1]
        j_min, j_max = ds_where.j.data[0], ds_where.j.data[-1]

        ds = ds.isel(
            face=ds_where.face.data,
            i=slice(i_min, i_max+1),
            j=slice(j_min, j_max+1),
            i_g=slice(i_min, i_max+1),
            j_g=slice(j_min, j_max+1))

    return ds
    

def delete_chunk_encoding(ds):
    """Delete `chunks` encoding.
    
    Mainly useful when rechunking the zarr storage
    """
    C = []
    if isinstance(ds, xr.Dataset):
        for c in ds:
            if "chunks" in ds[c].encoding:
                del ds[c].encoding["chunks"]
                C.append(c)
    elif isinstance(ds, xr.DataArray):
        if "chunks" in ds.encoding:
            del ds.encoding["chunks"]
            C.append(ds.name)
    #
    for c in ds.coords:
        if "chunks" in ds[c].encoding:
            del ds[c].encoding["chunks"]
            C.append(c)
    if C:
        print("enconding[\"chunks\"] deleted for: "+"/".join(C))
    return ds