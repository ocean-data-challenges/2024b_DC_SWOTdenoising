import xarray as xr
from typing import List, Union


def add_dataset_to_store(
    store: str,
    ds: xr.Dataset,
    variables: List[str] = [],
    ignore_coords: bool=True,
    **kwargs) -> xr.Dataset:


    if len(variables) == 0:
        variables = [v for v in ds]

    ds_old = ds
    if ignore_coords:
        ds = ds.reset_coords()


    ds[variables].to_zarr(store, mode="a", **kwargs)

    stored = xr.open_zarr(store)[variables]
    return ds_old.update(stored)


def add_array_to_store(store: str, da: xr.DataArray, **kwargs) -> xr.DataArray:
    try:
        variable_name = da.name
    except ValueError:
        raise ValueError("Input dataarray must have a valid name to be in the dataset !")

    ds = da.to_dataset(name=variable_name)
    
    stored = add_dataset_to_store(store, ds, variables=[variable_name], **kwargs)
    return stored[variable_name]

