import xarray as xr
import zarr
import shutil
import os
from typing import Dict

class MesoscaleDiagnostics:

    def __init__(self, dataset: xr.Dataset, axis_mapping: Dict[str, str]):
        self._obj = dataset
        self._axis_mapping = axis_mapping

    def store(self, store, variables, **kwargs):
        """Write fields to disk and update the dataset."""
        self._obj.reset_coords()[variables].to_zarr(store, mode="a", **kwargs)
        self._obj = self._obj.update(xr.open_zarr(store)[variables])
        return self

    def clean(self, store, variables, **kwargs):
        """Remove fields from storage.
        
        Use this to delete temporary fields written to disk after the
        final results have been produced.
        """
        for variable in variables:
            try:
                shutil.rmtree(os.path.join(store, variable))
            except FileNotFoundError:
                pass
        
        # Rebuild the zarr metadata
        zarr.consolidate_metadata(store)
        
        self._obj = self._obj[[v for v in self._obj if v not in variables]]
        return self

    def __repr__(self):
        return self._obj.__repr__()

    def _repr_html_(self):
        return self._obj._repr_html_()