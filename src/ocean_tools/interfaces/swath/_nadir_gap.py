import pyinterp
import xarray as xr
from typing import Optional, Tuple
import numpy as np

def fill_nadir_gap(
    ds: xr.Dataset,
    variable: str,
    method: str = "pad",
    fill_value: float = 0.0,
    nadir_gap: Optional[Tuple[int, int]] = (-10, 10)):
    """Fill the nadir gap using one of the pre-set methods."""

    # If the dataset is chunked (dask arrays) , ensure we do not have chunking over the pixels.
    # This is necessary when processing the gap
    if len(ds.chunks) > 0 and len(
            ds.chunks["num_pixels"]) > 1:
        raise Exception(
            "Chunks overs num_pixels dimension are forbidden, rechunk the array prior to using calling nadir gap filling.")

    # Ensure the cross-track-distance is in kilometers. We set every variable to None to prevent quantifying
    # all arrays
    pint_variables = {
        v: None
        for v in list(ds.variables) if v != "cross_track_distance"
    }
    cross_track_distance = (ds.cross_track_distance.pint.quantify(
        **pint_variables).pint.to("km").data.magnitude)
    cross_track_distance = cross_track_distance[
        0] if cross_track_distance.ndim > 1 else cross_track_distance
    cross_track_distance = cross_track_distance.compute()

    # We expect the cross track distance not to have missing data. The dataset
    # should have all points, even the gap points that are set to np.nan
    A = np.diff(cross_track_distance)
    if not np.allclose(A, np.median(A)):
        raise Exception(
            "Dataset malformed. The cross_track_distance should be constant (ie. the nadir gap must be present"
            "in dataset with invalid values")

    # Find indexer for gap data
    gap = np.logical_and(nadir_gap[0] < cross_track_distance,
                            cross_track_distance < nadir_gap[1])

    # Fill the gap using one of the available methods
    new_name = variable + "_filled"
    da = ds[variable]
    if method == "pad":
        ds[new_name] = da.copy()
        ds[new_name][:, gap] = fill_value
    elif method == "interp":
        ds[new_name] = xr.map_blocks(_fill_nadir_gap_interp,
                                da, (gap, ),
                                template=da)
    else:
        raise Exception(
            f"Unknown {method=} for filling swath gap at nadir, expected one of (pad, interp)"
        )

    return ds



def _fill_nadir_gap_interp(array: xr.DataArray,
                           gap: np.ndarray) -> xr.DataArray:
    # Fill nadir gap of one part of the swath using a regular grid interpolation
    grid2d = pyinterp.Grid2D(pyinterp.Axis(array.num_lines),
                             pyinterp.Axis(array.num_pixels[~gap]),
                             array[:, ~gap, ...])

    lines, pixels = np.meshgrid(array.num_lines, array.num_pixels[gap])
    result = array.copy()
    result[:, gap,
           ...] = pyinterp.bivariate(grid2d, lines.ravel(),
                                     pixels.ravel()).reshape(lines.shape).T

    return result