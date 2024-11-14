import logging
from enum import Enum, auto
import numpy as np
import cartopy.feature as cfeature
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
from typing import Optional

logger = logging.getLogger(__name__)


class SwathSelection(Enum):
    DESCENDING = auto()
    ASCENDING = auto()
    ALL = auto()

    
def plot_swaths_on_cycle(
    ds,
    variable_name,
    cycle,
    vmin=None,
    vmax=None,
    grid_labels=[True, True, False, False],
    mask: Optional[str] = "ancillary_surface_classification_flag",
    selection=SwathSelection.ALL,
    pass_variable_name="pass_number",
    cycle_variable_name="cycle_number",
    cmap="bwr",
):
    # Select cycle
    cycle_selection = ds[cycle_variable_name] == cycle
    if not np.any(cycle_selection):
        raise Exception(f"Selection over cycle {cycle} is empty")
    ds = ds.where(cycle_selection, drop=True).squeeze()       

    # Prepare figure and create the basemap
    fig, ax = plt.subplots(1, 1, figsize=(7, 4), subplot_kw=dict(projection=ccrs.PlateCarree()))
    ax.set_extent([
        ds.longitude.min().values,
        ds.longitude.max().values,
        ds.latitude.min().values,
        ds.latitude.max().values,
    ])
    ax.add_feature(cfeature.LAND)
    ax.add_feature(cfeature.COASTLINE)  


    # Keep only ascending or descending swath if needed
    if selection == SwathSelection.ASCENDING:
        ds = ds.where(ds[pass_variable_name] % 2 == 1, drop=True)
    if selection == SwathSelection.DESCENDING:
        ds = ds.where(ds[pass_variable_name] % 2 == 0, drop=True)

    # Apply mask
    if mask is not None:
        ds = ds.where(ds[mask] == 0)
        
    # Auto-detect colorange
    vmin = ds[variable_name].min().values if vmin is None else vmin
    vmax = ds[variable_name].max().values if vmax is None else vmax


    for p, dd in ds.groupby(pass_variable_name):
        pcm = ax.pcolormesh(
            dd["longitude"].values,
            dd["latitude"].values,
            dd[variable_name].values,
            vmin=vmin,
            vmax=vmax,
            cmap=cmap,
            rasterized=True)
        
    # Set colorbar with proper size
    im_ratio = (
        (ds.longitude.max().values - ds.longitude.min().values) /
        (ds.latitude.max().values - ds.latitude.min().values)
    )
    plt.colorbar(pcm, ax=ax, fraction=0.045/im_ratio)
    
    # Normalize longitudes
    lons = ds["longitude"].values
    lons[lons > 180] -= 360
    
    # Add grid lines
    gl = ax.gridlines(draw_labels=True, linewidth=0.5, color='gray', alpha=0.2, linestyle='--')
    (gl.left_labels, gl.bottom_labels, gl.right_labels, gl.top_labels) = grid_labels

    # Add title
    try:
        short_name = ds[variable_name].attrs["short_name"]
    except KeyError:
        short_name = variable_name

    try:
        units = ds[variable_name].attrs["units"]
    except KeyError:
        logger.warning(f"Could not determine units for {variable_name} from its attributes")
        units = "Unknown unit"

    ax.title.set_text(f'{short_name} ({units})')
