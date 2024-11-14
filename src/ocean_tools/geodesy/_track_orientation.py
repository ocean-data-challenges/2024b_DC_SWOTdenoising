import numpy as np
from pyinterp.geodetic import System
from typing import Dict

from ocean_tools.utilities.reshape import slice_along_axis

def track_orientation(
    latitude: np.ndarray,
    longitude: np.ndarray,
    along_track_axis: int = 0,
    half_width: int = 1,
    spheroid: System = System()):
    """ Determine angle of satellite track with respect the meridian passing the track.
    
    This method relies on the approximation of the track direction using neighbour points.
    The better the localisation, the better precision for the angle. When latitudes and longitudes
    are not very robust, it is possible to increase the half-width to smoothen the speed direction.
    
    SWOT remark: there will be a field computed by the ground segment in the L2 products (although it
    is computed for nadir only)
    
    Parameters
    ----------
    latitude
        Latitudes of the nadir track in degrees
    longitude
        Longitudes of the nadir track in degree
    along_track_axis
        Axis for the along track direction
    half_width
        Half-width of the finite difference calculation.  Set higher
        to smooth the signal if lats and lons are not smooth.
    spheroid
        Earth representation (defaults to WGS84)

    Returns
    -------
    angles_zonal_along: np.ndarray
        Angle between the equator and the along track direction (in radians)
        Positive angles follow the anti-clockwise direction
    """ 
    longitudes = np.radians(longitude)
    latitudes = np.radians(latitude)
    earth_radius = spheroid.mean_radius()

    delta_lon = (
        slice_along_axis(longitudes, along_track_axis, slice(half_width, None)) -
        slice_along_axis(longitudes, along_track_axis, slice(0,-half_width))
    )
    
    # Normalizing the delta_lon between [-pi, pi] will ensure we take the shortest
    # of the two paths available for the distance computation
    # For retrograde orbit (lon goes from 0.5° to 359.5° -> delta_lon = +359° -> -1°)
    delta_lon[delta_lon > np.pi] = delta_lon[delta_lon > np.pi] - 2 * np.pi
    # For prograde orbit (lon goes from 359.5° to 0.5° -> delta_lon = -359° -> +1°)
    delta_lon[delta_lon < -np.pi] = delta_lon[delta_lon < -np.pi] + 2 * np.pi
    
    slice_after = slice_along_axis(latitudes, along_track_axis, slice(half_width, None))
    slice_before = slice_along_axis(latitudes, along_track_axis, slice(0,-half_width))
    delta_lat = slice_after - slice_before
    dy = earth_radius * delta_lat
    
    dx_before = earth_radius * delta_lon * np.cos(slice_after)
    dx_after = earth_radius * delta_lon * np.cos(slice_before)
    
    # return padded dx and dy
    padding = [(0, 0) for ii in range(latitudes.ndim)]
    padding_before = padding.copy()
    padding_before[along_track_axis] = (half_width, 0)
    padding_after = padding.copy()
    padding_after[along_track_axis] = (0, half_width)
    
    dx = np.pad(dx_before, pad_width=padding_before) + np.pad(dx_after, pad_width=padding_after)
    dy = np.pad(dy, pad_width=padding_before) + np.pad(dy, pad_width=padding_after)

    # This gives the angle relative to the equator. Arctan2 is needed to keep
    # the direction info (direction = sens in french)
    return np.arctan2(dy, dx)


def projection_zonal_meridional(
    v_along: np.ndarray,
    v_across: np.ndarray,
    angles_zonal_along: np.ndarray,
    angle_along_across: float = np.pi / 2,
    ) -> Dict[str, np.ndarray]:
    """ Vector projection from the swath into the zonal/meridional components.

    Parameters
    ----------
    v_along
        Vector component in the along track direction
    v_across
        Vector component in the across track direction
    angles_zonal_along
        Angle between the equator and the along track direction (in radians)
        Positive angles follow the anti-clockwise direction
    angle_along_across
        Angle between the along track and across track directions. A positive
        angle follows the anti-clockwise direction. Defaults to pi/2 which defines
        a direct (along, across) coordinate system

    Returns
    -------
    v_zonal: np.ndarray
        Vector component in the zonal direction
    v_meridional: np.ndarray
        Vector component in the meridional direction
    """
    # Source: (I, J) = (Along, Across)
    # Destination: (i, j) = (zonal, meridional)
    v_zonal, v_meridional = vector_projection(
        v_I=v_along,
        v_J=v_across,
        angles_i_I=angles_zonal_along,
        angle_I_J=angle_along_across,
        angle_i_j=np.pi/2
    )

    return dict(v_zonal=v_zonal, v_meridional=v_meridional)


def projection_track(
    v_zonal: np.ndarray,
    v_meridional: np.ndarray,
    angles_zonal_along: np.ndarray,
    angle_along_across: float = np.pi / 2) -> Dict[str, np.ndarray]:
    """Vector projection from zonal/meridional to along/across coordinates.

    Parameters
    ----------
    v_zonal
        Vector component in the zonal direction
    v_meridional
        Vector component in the meridional direction
    angles_zonal_along
        Angle between the equator and the along track direction (in radians)
        Positive angles follow the anti-clockwise direction
    angle_along_across
        Angle between the along track and across track directions. A positive
        angle follows the anti-clockwise direction. Defaults to pi/2 which defines
        a direct (along, across) coordinate system
    
    Returns
    -------
    v_along: np.ndarray
        Vector component in the along track direction
    v_across: np.ndarray
        Vector component in the across track direction
    """
    angles_along_zonal = -angles_zonal_along

    # Source: (I, J) = (zonal, meridional)
    # Destination: (i, j) = (along, across)
    v_along, v_across = vector_projection(
        v_I=v_zonal,
        v_J=v_meridional,
        angles_i_I=angles_along_zonal,
        angle_I_J=np.pi/2,
        angle_i_j=angle_along_across
    )

    return dict(v_along=v_along, v_across=v_across)

def vector_projection(v_I, v_J, angles_i_I, angle_I_J, angle_i_j):
    """Project a vector from (I, J) to (i, j) coordinates.

    v_I
        Vector component over the I direction
    v_J
        Vector component over the J direction
    angles_i_I
        Angles between (i, I) (radians)
    angle_I_J
        Angle between (I, J) (radians). Unique value to deduce the J
        axis from the I axis (+-pi/2)
    angle_i_j
        Angle between (i, j) (radians). Used to determine if the 
        (i, j) coordinate system is direct or not (+pi/2 or -pi/2)

    Returns
    -------
    v_i
        Vector component over the i direction
    v_j
        Vector component over the j direction
    """
    v_i = v_I * np.cos(angles_i_I) + v_J * np.cos(angles_i_I + angle_I_J)
    v_j = (
        np.sign(angle_i_j) * 
        (v_I * np.sin(angles_i_I) + v_J * np.sin(angles_i_I + angle_I_J))
    )
    return v_i, v_j
