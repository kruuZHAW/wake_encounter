from typing import Any, Tuple
import math
import pandas as pd
import pyproj


def xy_to_latlon(
    x: Any,
    y: Any,
    orig_proj: pyproj.Proj = pyproj.Proj("epsg:2056"),
    final_proj: pyproj.Proj = pyproj.Proj("epsg:4326"),
) -> Tuple[Any, Any]:
    """
    Transforms x/y coordinates from one projection to another. By default, the
    original projection is CH1903+(Swissgrid) and the final projection is WGS84.

    Parameters
    ----------
    x : Any
        X-coordinate(s) in the original projection. Accepted are all common numeric
        scalar or array types.
    y : Any
        Y-coordinate(s) in the original projection. Accepted are all common numeric
        scalar or array types.
    orig_proj : pyproj.Proj, optional
        Projection of the provided x/y coordinates, by default pyproj.Proj("epsg:2056")
    final_proj : pyproj.Proj, optional
        Projection to which the provided x/y coordinates will be transormed, by default
        pyproj.Proj("epsg:4326")

    Returns
    -------
    Tuple[Any, Any]
        latitudes and longitudes in the final projection (WGS84)
    """
    # Create transformer
    transformer = pyproj.Transformer.from_proj(
        orig_proj, final_proj, always_xy=True
    )
    # Return transformed coordinates
    return transformer.transform(x, y)


def rotate_and_translate(
    x: float, y: float, angle_deg: float, dx: float, dy: float
) -> Tuple[float, float]:
    """
    Rotates and translates a point in 2D space. The rotation is performed first, then
    the translation.

    Parameters
    ----------
    x : float
        X-coordinate
    y : float
        Y-coordinate
    angle_deg : float
        Angle for the rotation in degrees
    dx : float
        Translation in x-direction
    dy : float
        Translation in y-direction

    Returns
    -------
    Tuple[float, float]
        x and y coordinates after rotation and translation have been applied
    """
    # Convert angle to radians
    angle_rad = math.radians(angle_deg)

    # Apply rotation
    x_rot = x * math.cos(angle_rad) - y * math.sin(angle_rad)
    y_rot = -(x * math.sin(angle_rad) + y * math.cos(angle_rad))

    # Apply translation
    x_trans = x_rot + dx
    y_trans = y_rot + dy
    return x_trans, y_trans


def change_ref_xy(
    row: pd.Series, x: float, y: float, angle_deg: float
) -> pd.Series:
    """
    When provided with a row of a dataframe containing information about the vortice
    positions relative to their initial position, returns the row with added columns
    containing the absolute coordinates of the vortices. This includes the position of
    the vortex cores, as well as the 2s and 3s uncertainties.

    Parameters
    ----------
    row : pd.Series
        Row of a dataframe containing the vortices positions relative to their initial
        position. The dataframe must contain the columns "x", "yl", "yr", "y_2s_r",
        "y_2s_l", "y_3s_r", "y_3s_l".
    x : float
        X-coordinate (CH1903+) of the new reference point (aircraft posision at t=0)
    y : float
        Y-coordinate (CH1903+) of the new reference point (aircraft posision at t=0)
    angle_deg : float
        Angle/heading of the aircraft at t=0

    Returns
    -------
    pd.Series
        Provided row with added columns containing the absolute coordinates of the
        vortices (xl_ref, yl_ref, xr_ref, yr_ref, x_2s_r_ref, y_2s_r_ref, x_2s_l_ref,
        y_2s_l_ref, x_3s_r_ref, y_3s_r_ref, x_3s_l_ref, y_3s_l_ref)
    """
    # Rotate and translate left wake
    xn, yn = rotate_and_translate(row.x, row.yl, angle_deg, x, y)
    row["xl_ref"] = xn
    row["yl_ref"] = yn
    # Rotate and translate right wake
    xn, yn = rotate_and_translate(row.x, row.yr, angle_deg, x, y)
    row["xr_ref"] = xn
    row["yr_ref"] = yn
    # Rotate and translate 2s uncertainties
    xn, yn = rotate_and_translate(row.x, row.y_2s_r, angle_deg, x, y)
    row["x_2s_r_ref"] = xn
    row["y_2s_r_ref"] = yn
    xn, yn = rotate_and_translate(row.x, row.y_2s_l, angle_deg, x, y)
    row["x_2s_l_ref"] = xn
    row["y_2s_l_ref"] = yn
    # Rotate and translate 3s uncertainties
    xn, yn = rotate_and_translate(row.x, row.y_3s_r, angle_deg, x, y)
    row["x_3s_r_ref"] = xn
    row["y_3s_r_ref"] = yn
    xn, yn = rotate_and_translate(row.x, row.y_3s_l, angle_deg, x, y)
    row["x_3s_l_ref"] = xn
    row["y_3s_l_ref"] = yn

    return row
