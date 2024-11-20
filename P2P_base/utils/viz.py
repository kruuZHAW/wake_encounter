import os
import math
from typing import Tuple
import numpy as np
import pandas as pd
import shapely.geometry as shg
import shapely.affinity as aff
import svgpath2mpl
import plotly.express as px
import plotly.graph_objects as go
from traffic.core import Flight
from . import coords


def add_rwy_offset(flight: Flight) -> Flight:
    """
    Adjusts the runway offset from an absolute value to a relative value. This allows
    to visualise the aircraft position relative to the runway centerline in the yz plot
    and correct the positions of the wakes accordingly.

    Parameters
    ----------
    flight : Flight
        Flight to which the runway offset should be adjusted.

    Returns
    -------
    Flight
        Flight with adjusted runway offset.
    """
    df = flight.data
    ap = flight.landing_airport()
    rwy_n = flight.data.ILS.mode()[0]
    rwy_bearing = ap.runways.data[ap.runways.data.name == rwy_n][
        "bearing"
    ].iloc[0]
    df = df.assign(
        b_diff=lambda df: df.distance * np.radians(rwy_bearing - df.bearing)
    )
    # Correction for passing 0/360 degrees
    df.loc[(df.b_diff > 330) & (df.b_diff <= 360), "b_diff"] = (
        df.loc[(df.b_diff > 330) & (df.b_diff <= 360), "b_diff"] - 360
    )

    return Flight(df)


def visualise_vortex(wake_df: pd.DataFrame) -> go.Figure:
    """
    Visualise the wake position (and uncertaintites) on the y-z plane over time.

    Parameters
    ----------
    wake_df : pd.DataFrame
        DataFrame containing information of a wake class.
    Returns
    -------
    go.Figure
        Figure containing the visualisation.
    """
    # Create dataframe containing info for both vortex cores
    data_l = wake_df[["yl", "zl", "gam_l"]].rename(
        columns={"yl": "y", "zl": "z", "gam_l": "gam"}
    )
    data_l["t"] = wake_df.index
    data_l["side"] = "left"

    data_r = wake_df[["yr", "zr", "gam_r"]].rename(
        columns={"yr": "y", "zr": "z", "gam_r": "gam"}
    )
    data_r["t"] = wake_df.index
    data_r["side"] = "right"

    data = pd.concat([data_l, data_r])
    data["size"] = 30

    # Add information about uncertainty
    for index, row in wake_df.iterrows():
        t = index
        gam = 1
        size = 1
        x_vals_ellipse_2s, y_vals_ellipse_2s = get_ellipse_points(
            row["y_2s_l"],
            row["y_2s_r"],
            row["z_2s_lo"],
            row["z_2s_hi"],
            num_points=100,
        )
        x_vals_ellipse_3s, y_vals_ellipse_3s = get_ellipse_points(
            row["y_3s_l"],
            row["y_3s_r"],
            row["z_3s_lo"],
            row["z_3s_hi"],
            num_points=100,
        )
        df_2s = pd.DataFrame(
            {
                "t": t,
                "gam": gam,
                "y": x_vals_ellipse_2s,
                "z": y_vals_ellipse_2s,
                "side": "uncert_2s",
                "size": size,
            }
        )
        df_3s = pd.DataFrame(
            {
                "t": t,
                "gam": gam,
                "y": x_vals_ellipse_3s,
                "z": y_vals_ellipse_3s,
                "side": "uncert_3s",
                "size": size,
            }
        )
        data = pd.concat([data, df_2s, df_3s])

    # Create figure and plot cores and uncertainty
    fig = px.scatter(
        data,
        x="y",
        y="z",
        color="gam",
        animation_frame="t",
        size="size",
        width=1000,
        range_x=[
            data.y.min() - 0.05 * abs(data.y.min() - data.y.max()),
            data.y.max() + 0.05 * abs(data.y.min() - data.y.max()),
        ],
        range_y=[
            data.z.min() - 0.05 * abs(data.z.min() - data.z.max()),
            data.z.max() + 0.05 * abs(data.z.min() - data.z.max()),
        ],
        range_color=(data.gam.min(), data.gam.max()),
    )

    # Add line for path of left wake
    line_trace_r = go.Scatter(
        x=data_l["y"],
        y=data_l["z"],
        mode="lines",
        line=dict(color="black"),
        showlegend=False,
        opacity=0.3,
    )
    fig.add_trace(line_trace_r)

    # Add line for path of right wake
    line_trace_r = go.Scatter(
        x=data_r["y"],
        y=data_r["z"],
        mode="lines",
        line=dict(color="black"),
        showlegend=False,
        opacity=0.3,
    )
    fig.add_trace(line_trace_r)

    return fig


def wind_direction_arrow(degrees: float) -> Tuple[float, float]:
    """
    Calculate the x and y components of the arrow tail in pixels to draw the arrow in
    the wind direction and with a length of 40 pixels.

    Parameters
    ----------
    degrees : float
        wind direction in degrees

    Returns
    -------
    Tuple[float, float]
        tuple containing the x and y components in pixels to draw the arrow in the wind
        direction and with a length of 40 pixels.
    """
    # Convert wind direction from degrees to radians
    radians = math.radians(degrees)

    # Calculate the x and y components of the arrow tail
    ax = 40 * math.sin(radians)
    ay = -40 * math.cos(radians)

    # Return the x and y components
    return ax, ay


def marker(path_to_file: str = None) -> shg.Polygon:
    """
    Create a shapely polygon representing the svg file provided as input.

    Parameters
    ----------
    name : str, optional
        path to the svg file, by default "data/airport.svg"

    Returns
    -------
    shg.Polygon
        shapely polygon representing the svg
    """

    def to_shapely(mpl: shg.MultiPolygon, simplify: float = 0.5):
        p = shg.MultiPolygon([shg.Polygon(a).simplify(simplify) for a in mpl])
        p = aff.affine_transform(p, [1, 0, 0, -1, 0, 0])
        scale = 1
        p = aff.affine_transform(p, [1, 0, 0, 1, -p.centroid.x, -p.centroid.y])
        return aff.affine_transform(
            p, [scale, 0, 0, scale, -p.centroid.x, -p.centroid.y]
        )

    if path_to_file is None:
        root_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
        path_to_file = os.path.join(root_dir, "data", "airport.svg")

    svgpath = pd.read_xml(path_to_file).loc[0, "d"]
    return to_shapely(svgpath2mpl.parse_path(svgpath).to_polygons())


def marker_mapbox(
    df: pd.Series,
    size: float = 0.01,
    lat: str = "latitude",
    lon: str = "longitude",
    heading: str = "track",
) -> shg.Polygon:
    """
    Creates a shapely polygon representing the aircraft icon at the given location and
    heading.

    Parameters
    ----------
    df : pd.Series
        Pandas series containing aircraft loaction and heading
    size : float, optional
        Sizing of the icon, by default 0.01
    lat : str, optional
        Column name of the series containing latitude info, by default "latitude"
    lon : str, optional
        Column name of the series containing longitude info, by default "longitude"
    heading : str, optional
        Column name of the series containing heading info, by default "track"

    Returns
    -------
    shg.Polygon
        Shapely polygon representing the aircraft icon at the given location and heading
    """

    m = marker()
    p = aff.affine_transform(m, [size, 0, 0, size, df[lon], df[lat]])
    return aff.rotate(p, 360 - df[heading]).geoms[0]


def get_ellipse_points(
    x_min, x_max, y_min, y_max, num_points=100
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Returns an array of x and y points that define an ellipse
    with the given minimum and maximum values on the x and y axes.

    Parameters
    ----------
    x_min : _type_
        Minimum value on the x axis
    x_max : _type_
        Maximum value on the x axis
    y_min : _type_
        Minimum value on the y axis
    y_max : _type_
        Maximum value on the y axis
    num_points : int, optional
        Number of points to generate for the ellipse, by default 100

    Returns
    -------
    Tuple[np.ndarray, np.ndarray]
        Tuple of arrays containing the x and y points that define the ellipse
    """
    # Define the center point of the ellipse
    x_center = (x_min + x_max) / 2
    y_center = (y_min + y_max) / 2

    # Define the radii of the ellipse
    x_radius = (x_max - x_min) / 2
    y_radius = (y_max - y_min) / 2

    # Generate the angles for the ellipse
    angles = np.linspace(0, 2 * np.pi, num_points)

    # Calculate the x and y coordinates for each angle
    x_points = x_center + x_radius * np.cos(angles)
    y_points = y_center + y_radius * np.sin(angles)

    return x_points, y_points


def visualise_interaction_xy(
    leader: Flight,
    trailer: Flight,
    interaction: pd.DataFrame,
    transformer: callable = None,
) -> go.Figure:
    """
    Visualise the interaction between the vortex of a leader aircraft and the trailer on
    the xy plane.

    Parameters
    ----------
    leader : Flight
        Flight object of the leader aircraft
    trailer : Flight
        Flight object of the trailer aircraft
    interaction : pd.DataFrame
        DataFrame containing the interaction information, as returned by
        `find_interactions`
    transformer : callable, optional
        Function to transform the x and y coordinates of the interaction to lat and lon

    Returns
    -------
    go.Figure
        Plotly figure containing the interaction visualisation on the xy plane
    """

    if transformer is None:
        transformer = coords.xy_to_latlon

    df_trailer_reduced = trailer.data[
        trailer.data.timestamp <= interaction.loc["timestamp"][0]
    ]

    # Define X-Y view figure
    fig_xy = go.Figure()

    # Add leader and trailer trajectories
    fig_xy.add_trace(
        go.Scattermapbox(
            mode="lines",
            lon=leader.data.longitude,
            lat=leader.data.latitude,
            line=dict(color="cyan", width=1),
            name="Trajectory leader",
            legendgroup="Trajectories",
        )
    )
    fig_xy.add_trace(
        go.Scattermapbox(
            mode="lines",
            lon=df_trailer_reduced.longitude,
            lat=df_trailer_reduced.latitude,
            line=dict(color="magenta", width=1),
            name="Trajectory trailer",
            legendgroup="Trajectories",
        )
    )

    # Add vortex uncertainties
    for uncertainty, color in zip(["3s", "2s"], ["red", "orange"]):
        lon, lat = transformer(
            [
                interaction.loc[f"x_{uncertainty}_l_ref"],
                interaction.loc[f"x_{uncertainty}_r_ref"],
            ],
            [
                interaction.loc[f"y_{uncertainty}_l_ref"],
                interaction.loc[f"y_{uncertainty}_r_ref"],
            ],
        )

        fig_xy.add_trace(
            go.Scattermapbox(
                mode="lines + markers",
                lon=lon,
                lat=lat,
                marker=dict(color=color),
                name=f"Interaction -> Cores uncertainty {uncertainty[0]}σ",
                legendgroup="Interaction",
            )
        )

    # Plot aircraft icon
    geoms = marker_mapbox(
        df_trailer_reduced.iloc[-1],
        size=0.00005,
        lat="latitude",
        lon="longitude",
    )
    x, y = geoms.exterior.coords.xy
    fig_xy.add_trace(
        go.Scattermapbox(
            fill="toself",
            lon=list(x),
            lat=list(y),
            fillcolor="magenta",
            name="Interaction -> Trailer position",
            legendgroup="Interaction",
            marker=dict(size=0, color="magenta"),
        ),
    )

    # Plot wake cores
    lon_l, lat_l = transformer(
        interaction.loc["xl_ref"], interaction.loc["yl_ref"]
    )
    lon_r, lat_r = transformer(
        interaction.loc["xr_ref"], interaction.loc["yr_ref"]
    )
    fig_xy.add_trace(
        go.Scattermapbox(
            mode="markers",
            lon=lon_l,
            lat=lat_l,
            marker=dict(color="royalblue"),
            name="Interaction -> Left wake core",
            legendgroup="Interaction",
        )
    )
    fig_xy.add_trace(
        go.Scattermapbox(
            mode="markers",
            lon=lon_r,
            lat=lat_r,
            marker=dict(color="green"),
            name="Interaction -> Right wake core",
            legendgroup="Interaction",
        )
    )

    # Adjust figure layout
    fig_xy.update_layout(
        width=1000,
        margin={"l": 0, "b": 0, "t": 40, "r": 50},
        mapbox_zoom=14,
        mapbox_center_lat=df_trailer_reduced.latitude.iloc[-1],
        mapbox_center_lon=df_trailer_reduced.longitude.iloc[-1],
        showlegend=True,
        mapbox_style="carto-positron",
        title_text="Interaction horizontal view",
    )

    # Add annotations (circulation, wind direction and speed)
    fig_xy.add_annotation(
        x=1.02,
        y=0.55,
        text=(
            "Circulation left vortex: "
            f"{round(interaction.loc['gam_l'][0], 3)} m\u00b2/s"
        ),
        xref="paper",
        yref="paper",
        xanchor="left",
        showarrow=False,
        valign="middle",
    )
    fig_xy.add_annotation(
        x=1.02,
        y=0.50,
        text=(
            "Circulation right vortex: "
            f"{round(interaction.loc['gam_r'][0], 3)} m\u00b2/s"
        ),
        xref="paper",
        yref="paper",
        xanchor="left",
        showarrow=False,
        valign="middle",
    )
    fig_xy.add_annotation(
        x=1.02,
        y=0.40,
        text=(
            "Wind direction: " f"{round(interaction.loc['wind_dir'][0], 3)} °"
        ),
        xref="paper",
        yref="paper",
        xanchor="left",
        showarrow=False,
        valign="middle",
    )
    fig_xy.add_annotation(
        x=1.02,
        y=0.35,
        text=(
            "Wind speed: " f"{round(interaction.loc['wind_speed'][0], 3)} kts"
        ),
        xref="paper",
        yref="paper",
        xanchor="left",
        showarrow=False,
        valign="middle",
    )
    ax, ay = wind_direction_arrow(interaction.loc["wind_dir"][0])
    fig_xy.add_annotation(
        xref="paper",
        yref="paper",
        x=1.1,
        y=0.3,
        showarrow=True,
        arrowhead=1,
        arrowcolor="black",
        arrowwidth=2,
        ax=ax,
        ay=ay,
        axref="pixel",
        ayref="pixel",
        arrowside="end",
    )

    # Return figure
    return fig_xy


def visualise_interaction_yz(
    leader: Flight, trailer: Flight, interaction: pd.DataFrame
) -> go.Figure:
    """
    Visualise the interaction between the vortex of a leader aircraft and the trailer on
    the yz plane.

    Parameters
    ----------
    leader : Flight
        Flight object of the leader aircraft
    trailer : Flight
        Flight object of the trailer aircraft
    interaction : pd.DataFrame
        DataFrame containing the interaction information, as returned by
        `find_interactions`

    Returns
    -------
    go.Figure
        Plotly figure containing the interaction visualisation on the yz plane
    """

    # Y-Z view -------------------------------------------------------------------------

    # get airport elevation
    ap_elevation = leader.landing_airport().altitude

    # Plot Environment (runway and ground)
    fig_yz = px.scatter()
    fig_yz.add_trace(
        go.Scatter(
            x=[0, 0],
            y=[
                ap_elevation,
                max(
                    float(interaction.loc["z_3s_hi_ref"].iloc[0]),
                    float(interaction.loc["geoaltitude"].iloc[0]),
                ),
            ],
            mode="lines",
            line=dict(color="black"),
            name="Runway centerline axis",
            legendgroup="Environment",
        )
    )
    fig_yz.add_trace(
        go.Scatter(
            x=[-500, 500],
            y=[ap_elevation, ap_elevation],
            mode="lines",
            line=dict(color="mediumseagreen"),
            name="Ground",
            legendgroup="Environment",
        )
    )

    # Plot interaction (trailer, wake cores, 2s uncertainty and 3s uncertainty)

    # y_correction = trailer.at(interaction.loc["timestamp"][0]).b_diff
    y_correction = trailer.data[
        trailer.data.timestamp <= interaction.loc["timestamp"][0]
    ]["b_diff"].iloc[0]

    fig_yz.add_trace(
        go.Scatter(
            x=[float(interaction.loc["rwy_offset_m"].iloc[0])],
            y=[float(interaction.loc["geoaltitude"].iloc[0])],
            mode="markers",
            line=dict(color="magenta"),
            name="Trailer position",
            legendgroup="Interaction -> Trailer position",
        )
    )
    fig_yz.add_trace(
        go.Scatter(
            x=[float(interaction.loc["yr"].iloc[0]) + y_correction],
            y=[float(interaction.loc["zr_ref"].iloc[0])],
            mode="markers",
            line=dict(color="green"),
            name="Interaction -> Right wake core",
            legendgroup="Interaction",
        )
    )
    fig_yz.add_trace(
        go.Scatter(
            x=[float(interaction.loc["yl"].iloc[0]) + y_correction],
            y=[float(interaction.loc["zl_ref"].iloc[0])],
            mode="markers",
            line=dict(color="royalblue"),
            name="Interaction -> Left wake core",
            legendgroup="Interaction",
        )
    )
    x_vals_ellipse_2s, y_vals_ellipse_2s = get_ellipse_points(
        float(interaction.loc["y_2s_l"].iloc[0] + y_correction),
        float(interaction.loc["y_2s_r"].iloc[0] + y_correction),
        float(interaction.loc["z_2s_lo_ref"].iloc[0]),
        float(interaction.loc["z_2s_hi_ref"].iloc[0]),
        num_points=100,
    )
    fig_yz.add_trace(
        go.Scatter(
            x=x_vals_ellipse_2s,
            y=y_vals_ellipse_2s,
            mode="lines",
            line=dict(color="orange"),
            name="Interaction -> Cores uncertainty 2σ",
            legendgroup="Interaction",
        )
    )
    x_vals_ellipse_3s, y_vals_ellipse_3s = get_ellipse_points(
        float(interaction.loc["y_3s_l"].iloc[0] + y_correction),
        float(interaction.loc["y_3s_r"].iloc[0] + y_correction),
        float(interaction.loc["z_3s_lo_ref"].iloc[0]),
        float(interaction.loc["z_3s_hi_ref"].iloc[0]),
        num_points=100,
    )
    fig_yz.add_trace(
        go.Scatter(
            x=x_vals_ellipse_3s,
            y=y_vals_ellipse_3s,
            mode="lines",
            line=dict(color="red"),
            name="Interaction -> Cores uncertainty 3σ",
            legendgroup="Interaction",
        )
    )

    # Add circulation annotations
    fig_yz.add_annotation(
        x=1.02,
        y=0.55,
        text=(
            "Circulation left vortex: "
            f"{round(interaction.loc['gam_l'][0], 3)} m\u00b2/s"
        ),
        xref="paper",
        yref="paper",
        xanchor="left",
        showarrow=False,
        valign="middle",
    )
    fig_yz.add_annotation(
        x=1.02,
        y=0.50,
        text=(
            "Circulation right vortex: "
            f"{round(interaction.loc['gam_r'][0], 3)} m\u00b2/s"
        ),
        xref="paper",
        yref="paper",
        xanchor="left",
        showarrow=False,
        valign="middle",
    )

    # Update figure layout
    fig_yz.update_layout(
        width=1000,
        margin=dict(l=40, r=50, t=40, b=20),
        title_text="Interaction Y-Z view",
    )

    fig_yz.update_xaxes(title_text="Distance from centerline [m]")
    fig_yz.update_yaxes(title_text="Geoaltitude [ft]")

    # Return figure
    return fig_yz


def visualise_interaction_xz(
    leader: Flight, trailer: Flight, interaction: pd.DataFrame
) -> go.Figure:
    """
    Visualise the interaction between the vortex of a leader aircraft and the trailer on
    the xz plane.

    Parameters
    ----------
    leader : Flight
        Flight object of the leader aircraft
    trailer : Flight
        Flight object of the trailer aircraft
    interaction : pd.DataFrame
        DataFrame containing the interaction information, as returned by
        `find_interactions`

    Returns
    -------
    go.Figure
        Plotly figure containing the interaction visualisation on the xz plane
    """

    # get airport elevation
    ap_elevation = leader.landing_airport().altitude

    # get trailer data closest to interaction timestamp
    df_trailer_reduced = trailer.data[
        trailer.data.timestamp <= interaction.loc["timestamp"][0]
    ]

    # Plot environment (runway and ground)
    fig_xz = go.Figure()
    max_dist = min(
        leader.data.distance.max(), df_trailer_reduced.distance.max(), 10
    )
    fig_xz.add_trace(
        go.Scatter(
            x=[0, max_dist],
            y=[ap_elevation, ap_elevation],
            mode="lines",
            line=dict(color="mediumseagreen"),
            name="Ground",
            legendgroup="Environment",
        )
    )
    fig_xz.add_trace(
        go.Scatter(
            x=[0, -1],
            y=[ap_elevation, ap_elevation],
            mode="lines",
            line=dict(color="dimgray"),
            name="Runway",
            legendgroup="Environment",
        )
    )
    fig_xz.add_trace(
        go.Scatter(
            x=[0],
            y=[ap_elevation],
            mode="markers",
            marker_symbol="arrow-up",
            marker=dict(color="black"),
            name="Runway threshold",
            legendgroup="Environment",
        )
    )

    # Plot leader and trailer trajectories
    fig_xz.add_trace(
        go.Scatter(
            x=leader.data.distance,
            y=leader.data.geoaltitude,
            mode="lines",
            line=dict(color="cyan"),
            name="Trajectory leader",
            legendgroup="Trajectories",
            opacity=0.4,
        )
    )
    fig_xz.add_trace(
        go.Scatter(
            x=df_trailer_reduced.distance,
            y=df_trailer_reduced.geoaltitude,
            mode="lines",
            line=dict(color="magenta"),
            name="Trajectory trailer",
            legendgroup="Trajectories",
            opacity=0.4,
        )
    )

    # Plot interaction (trailer position, vortex cores, 2s uncertainties, 3s
    # uncertainties)
    fig_xz.add_trace(
        go.Scatter(
            x=[interaction.loc["x_rway"][0], interaction.loc["x_rway"][0]],
            y=[
                interaction.loc["z_3s_hi_ref"][0],
                interaction.loc["z_3s_lo_ref"][0],
            ],
            mode="lines + markers",
            line=dict(color="red"),
            legendgroup="Interaction",
            name="Interaction -> cores uncertainty 3σ",
        )
    )
    fig_xz.add_trace(
        go.Scatter(
            x=[interaction.loc["x_rway"][0], interaction.loc["x_rway"][0]],
            y=[
                interaction.loc["z_2s_hi_ref"][0],
                interaction.loc["z_2s_lo_ref"][0],
            ],
            mode="lines + markers",
            line=dict(color="orange"),
            legendgroup="Interaction",
            name="Interaction -> cores uncertainty 2σ",
        )
    )
    fig_xz.add_trace(
        go.Scatter(
            x=[interaction.loc["distance"][0]],
            y=[interaction.loc["geoaltitude"][0]],
            mode="markers",
            marker=dict(color="magenta"),
            name="Interaction -> trailer position",
            legendgroup="Interaction",
        )
    )
    fig_xz.add_trace(
        go.Scatter(
            x=interaction.loc["x_rway"],
            y=interaction.loc["zl_ref"],
            mode="markers",
            marker=dict(color="royalblue"),
            legendgroup="Interaction",
            name="Interaction -> left wake core",
        )
    )
    fig_xz.add_trace(
        go.Scatter(
            x=interaction.loc["x_rway"],
            y=interaction.loc["zr_ref"],
            mode="markers",
            marker=dict(color="green"),
            legendgroup="Interaction",
            name="Interaction -> right wake core",
        )
    )

    # Add wind direction and speed annotations
    fig_xz.add_annotation(
        x=1.02,
        y=0.34,
        text=(
            "Circulation left vortex: "
            f"{round(interaction.loc['gam_l'][0], 3)} m\u00b2/s"
        ),
        xref="paper",
        yref="paper",
        xanchor="left",
        showarrow=False,
        valign="middle",
    )

    # Add annotations (circulation, wind direction and speed)
    fig_xz.add_annotation(
        x=1.02,
        y=0.27,
        text=(
            "Circulation right vortex: "
            f"{round(interaction.loc['gam_r'][0], 3)} m\u00b2/s"
        ),
        xref="paper",
        yref="paper",
        xanchor="left",
        showarrow=False,
        valign="middle",
    )

    # Update figure layout
    fig_xz.update_layout(
        width=1000,
        margin=dict(l=40, r=50, t=40, b=20),
        title_text="Interaction X-Z view",
    )
    fig_xz.update_xaxes(
        title_text="Distance from threshold [nm]",
        # autorange="reversed",
        range=[max_dist, -1],
    )
    fig_xz.update_yaxes(title_text="Geoaltitude [ft]")

    # Return figure
    return fig_xz
