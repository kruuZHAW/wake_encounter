import math
import os
import subprocess
from datetime import datetime
from os.path import join as opj
from typing import Optional, Tuple, Union
from warnings import warn

import pandas as pd
from traffic.core import Flight

from .aircraft import Aircraft
from .meteo import EDR, Meteo
from utils import coords, viz


timelike = Union[str, datetime, pd.Timestamp]

KNTS2MPS = 0.5144444444444445
FT2M = 0.3048
NM2M = 1852


def generate_p2p_wake(
    aircraft: Aircraft,
    meteo: Meteo,
    edr: EDR,
    run_id: int = 1,
    path: str = None,
    verbose: bool = False,
) -> pd.DataFrame:
    """Generate wake using the P2P model.

    Parameters
    ----------
    aircraft : Aircraft
        Generator aircraft.
    meteo : Meteo
        Meteo data.
    edr : EDR
        Eddy dissipation rate data.
    run_id : int, optional
        Run ID, by default 1
    path : str, optional
        If set, an ABSOLUTE path to the data files. If not set the default directory
        is used (in this package, folder "p2p"), by default None
    verbose : bool, optional
        Verbose output from P2P, by default False

    Returns
    -------
    DataFrame
        Results from P2P.
    """

    if path is None:
        path_to_dir = "p2p"
    else:
        path_to_dir = path

    requred_dirs = [
        opj(path_to_dir, "inputs"),
        opj(path_to_dir, "results"),
        opj(path_to_dir, "inputs", f"{run_id}"),
        opj(path_to_dir, "results", f"{run_id}"),
    ]

    for req_dir in requred_dirs:
        if not os.path.exists(req_dir):
            os.makedirs(req_dir)

    aircraft.write_aircraft_dat_file(
        file_path=opj(path_to_dir, "inputs", f"{run_id}", "ac_init.dat")
    )
    meteo.write_meteo_dat_file(
        file_path=opj(path_to_dir, "inputs", f"{run_id}", "meteo.dat")
    )
    edr.write_edr_dat_file(
        file_path=opj(path_to_dir, "inputs", f"{run_id}", "edr.dat")
    )

    ret = run_p2p(run_id, path=path, verbose=verbose)
    if verbose:
        print(ret.stdout.decode("utf-8"))
    if ret.returncode != 0:
        raise Exception(
            f"P2P did not run successfully. Error:\n{ret.stderr.decode('utf-8')}"
        )

    res = read_results(run_id, opj(path_to_dir, "results"))
    return res


class Wake:
    """A class representing the wake of an aircraft.

    Attributes
    ----------
    df : pandas.DataFrame
        The wake dataframe.

    Methods
    -------
    generate(cls, aircraft, meteo, edr, at_time) -> Wake:
        Generates a Wake object from the aircraft, meteo and edr data.
    """

    def __init__(self, df):
        self.df = df

    def __repr__(self):
        return self.df.__repr__()

    def _repr_html_(self):
        return self.df._repr_html_()

    @classmethod
    def generate(
        cls,
        aircraft: Union[Flight, Aircraft, str],
        meteo: Union[Meteo, str],
        edr: Union[EDR, str],
        at_time: Optional["timelike"] = None,
        run_id: int = 1,
        path: str = None,
        alt_gnd: float = 0,
        verbose: bool = False,
    ) -> "Wake":
        """
        Generates a Wake object from the aircraft, meteo and edr data.

        Parameters
        ----------
        aircraft : Union[Flight, Aircraft, str]
            Flight object, Aircraft object or path to the aircraft data file.
        meteo : Union[Meteo, str]
            Meteo object or path to the meteo data file.
        edr : Union[EDR, str]
            EDR object or path to the EDR data file.
        at_time : Optional["timelike"], optional
            Time at which to generate the wake, by default None
        run_id : int, optional
            Run ID, by default 1
        path : str, optional
            If set, an ABSOLUTE path to the data files. If not set the default directory
            is used (in this package, folder "p2p"), by default None
        verbose : bool, optional
            Verbose output from P2P, by default False

        Returns
        -------
        Wake
            A Wake object.
        """
        if isinstance(aircraft, Flight):
            flight = aircraft
            aircraft = Aircraft.from_flight(aircraft, at_time=at_time)
        elif isinstance(aircraft, str):
            aircraft = Aircraft.from_dat_file(aircraft)

        if isinstance(meteo, str):
            meteo = Meteo.from_dat_file(meteo)
        if isinstance(edr, str):
            edr = EDR.from_dat_file(edr)

        df = generate_p2p_wake(
            aircraft, meteo, edr, run_id=run_id, path=path, verbose=verbose
        ).dropna()
        if ("flight" in locals()) and isinstance(flight, Flight):
            # if it's a flight, get relative position
            try:
                alt_gnd = flight.landing_airport().altitude
                fdata = flight.at(at_time)
                assert fdata is not None
                df["zl_ref"] = df["zl"] / FT2M + alt_gnd
                df["zr_ref"] = df["zr"] / FT2M + alt_gnd
                df["z_2s_lo_ref"] = df["z_2s_lo"] / FT2M + alt_gnd
                df["z_2s_hi_ref"] = df["z_2s_hi"] / FT2M + alt_gnd
                df["z_3s_lo_ref"] = df["z_3s_lo"] / FT2M + alt_gnd
                df["z_3s_hi_ref"] = df["z_3s_hi"] / FT2M + alt_gnd
                df["x_rway"] = fdata.distance - df["x"] / NM2M
                angle = 90 - fdata["track"]
                if angle < 0:
                    angle += 360
                df = df.apply(
                    lambda wrow: coords.change_ref_xy(
                        wrow, fdata.x, fdata.y, angle
                    ),
                    1,
                    result_type="expand",
                )
                df["wind_dir"] = meteo.get_wind_dir(10)
                df["wind_speed"] = meteo.get_wind_speed(10)
                df.index = pd.to_timedelta(df.index, unit="s")
                df.index = df.index + fdata.name
            except NameError as e:
                print(f"{e}")
                pass
        else:
            # # assume that the coordinates are already relative
            # df["zl_ref"] = df["zl"] / FT2M + alt_gnd
            # df["zr_ref"] = df["zr"] / FT2M + alt_gnd
            # df["z_2s_lo_ref"] = df["z_2s_lo"] / FT2M + alt_gnd
            # df["z_2s_hi_ref"] = df["z_2s_hi"] / FT2M + alt_gnd
            # df["z_3s_lo_ref"] = df["z_3s_lo"] / FT2M + alt_gnd
            # df["z_3s_hi_ref"] = df["z_3s_hi"] / FT2M + alt_gnd
            # df["x_rway"] = df["x"] / NM2M
            # angle = (90 - 280) * 0
            # x_origin = 500
            # y_origin = 0
            # df = df.apply(
            #     lambda wrow: coords.change_ref_xy(
            #         wrow, x_origin, y_origin, angle
            #     ),
            #     1,
            #     result_type="expand",
            # )
            # df["wind_dir"] = meteo.dir[0]
            # df.index = pd.to_timedelta(df.index, unit="s")
            # if at_time is not None:
            #     df.index = df.index + at_time
            pass

        return cls(df)

    def interaction_with_trailer(self, trailer_flight: Flight) -> pd.DataFrame:
        """
        Determines the timestamp of interaction (closest point) between the wake and
        the trailer aircraft. Information about the interaction is returned in a
        dataframe.

        Parameters
        ----------
        trailer_flight : Flight
            Flight object of the trailer aircraft.

        Returns
        -------
        pd.DataFrame
            A dataframe containing information about the interaction between the wake
            and trailer aircraft.
        """
        d_wake = self.df.x_rway.mean()
        tdata = trailer_flight.data
        tdata["rwy_offset_m"] = tdata.b_diff * NM2M
        # return row with closest value to d_wake
        trailer_row = tdata.loc[tdata["distance"].sub(d_wake).abs().idxmin()]
        wrs = (
            self.df.resample("0.001S")
            .interpolate()
            .resample("1S")
            .interpolate()
        )
        row = wrs.loc[trailer_row.timestamp]
        return pd.DataFrame(pd.concat([row, trailer_row]))

    def plot(self):
        """
        Visualises the wake of the aircraft on the y-z plane.
        """
        return viz.visualise_vortex(self.df)


def run_p2p(
    run_id: int, path: str = None, verbose: bool = False
) -> subprocess.CompletedProcess:
    # TODO:
    # I think the whole assembly of the fort.13 file could be mainly skipped since, as
    # far as I can tell, the options are not used in the P2P code. (mora, 07.09.2023)

    from sys import platform

    if not platform.startswith("linux"):
        # different OS than Linux were not tested, if it works, remove this
        raise Exception("P2P is intended to run on Linux.")

    path_lidar = "Lidar_110407_054603UTC_A320.dat"
    first = "1"
    campaign = "ZHAW"
    rwy = "14"
    path_edr = opj("inputs", f"{run_id}", "EDR.dat")
    path_ac = opj("inputs", f"{run_id}", "ac_init.dat")
    # path_meteo = opj("inputs", f"{run_id}", "meteo.dat")
    # path_lidar contains a bunch of information on the case, unpack it
    date = path_lidar[6:12]
    ac = path_lidar[23:27]
    h = path_lidar[13:15]
    min = path_lidar[15:17]
    s = path_lidar[17:19]

    # assemble file "fort.13"
    if path is None:
        # use the default, relative path
        fort13_path = opj("p2p", "inputs", f"{run_id}", "fort.13")
    else:
        fort13_path = opj(path, "inputs", f"{run_id}", "fort.13")

    with open(fort13_path, "w") as fort13:
        fort13.write(path_lidar + "\n")
        fort13.write(first + "\n")
        fort13.write(ac + "\n")  # A/C
        fort13.write(date + "\n")  # d
        fort13.write(h + "\n")  # h
        fort13.write(min + "\n")  # min
        fort13.write(s + "\n")  # s
        fort13.write(campaign + "\n")  # camp
        fort13.write(rwy + "\n")  # rwy
        fort13.write("" + "\n")  # uac
        fort13.write(path_edr + "\n")  # edrfile
        fort13.write(path_ac + "\n")  # acfile

    # assemble command line call
    root_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    p2p_path = opj(root_dir, "p2p")

    cmd_str = f"cd {p2p_path}; ./P2P -o {run_id}"
    if path is not None:
        cmd_str += f" -p {path}"
    if verbose:
        cmd_str += " -v"

    # run P2P and return the command line output
    ret = subprocess.run(cmd_str, shell=True, capture_output=True)
    return ret


def read_results(run_id, pathres: str = opj("p2p", "results")) -> pd.DataFrame:
    pathres = opj(pathres, f"{run_id}")
    df_traj = pd.read_csv(
        opj(pathres, "TRAJEC.dat"),
        sep=r"\s+",  # whitespace
        header=None,
        names=["t", "gam_l", "yl", "zl", "gam_r", "yr", "zr", "x"],
    ).set_index("t")
    df_traj["gam_l"] = df_traj["gam_l"].abs()
    df_traj["gam_r"] = df_traj["gam_r"].abs()
    df_z = pd.read_csv(
        opj(pathres, "P2P_lev_z.dat"),
        sep=r"\s+",
        header=None,
        names=["t", "z_2s_lo", "z_2s_hi", "z_3s_lo", "z_3s_hi"],
    ).set_index("t")
    df_y = pd.read_csv(
        opj(pathres, "P2P_lev_y.dat"),
        sep=r"\s+",
        header=None,
        names=["t", "y_2s_l", "y_2s_r", "y_3s_l", "y_3s_r"],
    ).set_index("t")

    df_g = pd.read_csv(
        opj(pathres, "P2P_lev_g.dat"),
        sep=r"\s+",
        header=None,
        names=["t", "g_2s_lo", "g_2s_hi", "g_3s_lo", "g_3s_hi"],
    ).set_index("t")
    return pd.concat([df_traj, df_z, df_y, df_g], axis=1)


def rotate_and_translate(
    x: float, y: float, angle_deg: float, dx: float, dy: float
) -> Tuple[float, float]:
    # Convert angle to radians
    angle_rad = math.radians(angle_deg)

    # Apply rotation
    x_rot = x * math.cos(angle_rad) - y * math.sin(angle_rad)
    y_rot = x * math.sin(angle_rad) + y * math.cos(angle_rad)

    # Apply translation
    x_trans = x_rot + dx
    y_trans = y_rot + dy
    return x_trans, y_trans
