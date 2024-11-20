import os
from os.path import join as opj
import pandas as pd
from wake.aircraft import Aircraft
from wake.meteo import EDR, Meteo


def run_p2p(path="p2p") -> None:
    """
    Runs the p2p model.

    Parameters
    ----------
    path : str, optional
        Path to the p2p model, by default "p2p"
    """
    path_lidar = "Lidar_110407_054603UTC_A320.dat"
    first = "1"
    campaign = "ZHAW"
    rwy = "14"
    path_edr = "EDR.dat"
    path_ac = "ac_init.dat"

    # path_lidar contains a bunch of information on the case, unpack it
    date = path_lidar[6:12]
    ac = path_lidar[23:27]
    h = path_lidar[13:15]
    min = path_lidar[15:17]
    s = path_lidar[17:19]

    # assemble file "fort.13"
    with open(opj(path, "fort.13"), "w") as fort13:
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

    # run program
    os.system("cd p2p; ./P2P")


def read_results(pathres: str = "p2p") -> pd.DataFrame:
    """
    Reads the results of the p2p model and returns them in an aggregated form as a
    pandas dataframe.

    Parameters
    ----------
    pathres : str, optional
        Path to the folder containing the output files of the p2p model (TRAJEC.DAT,
        P2P_lev_z.dat, P2P_lev_y.dat, P2P_lev_g.dat), by default "p2p"

    Returns
    -------
    pd.DataFrame
        A dataframe containing the results of the p2p model in an aggregated form
    """
    # Process file containing deterministic values
    df_traj = pd.read_csv(
        opj(pathres, "TRAJEC.DAT"),
        sep=r"\s+",  # whitespace
        header=None,
        names=["t", "gam_l", "yl", "zl", "gam_r", "yr", "zr", "x"],
    ).set_index("t")
    df_traj["gam_l"] = df_traj["gam_l"].abs()
    df_traj["gam_r"] = df_traj["gam_r"].abs()
    # Process files containing probabilistic envelopes for the vortex heights
    df_z = pd.read_csv(
        opj(pathres, "P2P_lev_z.dat"),
        sep=r"\s+",
        header=None,
        names=["t", "z_2s_lo", "z_2s_hi", "z_3s_lo", "z_3s_hi"],
    ).set_index("t")
    # Process files containing probabilistic envelopes for the vortex lateral positions
    df_y = pd.read_csv(
        opj(pathres, "P2P_lev_y.dat"),
        sep=r"\s+",
        header=None,
        names=["t", "y_2s_l", "y_2s_r", "y_3s_l", "y_3s_r"],
    ).set_index("t")
    # Process files containing probabilistic envelopes for the vortex circulation
    df_g = pd.read_csv(
        opj(pathres, "P2P_lev_g.dat"),
        sep=r"\s+",
        header=None,
        names=["t", "g_2s_lo", "g_2s_hi", "g_3s_lo", "g_3s_hi"],
    ).set_index("t")
    # return combined dataframe
    return pd.concat([df_traj, df_z, df_y, df_g], axis=1)


def generate_p2p_wake(
    aircraft: Aircraft, meteo: Meteo, edr: EDR
) -> pd.DataFrame:
    """
    For a given aircraft, meteo and edr, generate the wake using the p2p model. The
    results are returned as a dataframe.

    Parameters
    ----------
    aircraft : Aircraft
        Aircraft object corresponding to the wake generator
    meteo : Meteo
        Meteo object matching the conditions at the time of the wake generation
    edr : EDR
        EDR object matching the conditions at the time of the wake generation

    Returns
    -------
    pd.DataFrame
        A dataframe containing the results of the p2p model in an aggregated form
    """
    # write files
    aircraft.write_aircraft_dat_file()
    meteo.write_meteo_dat_file()
    edr.write_edr_dat_file()
    # run p2p model
    run_p2p()
    # read results and return
    return read_results()
