from dataclasses import dataclass
from datetime import datetime
from typing import Iterable, Optional, Union
import os
import numpy as np
import pandas as pd
import scipy.stats as st
from traffic.data.weather.metar import METAR

timelike = Union[str, datetime, pd.Timestamp]

KNTS2MPS = 0.5144444444444445
FT2M = 0.3048


@dataclass
class Meteo:
    """
    A class representing meteorological data.

    Attributes:
    -----------
    zm : Iterable[float]
        Iterable of floats representing the height in meters.
    um : Iterable[float]
        Iterable of floats representing the eastward wind component in m/s.
    vm : Iterable[float]
        Iterable of floats representing the northward wind component in m/s.
    wm : Iterable[float]
        Iterable of floats representing the vertical wind component in m/s.
    qq : Iterable[float]
        Iterable of floats representing the specific humidity in kg/kg.
    temp : Iterable[float]
        Iterable of floats representing the temperature in Kelvin.
    th : Iterable[float]
        Iterable of floats representing the potential temperature in Kelvin.
    tke : Iterable[float]
        Iterable of floats representing the turbulence kinetic energy in m^2/s^2.

    Methods:
    --------
    from_dat_file(filepath):
        Reads meteo data from a dat file and returns an instance of the Meteo
        class.
    write_meteo_dat_file(file_path):
        Writes meteorological data to a dat file.


    """

    zm: Iterable[float]  # height (m)
    um: Iterable[float]  # eastward wind component (m/s)
    vm: Iterable[float]  # northward wind component (m/s)
    wm: Iterable[float]  # vertical wind component (m/s)
    qq: Iterable[float]  # specific humidity (kg/kg)
    temp: Iterable[float]  # temperature (K)
    dir: Iterable[float]  # wind direction (deg)
    th: Iterable[float]  # potential temperature (K)
    tke: Iterable[float]  # turbulence kinetic energy (m^2/s^2)

    def __repr__(self):
        return pd.DataFrame(self.__dict__).__repr__()

    def _repr_html_(self):
        return pd.DataFrame(self.__dict__)._repr_html_()

    @classmethod
    def from_idaweb_klo(
        cls,
        timestamp: timelike,
        bearing: float,
        qq: float = 0.05,
        idaweb_file_path: Optional[str] = None,
        tke: Optional[float] = None,
        extrapolate: Optional[bool] = True,
    ):
        """
        Reads meteorological data from idaweb and returns an instance of the
        Meteo class.

        Parameters:
        -----------
        cls: Meteo
            Class name.
        timestamp : timelike
            Timestamp for which meteo data is required.
        bearing : float
            Bearing of the flight.
        qq : float, optional (default=0.05)
            Mean turbulence velocity.
        idaweb_file_path : str, optional (default=None)
            Path to the idaweb file. This can either be a CSV or a pickle and the file
            type will be determined based on the extension. If None, the default path is
            used.
        tke : float, optional (default=None)
            Turbulent kinetic energy.
        extrapolate : bool, optional (default=True)
            If True, the wind speed and direction are extrapolated to different
            altitudes.

        Returns:
        --------
        Meteo
            An instance of the Meteo class.
        """

        if idaweb_file_path is None:
            idaweb_file_path = "data/meteo_KLO.pkl"

        if idaweb_file_path.endswith(".csv"):
            # types = defaultdict(lambda: float, time=str, stn=str)
            # df = pd.read_csv(idaweb_file_path, sep=",", dtype=types)
            df = pd.read_csv(idaweb_file_path, sep=",", dtype=str)
            df["time"] = pd.to_datetime(
                df["time"], format="%Y-%m-%d %H:%M:%S%z", utc=True
            )
            no_conv_cols = ["time", "stn"]
            for col in df.columns:
                if col not in no_conv_cols:
                    df[col] = pd.to_numeric(df[col], errors="coerce")

            df.set_index("time", inplace=True)

        elif idaweb_file_path.endswith(".pkl") or idaweb_file_path.endswith(
            ".pickle"
        ):
            df = pd.read_pickle(idaweb_file_path)
        else:
            raise ValueError(
                "idaweb_file_path should be either a csv or a pickle file"
            )
        label_wind_dir = "dkl010z0"  # deg
        label_wind_vel = "fkl010z0"  # m/s
        label_pot_temp = "tpo200s0"  # deg C
        label_temp = "tre005s0"  # deg C
        label_wind_vel_std = "fkl010za"
        label_wind_dir_std = "dkl010za"
        meteo_time = df.iloc[
            df.index.get_indexer([timestamp], method="nearest")[0]
        ]
        wind_vel = meteo_time[label_wind_vel]
        wind_dir = meteo_time[label_wind_dir]
        wind_vel_std = float(meteo_time[label_wind_vel_std])
        wind_dir_std = float(meteo_time[label_wind_dir_std])
        print(
            f"wind speed: {wind_vel} m/s std:{wind_vel_std}",
            f", Wind direction: {wind_dir} deg std: {wind_dir_std} deg",
        )
        pot_temp = meteo_time[label_pot_temp] + 273.15
        temp = meteo_time[label_temp] + 273.15
        alpha = np.deg2rad(wind_dir - bearing)
        crosswind = -wind_vel * np.sin(alpha)
        headwind = wind_vel * np.cos(alpha)
        if tke is None:
            tke = compute_tke(
                wind_vel,
                wind_vel_std,
                wind_dir,
                wind_dir_std,
                n=10000,
            )
        if extrapolate:
            headwinds = []
            crosswinds = []
            tkes = []
            alts = [2, 5, 10, 20, 50, 100, 200, 300, 400, 500]
            n = len(alts)
            for alt in alts:
                headwinds.append(extrapolate_wind_at_alt(alt, 10, headwind))
                crosswinds.append(extrapolate_wind_at_alt(alt, 10, crosswind))
                tkes.append(
                    compute_tke(
                        extrapolate_wind_at_alt(alt, 10, wind_vel),
                        wind_vel_std,
                        wind_dir,
                        wind_dir_std,
                        n=10000,
                    )
                )
            return cls(
                alts,
                headwinds,
                crosswinds,
                [0] * n,
                [qq] * n,
                [temp] * n,
                [wind_dir] * n,
                [pot_temp] * n,
                tkes,
            )

        return cls(
            [10],
            [headwind],
            [crosswind],
            [0],
            [qq],
            [temp],
            [wind_dir],
            [pot_temp],
            [tke],
        )

    @classmethod
    def from_metar(
        cls,
        airport: str,
        timestamp: timelike,
        bearing: float,
        pot_temp: float = 300,
        qq: float = 0.05,
        tke: float = 0.01,
        extrapolate=True,
    ):
        df = METAR(airport).get(timestamp).set_index("time")
        row = df.iloc[df.index.get_indexer([timestamp], method="nearest")[0]]
        wind_vel = row.wind_speed.value() * KNTS2MPS
        wind_dir = row.wind_dir.value()
        temp = row.temp.value() + 273.15
        alpha = np.deg2rad(wind_dir - bearing)
        crosswind = -wind_vel * np.sin(alpha)
        headwind = wind_vel * np.cos(alpha)
        if extrapolate:
            headwinds = []
            crosswinds = []
            alts = [2, 5, 10, 20, 50, 100, 200, 300, 400, 500]
            n = len(alts)
            for alt in alts:
                headwinds.append(extrapolate_wind_at_alt(alt, 10, headwind))
                crosswinds.append(extrapolate_wind_at_alt(alt, 10, crosswind))

            return cls(
                alts,
                headwinds,
                crosswinds,
                [0] * n,
                [qq] * n,
                [temp] * n,
                [wind_dir] * n,
                [pot_temp] * n,
                [tke] * n,
            )

        return cls(
            [10], [headwind], [crosswind], [0], [qq], [temp], [pot_temp], [tke]
        )

    @classmethod
    def from_dat_file(cls, filepath: str):
        """
        Reads meteorological data from a dat file and returns an instance of the
         Meteo class.

        Parameters:
        -----------
        filepath : str
            The path to the dat file.

        Returns:
        --------
        Meteo:
            An instance of the Meteo class.
        """
        met = pd.read_csv(filepath, sep=r"\s+")
        return cls(
            *(
                met[c].values
                for c in ["z", "u", "v", "w", "q", "T", "dir", "th", "tke"]
            )
        )

    def write_meteo_dat_file(self, file_path: str = "p2p/inputs/meteo.dat"):
        with open(file_path, "w") as f:
            f.write(
                f"{'z':15}"
                f"{'u':15}"
                f"{'v':15}"
                f"{'w':15}"
                f"{'q':15}"
                f"{'T':15}"
                f"{'dir':15}"
                f"{'th':15}"
                f"{'tke':15}\n"
            )
            for z, u, v, w, q, T, dir, th, tke in zip(
                self.zm,
                self.um,
                self.vm,
                self.wm,
                self.qq,
                self.temp,
                self.dir,
                self.th,
                self.tke,
            ):
                f.write(
                    f"{z:.3f}"
                    f"{u:15.3f}"
                    f"{v:15.3f}"
                    f"{w:15.3f}"
                    f"{q:15.3f}"
                    f"{T:15.3f}"
                    f"{dir:15}"
                    f"{th:15.3f}"
                    f"{tke:15.3f}\n"
                )

    @classmethod
    def from_sensor_file(
        cls,
        sensor_file_path: str,
        timestamp: timelike,
        bearing: float,
        qq: float = 0.05,
        idaweb_file_path: Optional[str] = None,
        tke: Optional[float] = None,
        extrapolate: Optional[bool] = True,
    ):
        """
        Reads meteorological data from sensor file and returns an instance of the
        Meteo class.

        Parameters:
        -----------
        cls: Meteo
            Class name.
        sensor_file_path : str
            Path to the sensor file.
        timestamp : timelike
            Timestamp for which meteo data is required.
        bearing : float
            Bearing of the flight.
        qq : float, optional (default=0.05)
            Mean turbulence velocity.
        idaweb_file_path : str, optional (default=None)
            Path to the idaweb file. This can either be a CSV or a pickle and the file
            type will be determined based on the extension. If None, the default path is
            used.
        tke : float, optional (default=None)
            Turbulent kinetic energy.
        extrapolate : bool, optional (default=True)
            If True, the wind speed and direction are extrapolated to different
            altitudes.

        Returns:
        --------
        Meteo
            An instance of the Meteo class.
        """
        from_idaweb = Meteo.from_idaweb_klo(
            timestamp,
            bearing,
            qq,
            idaweb_file_path,
            tke,
            extrapolate,
        )
        df = pd.read_parquet(sensor_file_path)
        df.index = df.index.tz_localize(timestamp.tz)

        start_time = timestamp - pd.Timedelta(minutes=2)
        stop_time = timestamp + pd.Timedelta(minutes=2)

        df["wind_speed"] = df["wind_speed"] * KNTS2MPS
        meteo_time = df.query("@start_time <= index <= @stop_time")
        wind_vel = meteo_time["wind_speed"].median()
        wind_dir = meteo_time["wind_direction"].median()
        wind_vel_std = meteo_time["wind_speed"].std()
        wind_dir_std = meteo_time["wind_direction"].std()

        alpha = np.deg2rad(wind_dir - bearing)
        crosswind = -wind_vel * np.sin(alpha)
        headwind = wind_vel * np.cos(alpha)
        if tke is None:
            tke = compute_tke(
                wind_vel,
                wind_vel_std,
                wind_dir,
                wind_dir_std,
                n=10000,
            )
        if extrapolate:
            headwinds = []
            crosswinds = []
            tkes = []
            alts = [2, 5, 10, 20, 50, 100, 200, 300, 400, 500]
            n = len(alts)
            for alt in alts:
                headwinds.append(extrapolate_wind_at_alt(alt, 10, headwind))
                crosswinds.append(extrapolate_wind_at_alt(alt, 10, crosswind))
                tkes.append(
                    compute_tke(
                        extrapolate_wind_at_alt(alt, 10, wind_vel),
                        wind_vel_std,
                        wind_dir,
                        wind_dir_std,
                        n=10000,
                    )
                )
            from_idaweb.um = headwinds
            from_idaweb.vm = crosswinds
            from_idaweb.dir = [wind_dir] * n
            from_idaweb.tke = tkes
            return from_idaweb

        from_idaweb.um = headwind
        from_idaweb.vm = crosswind
        from_idaweb.dir = wind_dir
        from_idaweb.tke = tke
        return from_idaweb

    def get_wind_speed(self, alt: float):
        um = np.interp(alt, self.zm, self.um)
        vm = np.interp(alt, self.zm, self.vm)
        wm = np.interp(alt, self.zm, self.wm)
        return np.sqrt(um**2 + vm**2 + wm**2)

    def get_wind_dir(self, alt: float):
        return np.interp(alt, self.zm, self.dir)


@dataclass
class EDR:
    """
    Class representing the edr values for different heights.

    Parameters:
    -----------
    zm : Union[float, Iterable[float]]
        Height(s) at which edr values are measured.
    edr : Union[float, Iterable[float]]
        EDR values at corresponding height(s).

    Methods:
    --------
    from_dat_file(cls, filepath):
        Class method to read edr data from a file and create an instance of EDR
        class.

    write_edr_dat_file(self, file_path):
        Method to write the edr data to a file.
    """

    zm: Union[float, Iterable[float]]
    edr: Union[float, Iterable[float]]

    def __repr__(self):
        return pd.DataFrame(self.__dict__).__repr__()

    def _repr_html_(self):
        return pd.DataFrame(self.__dict__)._repr_html_()

    @classmethod
    def from_meteo(cls, meteo_data: Meteo):
        """
        Class method to create an instance of EDR class from a Meteo instance.

        Parameters:
        -----------
        cls: EDR
            Class name.
        meteo: Meteo
            An instance of Meteo class.

        Returns:
        --------
        EDR
            An instance of EDR class.
        """
        return cls(meteo_data.zm, compute_edr(meteo_data.tke))

    @classmethod
    def from_dat_file(cls, filepath):
        """
        Class method to read edr data from a file and create an instance of EDR
        class.

        Parameters:
        -----------
        cls: EDR
            Class name.
        filepath : str
            Path of the dat file to read the edr data.

        Returns:
        --------
        EDR
            An instance of EDR class.
        """
        edr = pd.read_csv(
            filepath,
            sep=r"\s+",
            header=None,
            names=["z", "edr"],
        )
        return cls(*(edr[c].values for c in ["z", "edr"]))

    def write_edr_dat_file(self, file_path: str = "p2p/inputs/EDR.dat"):
        """
        Method to write the edr data to a file.

        Parameters:
        -----------
        self: EDR
            An instance of EDR class.
        file_path : str, optional (default="p2p/EDR.dat")
            Path of the dat file to write the edr data.
        """

        # if not os.path.exists(file_path):
        #     os.makedirs(file_path)

        with open(file_path, "w") as f:
            if isinstance(self.zm, float):
                f.write(f"{self.zm:15.3f} {self.edr:15.3f}\n")
            elif isinstance(self.edr, Iterable):
                for z, edr in zip(
                    self.zm,
                    self.edr,
                ):
                    f.write(f"{z:.3f} {edr:15.8f}\n")
            else:
                raise ValueError(
                    "zm and edr should be either float or Iterable"
                )


def extrapolate_wind_at_alt(
    z: float,
    z1: float,
    v1: float,
    z0: float = 0.03,
) -> float:
    """Extrapolate wind at altitude `z` from wind at altitude `z1`.

    Parameters
    ----------
    z: float
        Altitude at which to extrapolate wind speed.
    z1: float
        Altitude at which wind speed is known.
    v1: float
        Wind speed at altitude `z1`.
    z0: float, default: 0.03
        Roughness length.

    Returns
    -------
    float
        Extrapolated wind speed at altitude `z`.
    """
    # can be approximated as 2/3 to 3/4 of the average height of the obstacles
    d = 0.1  # displacement height (m)
    num = np.log((z - d) / z0)
    den = np.log((z1 - d) / z0)
    return v1 * num / den


def compute_tke(
    wind_vel_mean: float,
    wind_vel_std: float,
    wind_dir_mean: float,
    wind_dir_std: float,
    n: int = 10000,
) -> float:
    """Compute turbulence kinetic energy using bruteforce.

    I am not even sure this is valid but I tried to implement something...

    Parameters
    ----------
    wind_vel_mean : float
        Mean wind velocity [m/s]
    wind_vel_std : float
        Standard deviation of wind velocity [m/s]
    wind_dir_mean : float
        Mean wind direction [deg]
    wind_dir_std : float
        Standard deviation of wind direction [deg]
    n : int, optional
        Number of samples to draw, by default 10000

    Returns
    -------
    float
        Turbulence kinetic energy [m^2/s^2]
    """

    wind_vel = st.norm(loc=wind_vel_mean, scale=wind_vel_std).rvs(n)
    wind_dir = st.norm(loc=wind_dir_mean, scale=wind_dir_std).rvs(n)
    alpha_0 = np.deg2rad(wind_dir_mean)
    alpha = np.deg2rad(wind_dir)
    cw0 = wind_vel_mean * np.sin(alpha_0)
    hw0 = wind_vel_mean * np.cos(alpha_0)
    crosswind = -wind_vel * np.sin(alpha)
    headwind = wind_vel * np.cos(alpha)
    return (
        0.5
        * (((crosswind - cw0) ** 2).sum() + ((headwind - hw0) ** 2).sum())
        / n
    )


def compute_edr(tke: Union[Iterable, float], L: float = 311):
    """Compute the eddy dissipation rate.

    Parameters
    ----------
    tke : float
        Turbulent kinetic energy.
    L : float
        Length scale.

    Returns
    -------
    edr : float
        Eddy dissipation rate.

    """
    # Frech et al. (2007) => L = 311m
    # Ringley et al. (2007) => L = 336m

    if isinstance(tke, Iterable):
        return np.array([t ** (3 / 2) / L for t in tke])
    return tke ** (3 / 2) / L
