from dataclasses import dataclass
from datetime import datetime
from typing import Iterable, Optional, Tuple, Union
import numpy as np
import pandas as pd
import os
from os.path import join as opj
from traffic.core import Flight

timelike = Union[str, datetime, pd.Timestamp]

KNTS2MPS = 0.5144444444444445
FT2M = 0.3048


@dataclass
class Aircraft:
    """
    Contains information about an aircraft at a given time.

    Attributes
    ----------
    x: Iterable[float]
        x-coordinate of the aircraft's position
    y: Iterable[float]
        y-coordinate of the aircraft's position
    z: Iterable[float]
        aircraft's altitude above ground level
    speed_mps: Iterable[float]
        speed of the aircraft in meters per second
    gam0: Iterable[float]
        flight path angle of the aircraft in radians
    b0: Iterable[float]
        bank angle of the aircraft in radians
    gpa: Iterable[float]
        glide path angle of the aircraft in radians

    Methods
    -------
    from_dat_file(cls, filepath)
        Creates an Aircraft object from a dat file.
    from_flight(cls, flight, at_time=None)
        Creates an Aircraft object from a flight object.
    write_aircraft_dat_file(self, file_path)
        Writes aircraft data to a dat file.
    """

    x: Iterable[float]
    y: Iterable[float]
    z: Iterable[float]
    speed_mps: Iterable[float]
    gam0: Iterable[float]
    b0: Iterable[float]
    gpa: Iterable[float]

    def __repr__(self):
        return pd.DataFrame(self.__dict__).__repr__()

    def _repr_html_(self):
        return pd.DataFrame(self.__dict__)._repr_html_()

    @classmethod
    def from_dat_file(cls, filepath):
        """
        Create an Aircraft object from a dat file.

        Parameters
        ----------
        filepath : str
            The file path to the dat file.

        Returns
        -------
        Aircraft
            An Aircraft object with data from the dat file.
        """
        ac = pd.read_csv(filepath, sep=r"\s+")
        return cls(
            *(
                ac[c].values
                for c in ["x0", "y0", "z0", "uac", "gam0", "b0", "gpa"]
            )
        )

    @classmethod
    def from_flight(
        cls,
        flight: Flight,
        at_time: Optional[timelike] = None,
        **kwargs,
    ):
        """
        Create an Aircraft object from a flight object.

        Parameters
        ----------
        flight : Flight
            The flight object to use.
        at_time : str, datetime, pd.Timestamp, optional
            The time at which to create the Aircraft object.
            If None, uses the last recorded time.

        Returns
        -------
        Aircraft
            An Aircraft object with data from a Flight.
        """
        fdata = flight.at(at_time)
        assert fdata is not None
        alt_gnd = flight.landing_airport().altitude
        try:
            alt_agl = (fdata["corrected_altitude"] - alt_gnd) * FT2M
        except:
            alt_agl = (fdata["geoaltitude"] - alt_gnd) * FT2M

        speed_mps = fdata["groundspeed"] * KNTS2MPS
        gam0, b0 = get_gamma_n_b(speed_mps, flight.typecode, **kwargs)
        # TODO: compute gpa (glide path angle)
        return cls(
            [fdata.x], [fdata.y], [alt_agl], [speed_mps], [gam0], [b0], [3]
        )

    @classmethod
    def from_observation(
        cls,
        x: float,
        y: float,
        z: float,
        speed_mps: float,
        typecode: str,
        mass_kg: Optional[float] = None,
        **kwargs,
    ):
        """
        Create an Aircraft object from an observation.

        Parameters
        ----------
        x: float
            x-coordinate of the aircraft's position
        y: float
            y-coordinate of the aircraft's position
        z: float
            aircraft's altitude above ground level
        speed_mps: float
            speed of the aircraft in meters per second
        typecode: str
            The type code of the aircraft
        at_time : str, datetime, pd.Timestamp, optional
            The time at which to create the Aircraft object.
            Only required when an interaction with another aircraft is to be simulated.
        kwargs:
            Additional keyword arguments to pass to get_gamma_n_b

        Returns
        -------
        Aircraft
            An Aircraft object with data from the observation.
        """
        # TODO: compute gpa (glide path angle)
        gam0, b0 = get_gamma_n_b(
            speed_mps=speed_mps, typecode=typecode, mass_kg=mass_kg, **kwargs
        )
        return cls([x], [y], [z], [speed_mps], [gam0], [b0], [3])

    def write_aircraft_dat_file(
        self, file_path: str = "p2p/inputs/ac_init.dat"
    ):
        pos_x, pos_y = 0, 0
        with open(file_path, "w") as f:
            f.write(
                f"{'z0':>15}"
                f"{'gam0':>15}"
                f"{'b0':>15}"
                f"{'uac':>15}"
                f"{'x0':>15}"
                f"{'y0':>15}"
                f"{'gpa':>15}\n"
            )
            for z, speed_mps, gam0, b, gpa in zip(
                self.z,
                self.speed_mps,
                self.gam0,
                self.b0,
                self.gpa,
            ):
                f.write(
                    f"{z:15.3f}"
                    f"{gam0:15.3f}"
                    f"{b:15.3f}"
                    f"{speed_mps:15.3f}"
                    f"{pos_x:15.3f}"
                    f"{pos_y:15.3f}"
                    f"{gpa:15.3f}\n"
                )

    def to_args(
        self, args_order=["z", "gam0", "b0", "speed_mps", "x", "y", "gpa"]
    ):
        return " ".join(
            str(v) for v in [getattr(self, arg)[0] for arg in args_order]
        )


def get_gamma_n_b(
    speed_mps: float,
    typecode: Optional[str],
    mtow_pct: float = 1,
    mass_kg: Optional[float] = None,
    b_m: Optional[float] = None,
) -> Tuple[float, float]:
    """
    Returns the wake vortex strength (gamma) and wingspan (b) of the aircraft

    Parameters
    ----------
    speed_mps : float
        The true airspeed of the aircraft (in m/s)
    typecode : str
        The type code of the aircraft
    mtow_pct : float, optional
        The percentage of the MTOW to be used as the aircraft mass (default=0.7)
    mass_kg : float, optional
        The mass of the aircraft (in kg). If not provided, the mass will be
        calculated from the MTOW and the mtow_pct.
    b_m : float, optional
        The wingspan of the aircraft (in m). If not provided, the wingspan
        will be calculated from the typecode.

    Returns
    -------
    gam0 : float
        The wake vortex strength (gamma) of the aircraft
    b0 : float
        The wingspan (b) of the aircraft
    """

    # Load aircraft database
    root_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    ac_db = pd.read_excel(
        opj(
            root_dir,
            "data",
            "FAA-Aircraft-Char-DB-AC-150-5300-13B-App-2023-09-07.xlsx",
        ),
        sheet_name="ACD_Data",
        header=0,
    )
    if (typecode is None) & (b_m is None or mass_kg is None):
        raise ValueError(
            "Either typecode or both b_m and mass should be provided"
        )
    # Get wingspan of aircraft
    if b_m is None:
        b_m = (
            ac_db[ac_db["ICAO_Code"] == typecode][
                "Wingspan_ft_without_winglets_sharklets"
            ].values[0]
            * FT2M
        )
        # If no wingspan without winglets available, use with winglets
        if np.isnan(b_m):
            b_m = (
                ac_db[ac_db["ICAO_Code"] == typecode][
                    "Wingspan_ft_with_winglets_sharklets"
                ].values[0]
                * FT2M
            )

    if mass_kg is None:
        # Get MTOW of aircraft (in kg)
        mtow = (
            ac_db[ac_db["ICAO_Code"] == typecode]["MTOW_lb"].values[0]
            * 0.453592
        )

        # Calculate mass of aircraft (in kg)
        mass_kg = mtow * mtow_pct

    # Set air density (in kg/m^3)
    rho = 1

    # Set gravitational acceleration (in m/s^2)
    g = 9.81

    # Set spanwise load factor (assume elliptically loaded wing)
    s_l = np.pi / 4

    # Calculate gamma
    gam0 = (mass_kg * g) / (rho * speed_mps * b_m * s_l)

    # Calculate wingspan
    b0 = b_m * s_l

    return gam0, b0
