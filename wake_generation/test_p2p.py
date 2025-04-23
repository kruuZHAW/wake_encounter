import os
import sys
import click
import shutil

import numpy as np
import pandas as pd
from traffic.data import opensky, airports

parent_dir = os.path.abspath("/home/kruu/git_folder/")
if parent_dir not in sys.path:
    sys.path.append(parent_dir)

root_dir = os.path.join(parent_dir, "P2P_base")
print(f"Parent dir: {parent_dir}")
print(f"Root dir: {root_dir}")
print(f"CWD: {os.getcwd()}")
    
import warnings
warnings.simplefilter(action="ignore", category=FutureWarning)

from P2P_base.wake import aircraft, meteo, wake

def main():
    flight = (
        opensky.history(
            "2018-07-04 08:58:28",
            stop="2018-07-04 09:06:36",
            icao24="4b1887",
            return_flight=True,
        )
        .distance(airports["LSZH"])
        .compute_xy("epsg:2056")
    )
    meteo_data = meteo.Meteo.from_metar('LSZH', timestamp=flight.stop, bearing=flight.at().track, extrapolate=True)
    fpath_edr = os.path.join(root_dir, "p2p", "EDR.dat")
    edr_data = meteo.EDR.from_dat_file(fpath_edr)
    
    wakes = wake.Wake.generate(aircraft=flight,
                                meteo=meteo_data,
                                edr=edr_data,
                                at_time=(flight.stop - pd.Timedelta('300s')),
                                path= os.path.abspath("/home/kruu/git_folder/wake_encounter/data/test_wakes"),
                                verbose=True) 
    
if __name__ == "__main__":
    main()