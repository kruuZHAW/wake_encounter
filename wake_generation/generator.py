"""Script generating the wake scenario"""

import os
import sys
import click
import shutil

import numpy as np
import pandas as pd

#Give directory where the module P2P_base is.
dir_p2p = os.path.abspath(os.path.join(os.getcwd(), os.pardir))
if  dir_p2p not in sys.path:
    sys.path.append(dir_p2p)
    
import warnings
warnings.simplefilter(action="ignore", category=FutureWarning)

from P2P_base.wake import aircraft, meteo, wake



CEL2K=273.15 # Celsuis to Kelvin
L = 0.0065 # Temperature lapse rate in K/m accroding to ISA: https://en.wikipedia.org/wiki/International_Standard_Atmosphere
R = 8.3144598 # Universal gas constant in J/(mol.K)
g = 9.81 # Gravitational acceleration in m/s^2
M = 0.028964 # Molar mass of Earth's air in kg/mol

def extrapolate_wind(
    alts: list[float],
    wind_vel: float,
    wind_dir: float, 
    vertical_vel: float,
    alt_sensor: float,
) -> tuple[list[float], list[float], list[float], list[float]]:

    headwinds = []
    crosswinds = []
    verticals = []
    
    crosswind = -wind_vel * np.sin(np.deg2rad(wind_dir)) 
    headwind = wind_vel * np.cos(np.deg2rad(wind_dir))

    for alt in alts:
        headwinds.append(meteo.extrapolate_wind_at_alt(alt, alt_sensor, headwind)) 
        crosswinds.append(meteo.extrapolate_wind_at_alt(alt, alt_sensor, crosswind)) 
        verticals.append(meteo.extrapolate_wind_at_alt(alt, alt_sensor, vertical_vel))

    directions = [wind_dir]*len(alts) #direction of the wind is constant
    
    return headwinds, crosswinds, verticals, directions
    
def extrapolate_temp(
    alts: list[float],
    temp: float,
    ) -> list[float]:
    
    temps = []
    
    for alt in alts:
        if alt <= 11000: 
            temps.append(temp - L*alt)
        elif alt > 11000 and alt <= 20000:
            temps.append(temp - L*11000)
        else:
            raise ValueError(
                        "You're flying too high"
                    )
    return temps
    
def calculate_potential_temp(
    alts: list[float],
    temp: float,
    temps: list[float],
    p_ref: float,
    ) -> tuple[list[float], list[float]]:
    
    p = p_ref*(temp / (temp + L*alts))**((g*M)/(R*L))
    th = np.array(temps)*((p_ref/p) ** 0.286)
    
    return p, th
    
    


@click.command()
@click.argument("run_id", type=int) # simulation id
@click.argument("alt_aircraft", type=float) # altitude above ground in m
@click.argument("spread", type=int) # spread of the range of altitudes around the aircraft in m
@click.argument("step", type=int) # step size for the altitude range in m
@click.argument("wind_vel", type=float) # measured wind velocity at the airport in m/s
@click.argument("wind_dir", type=float) # measured wind direction at the airport 
@click.argument("temp", type=float) # mesured temperature at the airport in °C
@click.argument("p_ref", type=float) # pressure at sea level measured at the airport in hPa
@click.argument("tke", type=float) # turbulent kinetic energy in m^2/s^2 (typically between 0 and 50)
@click.argument("speed", type=float) # speed of the aircraft in m/s
@click.argument("mass", type=float) # mass of the aircraft in kg
@click.argument("wingspan", type=float) # wingspan of the aircraft in m

@click.option('--wind_vertical_vel', default=0, type=float) # wind vertical velocity in m/s
@click.option('--alt_sensor', default=10, type=float) # meteo sensor altitude above ground in m
@click.option('--qq', default=0.05, type=float) # mean turbulent velocity in m/s
@click.option('--gpa', default=0, type=float) # glide path angle of the aircraft



def main(
    run_id: int,
    alt_aircraft: float,
    spread: int,
    step: int, 
    wind_vel: float,
    wind_dir: float,
    temp: float, 
    p_ref: float, 
    tke: float, 
    speed: float,
    mass: float,
    wingspan: float, 
    wind_vertical_vel: float, 
    alt_sensor: float,
    qq: float, 
    gpa: float, 
):
    
    click.echo("Building meteo data...")
    temp += CEL2K
    p_ref = p_ref*100
    
    alts = np.arange(alt_aircraft-spread,alt_aircraft+spread+step,step)
    headwinds, crosswinds, verticals, directions = extrapolate_wind(alts, wind_vel, wind_dir, wind_vertical_vel, alt_sensor)
    temps = extrapolate_temp(alts, temp)
    ps, ths = calculate_potential_temp(alts, temp, temps, p_ref)
    qqs = [qq]*len(alts)
    tkes = [tke]*len(alts)
    
    fpath_meteo = os.path.abspath(os.path.join(os.getcwd(), os.pardir, "data", "simulations", str(run_id), "wakes", "meteo.dat"))
    os.makedirs(os.path.dirname(fpath_meteo), exist_ok=True)
    
    with open(fpath_meteo, "w") as f:
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
            alts,
            headwinds,
            crosswinds,
            verticals,
            qqs,
            temps,
            directions,
            ths,
            tkes,
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
            
    click.echo("Building aircraft data...")
    s_l = np.pi/4 # spanwise load factor
    
    pos_x, pos_y, pos_z = 0, 0, alt_aircraft
    p_aircraft = ps[np.where(alts == alt_aircraft)[0][0]]
    temp_aircraft = temps[np.where(alts == alt_aircraft)[0][0]] 
    rho_aircraft = (p_aircraft * M) / (R * temp_aircraft) # air density at aircraft altitude
    gam0 = (mass * g)/(rho_aircraft * s_l * wingspan * speed)
    b0 = s_l * wingspan
    
    
    fpath_aircraft = os.path.abspath(os.path.join(os.getcwd(), os.pardir, "data", "simulations", str(run_id), "wakes", "ac_init.dat"))
    os.makedirs(os.path.dirname(fpath_aircraft), exist_ok=True)

    with open(fpath_aircraft, "w") as f:
        f.write(
            f"{'z0':>15}"
            f"{'gam0':>15}"
            f"{'b0':>15}"
            f"{'uac':>15}"
            f"{'x0':>15}"
            f"{'y0':>15}"
            f"{'gpa':>15}\n"
        )
        
        f.write(
            f"{pos_z:15.3f}"
            f"{gam0:15.3f}"
            f"{b0:15.3f}"
            f"{speed:15.3f}"
            f"{pos_x:15.3f}"
            f"{pos_y:15.3f}"
            f"{gpa:15.3f}\n"
        )
    
    meteo_data = meteo.Meteo.from_dat_file(fpath_meteo)
    aircraft_data = aircraft.Aircraft.from_dat_file(fpath_aircraft)
    edr_data = meteo.EDR.from_meteo(meteo_data)
    
    wakes = wake.Wake.generate(aircraft=aircraft_data,
                              meteo=meteo_data,
                              edr=edr_data,
                              path= os.path.abspath(os.path.join(os.getcwd(), os.pardir, "data", "simulations", str(run_id), "wakes")),
                            #   run_id=run_id,
                              verbose=True) 
    
    wakes.df.to_parquet(os.path.abspath(os.path.join(os.getcwd(), "..", "data", "simulations", str(run_id), "wakes", "wakes_df.parquet")))
    shutil.rmtree(os.path.abspath(os.path.join(os.getcwd(), os.pardir, "data", "simulations", str(run_id), "wakes", "inputs")))
    shutil.rmtree(os.path.abspath(os.path.join(os.getcwd(), os.pardir, "data", "simulations", str(run_id), "wakes", "results")))
    
    inputs = {
    "run_id": [run_id],                
    "alt_aircraft": [alt_aircraft],     
    "spread": [spread],                
    "step": [step],                  
    "wind_vel": [wind_vel],          
    "wind_dir": [wind_dir],           
    "temp": [temp],             
    "p_ref": [p_ref],            
    "tke": [tke],                 
    "speed": [speed],              
    "mass": [mass],             
    "wingspan": [wingspan],            
    "wind_vertical_vel": [wind_vertical_vel],   
    "alt_sensor": [alt_sensor],       
    "qq": [qq],                  
    "gpa": [gpa]                  
    }
    pd.DataFrame(inputs).to_parquet(os.path.join(os.getcwd(), os.pardir, "data", "simulations", str(run_id), "wakes", "param.parquet"))
    
    
    click.echo("Done !")

if __name__ == "__main__":
    main()
