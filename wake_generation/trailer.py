import os
import sys

#Give directory where the module P2P_base is.
parent_dir = os.path.abspath("/home/kruu/git_folder/")
root_dir = os.path.join(parent_dir, "P2P_base")
if parent_dir not in sys.path:
    sys.path.append(root_dir)
    
from random import sample

import utils.viz as viz

import pandas as pd
import numpy as np
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Union
import click

import warnings
warnings.simplefilter(action="ignore", category=FutureWarning)

def calculate_wake_locations(wake_df: pd.DataFrame) -> pd.DataFrame:
    
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

    data_wakes = pd.concat([data_l, data_r])
    data_wakes["size"] = 30

    # Add information about uncertainty
    for index, row in wake_df.iterrows():
        t = index
        gam = 1
        size = 1
        x_vals_ellipse_2s, y_vals_ellipse_2s = viz.get_ellipse_points(
            row["y_2s_l"],
            row["y_2s_r"],
            row["z_2s_lo"],
            row["z_2s_hi"],
            num_points=100,
        )
        x_vals_ellipse_3s, y_vals_ellipse_3s = viz.get_ellipse_points(
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
        
        data_wakes = pd.concat([data_wakes, df_2s, df_3s])
        
        return data_wakes
    

def generate_encounter(v: float, # m/s, velocity of the aircraft
                       t_target: int, # s, time when the encounter hit the wake gate
                       time_range: int, # s, duration of the wakes
                       timestep: float, # s, time resolution of wakes
                       x_target: float, # m, target x-coordinate
                       y_target: float, # m, target y-coordinate
                       z_target: float, # m, target altitude-coordinate
                       theta: float, # Lateral angle (w.r.t. x_axis i.e. generator speed vector)
                       phi: float, # vertical angle (w.r.t. the y-axis)
                       ) -> pd.DataFrame:
    
    # Velocity components
    v_x = v * np.cos(theta)
    v_y = v * np.sin(theta) * np.cos(phi)
    v_z = v * np.sin(theta) * np.sin(phi)
    
    # Initial position
    x_0 = x_target - v_x * t_target
    y_0 = y_target - v_y * t_target
    z_0 = z_target - v_z * t_target
    
    # Time range of the wake simulation
    extended_times = np.arange(1, time_range + 1, step= timestep) #Step to have refined timestamps for the trailer trajectory
    extended_times = np.round(extended_times, 6) # Same rounding as wakes_df index
    
    # Hit gate marker
    hit_mask = [False] * len(extended_times)
    hit_mask[t_target-1] = True
    
    # Positions for the extended time range
    positions = {
        'Time': extended_times,
        'X': x_0 + v_x * extended_times,
        'Y': y_0 + v_y * extended_times,
        'Z': z_0 + v_z * extended_times, 
        'hit_gate': hit_mask,
    }
    
    # Create the extended DataFrame
    encounter_trajectory = pd.DataFrame(positions).set_index('Time')
    
    return encounter_trajectory


@click.command()
@click.argument("out_path", type=str) 
@click.argument("run_id", type=int) # simulation id
@click.argument("timestep", type=float) # simulation id

@click.option('--wake_id', default=0, type=int) #Wake scenario to be used
@click.option('--aircraft_type', default="A320", type=str) # aircraft type of the trailer
@click.option('--crop_distance', default=2000, type=float) # distance from the closest wake from which we crop the trailer trajectory in m

def main(
    out_path:str,
    run_id:int,
    timestep:float,
    wake_id:int,
    aircraft_type:str,
    crop_distance:int,
):
    # fpath_wakes = os.path.abspath(os.path.join(os.getcwd(), "..", "data", "simulations", "wakes", str(wake_id), "wakes_df.parquet"))
    fpath_wakes = os.path.abspath(os.path.join(out_path, "wakes", str(wake_id), "wakes_df.parquet"))
    
    click.echo("Loading wakes...")
    wakes_df = pd.read_parquet(fpath_wakes)
    
    click.echo("Extracting gate...")
    data_wakes = calculate_wake_locations(wakes_df)

    click.echo("Generating random encounter...")
    t_range = wakes_df.index.max()
    
    # v = np.random.randint(60, 180) #in m/s: Operating speeds of an A320  
    # t_target = np.random.randint(1,t_range)
    # theta = np.random.randint(-90, 90)
    theta = sample([-90, -80, -60, -45, -30, -20, -10, 0, 10, 20, 30, 45, 60, 80, 90], 1)[0]
    # phi = np.random.randint(-10,10)
    
    v = 80
    t_target = 20
    # theta = 0
    phi = 0
    
    print(f"Speed: {v:.2f} m/s.")
    print(f"Target time: {t_target:.2f} s.")
    print(f"Lateral angle: {theta:.2f} degrees.")
    print(f"Vertical angle: {phi:.2f} degrees.")
    
    
    # Target coordinates
    # Left Wake at t = 20 : x = 0, y_l = -77.225, z_l = 1969.740
    # Right Wake at t = 20, x = 0, y_r = -40.082, z_r = 1969.740
    
    # x_target, y_target, z_target = wakes_df.x.iloc[0], np.random.uniform(data_wakes.y.min(), data_wakes.y.max()), np.random.uniform(data_wakes.z.min(), data_wakes.z.max())
    x_target, y_target, z_target = wakes_df.x.iloc[0], -40.082, 1969.740
    # x_target, y_target, z_target = wakes_df.x.iloc[0], np.random.uniform(-300, 100), np.random.uniform(1900, 2100)
    
    
    params = {
    "wake_id": [wake_id],
    "aircraft_type": [aircraft_type],
    "crop_distance": [crop_distance],
    "speed": [v],                
    "t_target": [t_target],     
    "theta": [theta],                
    "phi": [phi],                  
    "x_target": [x_target],          
    "y_target": [y_target],           
    "z_target": [z_target],                          
    }
    encounter_name = f"encounter_df.parquet"
    save_path = os.path.join(out_path, "encounters", str(run_id))
    os.makedirs(os.path.dirname(os.path.join(save_path, "param.parquet")), exist_ok=True)
    pd.DataFrame(params).to_parquet(os.path.join(save_path, "param.parquet"))
    
    print(f"X target: {x_target:.2f} m.")
    print(f"Y target: {y_target:.2f} m.")
    print(f"Z target: {z_target:.2f} m.")
    
    theta = np.radians(theta)
    phi = np.radians(phi)
    
    #generate encounter
    encounter = generate_encounter(v,
                       t_target,
                       t_range,
                       timestep,
                       x_target,
                       y_target,
                       z_target,
                       theta,
                       phi,)
    
    #calulate distance between trajectory and the center line of the wake (no matter the size of the tube)
    encounter["dist_left_wake"] = np.sqrt((encounter.Y - data_wakes.query("side == 'left'").y)**2 
        + (encounter.Z - data_wakes.query("side == 'left'").z)**2) 
    
    encounter["dist_right_wake"] = np.sqrt((encounter.Y - data_wakes.query("side == 'right'").y)**2 
        + (encounter.Z - data_wakes.query("side == 'right'").z)**2) 
    
    #Crop when the distance to the wake is less than 1000m for left wake or right wake
    encounter = encounter.query(f"(dist_right_wake <= {crop_distance}) | (dist_left_wake <= {crop_distance})")
    params["time_range"] = [encounter.index[0], encounter.index[-1]]

    encounter.to_parquet(os.path.abspath(os.path.join(save_path, encounter_name)))
    
    click.echo("Done !")

if __name__ == "__main__":
    main()
        
    
