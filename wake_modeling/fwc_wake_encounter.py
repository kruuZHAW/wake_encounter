import os
import click
import shutil

from pathlib import Path
from typing import Dict, List, Optional, Tuple, Union

import pandas as pd

import fwc_velocity_field
from fwc_get_settings import define_simulation_settings
from fuselage_wing_configuration import Fuselage_Wing_Configuration

rho = 1.225 #air density

# run_id              = "random_traj"
# typecode            = 'A320'
# calc_new_vl_field   = True
# wake                = '../data/modeling_inputs/wakes_df 3.parquet'
# trajectory          = '../data/modeling_inputs/encounter_ID_1_101lat_3vert_250speed_22s.parquet'
# velocity_input_path = '../data/modeling_inputs/velocity_field/'
# aircraft_db         = pd.read_parquet('../data/aircraft_db.parquet')

#Example of use:
# python fwc_wake_encounter.py 1 

@click.command()
@click.argument("out_path", type=str) 
@click.argument("run_id", type=int) # simulation id

@click.option('--wake_id', default=0, type=int) #Wake scenario to be used
@click.option('--v_field', default=True, type=bool) # calculation of a new velocity field
@click.option('--typecode', default="A320", type=str) # typecode of trailer aircraft

def main(
    out_path: str,
    run_id: int,
    wake_id:int,
    v_field: bool,
    typecode: str,
):
    
    #Data path
    wake_path = os.path.join(out_path, "wakes", str(wake_id), "wakes_df.parquet")
    traj_path = os.path.join(out_path, "encounters", str(run_id), "encounter_df.parquet")
    v_field_path = os.path.join(out_path, "encounters", str(run_id), "velocity_field")
    aircraft_db_path = os.path.join(os.getcwd(), os.pardir, "data", "aircraft_database", "aircraft_db.parquet")
    
    aircraft_db = pd.read_parquet(aircraft_db_path)
    
    if typecode not in aircraft_db.index:
        raise ValueError(f"Typecode '{typecode}' is not found in the aircraft database.")

    row = aircraft_db.loc[typecode]
    s_ref=row.wing_span*row.wing_mac
    c_ref=row.wing_mac
    
    ### MODIFY TO GO TO RIGHT LOCATION ###
    # Define basic case parameters
    case_name = typecode + "_" + str(run_id)
    case_route = os.path.join(out_path, "encounters", str(run_id))
    output_route =  os.path.join(out_path, "encounters", str(run_id))
    post_route = os.path.join(output_route, case_name, "savedata", case_name + '.data.h5')

    # velocity field
    if v_field:
        if os.path.exists(v_field_path):
            shutil.rmtree(v_field_path)
        os.makedirs(v_field_path)
        V_inf,time = fwc_velocity_field.main(wake_path, traj_path, v_field_path, v_field)
    else:
        V_inf,time = fwc_velocity_field.main(wake_path, traj_path, v_field_path, v_field)

    # Create an instance of Fuselage_Wing_Configuration
    fuselage_wing_config = Fuselage_Wing_Configuration(case_name, case_route, output_route)

    fuselage_wing_config.init_aeroelastic(
        half_wingspan=row.wing_span/2,  # Half of the wingspan
        half_htailspan=row.htail_span/2,
        vtailspan=row.vtail_span,
        fuselage_length=row.tail_position_x,  # Length of the fuselage/lever arm
        chord_wing=row.wing_mac,
        chord_htail=row.htail_mac,
        chord_vtail=row.vtail_mac,
        num_chordwise_panels=4,  # Number of chordwise panels for the wing
        elastic_axis=0.5,  # Elastic axis
        sweep_wing = row.wing_sweep,
        sweep_htail = row.htail_sweep,
        sweep_vtail = row.vtail_sweep,
        dihed_wing = 0.,
        dihed_htail = 0.,
        alpha_zero_deg_wing=3.0,  # Initial angle of attack
        alpha_zero_deg_htail=-3.0,
        max_radius=1.5,  # Max radius of the fuselage (for cylindrical shape)
        lifting_only=True,  # Include fuselage interactions
        airfoil_wing=row.wing_airfoil,
        airfoil_htail='naca23012_neg',
        vertical_wing_position=0.,
        vertical_htail_position=row.tail_position_z - row.wing_position_z
    )

    fuselage_wing_config.generate()

    # Define the flow sequence for SHARPy and the simulation settings
    flow_sequence = ['BeamLoader', 'AerogridLoader', 'StaticCoupled', 'DynamicCoupled', 'AerogridPlot', 'BeamPlot', 'AeroForcesCalculator', 'SaveData']

    settings = define_simulation_settings(
        flow=flow_sequence,
        model=fuselage_wing_config,  # The configuration instance
        alpha_deg=2.0,  # Angle of attack
        u_inf=V_inf,  # Freestream velocitys
        lifting_only=True,  # Include fuselage
        dt=0.1,  # Time step for dynamic simulations (if applicable)
        n_tsteps=int(time/1), # division by dt
        n_step=15,  # Number of load steps for static coupled simulations
        velocity_field_route = os.path.join(v_field_path, "velocity_field.xdmf")
    )

    # Create the SHARPy input file (a .sharpy file)
    fuselage_wing_config.create_settings(settings)

    fuselage_wing_config.run()

    results = fuselage_wing_config.calculate_aero_coefficients(post_route, V_inf, rho, s_ref, c_ref)
    results = pd.DataFrame(results)
    results.to_parquet(os.path.join(output_route, "results.parquet"))

    fuselage_wing_config.clean()

    #Delete simulation files
    # shutil.rmtree(os.path.join(output_route, case_name))

if __name__ == "__main__":
    main()
