from fuselage_wing_configuration import Fuselage_Wing_Configuration
import pandas as pd
import fwc_velocity_field
import os
import shutil

run_id              = "random_traj"
typecode            = 'A320'
calc_new_vl_field   = True
wake                = '../data/modeling_inputs/wakes_df 3.parquet'
trajectory          = '../data/modeling_inputs/encounter_ID_1_101lat_3vert_250speed_22s.parquet'
velocity_input_path = '../data/modeling_inputs/velocity_field/'
aircraft_db         = pd.read_parquet('../data/aircraft_db.parquet')

row = aircraft_db.loc[typecode]
s_ref=row.wing_span*row.wing_mac
c_ref=row.wing_mac

# Define basic case parameters
case_name = typecode + "_" + str(run_id)
case_route = "../data/output"
output_route = "../data/output"
post_route = output_route + "/" + case_name + '/savedata/' + case_name + '.data.h5'

# velocity field
rho = 1.225
if calc_new_vl_field:
    if os.path.exists(velocity_input_path):
        shutil.rmtree(velocity_input_path)
    os.makedirs(velocity_input_path)
    V_inf,time = fwc_velocity_field.main(wake, trajectory, velocity_input_path, calc_new_vl_field)
else:
    V_inf,time = fwc_velocity_field.main(wake, trajectory, velocity_input_path, calc_new_vl_field)

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

from fwc_get_settings import define_simulation_settings

# Define the flow sequence for SHARPy and the simulation settings
flow_sequence = ['BeamLoader', 'AerogridLoader', 'StaticCoupled', 'DynamicCoupled', 'AerogridPlot', 'BeamPlot', 'AeroForcesCalculator', 'SaveData']

settings = define_simulation_settings(
    flow=flow_sequence,
    model=fuselage_wing_config,  # The configuration instance
    alpha_deg=2.0,  # Angle of attack
    u_inf=V_inf,  # Freestream velocitys
    lifting_only=True,  # Include fuselage
    dt=1,  # Time step for dynamic simulations (if applicable)
    n_tsteps=int(time/1), # division by dt
    n_step=15  # Number of load steps for static coupled simulations
)

# Create the SHARPy input file (a .sharpy file)
fuselage_wing_config.create_settings(settings)

fuselage_wing_config.run()

results = fuselage_wing_config.calculate_aero_coefficients(post_route, V_inf, rho, s_ref, c_ref)
results = pd.DataFrame(results)
results.to_parquet(os.path.join(output_route, case_name + "_results.parquet"))

fuselage_wing_config.clean()

#Delete simulation files
shutil.rmtree(os.path.join(output_route, case_name))