# Wake Encounter
This repo contains several scripts that work together to setup a simulation for the VLM SHARPy. It is derived from the fuselage wing configuration (fwc) template.
## Contents
- "database": contains the aircraft data, as well as wake and encounter trajectorys.
- "fwc_velocity_field.py": creates a velocity field based on a saved wake and encounter trajectory. Also outputs the mean velocity along the trajectory.
- "fuselage_wing_configuration.py": top-level class interacting with the three sub-classes aeroelastic, aero and structure, stores methods for running the simulation. Produces the .sharpy file for the simulation.
- "fwc_structure.py": Sets up a beam representation of the aircraft including connectivity matrix. Structural parameters are set to represent a completely rigid aircraft, to omit the structural calculation in the VLM.Produces the .fem.h5 file for the simulation.
- "fwc_aero.py": Adds the vortex lattice based on the geometry stored in the structure class. Produces the .aero.h5 file for the simulation.
- "fwc_fuselage.py": Holds a class designed to set-up the non-lifting body. Currently not used to simplify the setup of the simulation and speed up computation time.
- "fwc_wake_encounter.py": Input script, calls all necessary commands to set up a simulation. Details below.
- "fwc_get_settings.py": Holds default settings for all SHARPy solvers involved in the simulation. Inside this script, the path to the velocity field is set. The number of cores is hardcoded here as well, and set to 4. If a cluster is used, the number of cores can be set here to improve computation time.

## fwc_wake_encounter
This script serves as the input/output for the simulation. It takes several inputs and outputs an array holding the force and moment coefficients with the size n_timesteps x 6.
The following default values are set:
```
typecode            = 'A320'
calc_new_vl_field   = True
wake                = './database/input/wakes_df 3.parquet'
trajectory          = './database/input/encounter_90deg_left_wake_75s.parquet'
velocity_input_path = './database/input/velocity_field/'
```
- typecode refers to the follower aircraft that experiences the wake. The following aircraft are implemented as of now: A319, A320, A333, A359, A380, A221, A223 and E170.
- calc_new_vl_field is a flag. If set to true, a new velocity field is calculated, otherwise only the infinite velocity is calculated.
- wake refers to the path where the wake-parquet is saved.
- trajectory refers to the path where the trajectory-parquet is saved.
- velocity_input_path refers to the path where the velocity field will be stored.

## Installation
To get the scripts to work, follow the instructions to install SHARPy (can be found here: https://ic-sharpy.readthedocs.io/en/latest/content/installation.html). Some additional packages are required to process the scripts, please refer to the environment.yml file for all required packages.
