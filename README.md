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

# Limitations/Simplifications/Approximations
## Velocity field
The velocity field caused by the generator aircraft is calculated using the Biot-Savart Law for a constant strength finite length vortex filament. To retrieve this filament at each timestep, two points from the wake trajectory are taken with the filament representing the shortest connection between those two points. The first point represents the wake position at the current timestep, while the second point refers to the wake position at the last timestep. The calculated velocity field is therefore one timestep shorter than the follower trajectory and generator wake trajectory data fed into the function. Since the Biot-Savart Law produces singularities when the distance between vortex and the point of interest shrinks to zero. For that reason, a cutoff distance is set to the assumed wake vortex radius. For all points that meet this criterion the cutoff distance replaces the actual distance, limiting the induced velocities.
## Follower geometry
The follower aircraft is modeled through its lifting surfaces. The fuselage is not fed into the vortex lattice method. The lifting surfaces are modeled through their geometric properties (span, mean aerodynamic chord, geometric incidence angle). The geometric incidence angle is set to a constant three degrees for the wing and to negative three degrees for the horizontal stabilizer. A trim calculation to determine the correct incidence angle cannot be performed without the detailed mass distribution. The chord, sweep and airfoil distribution is also assumed to be a constant value. The chord is calculated from the geometric root and tip chord by the formula for the mean aerodynamic chord. The twodimensional airfoil is represented through its camber line. The airfoil thickness is not modeled in vortex lattice methods. For all aircraft, a generic supercritical airfoil (NASA SC(2)-0712) is used for the main wing while a inverted NACA23012 airfoil is assigned ot the horizontal stabilizer. The vertical stabilizer airfoil is a symmetric airfoil, similar to a flat plate (because of the camber line representation). Control surfaces are not included. Due to stability issues with the vortex lattice method, the dihedral angle is not included in the calculation.  
## Vortex Lattice Method
Vortex Lattice Methods build on potential theory. For this reason, three main approximations are passed down. Firstly, the flow field is assumed to be irrotational, i.e. the curl of the velocity vector field is zero. This can typically be assumed in free flow problems. Additionally, inviscid flow is assumed. The absence of viscosity also eliminates the no slip condition on surfaces. Without viscosity and rotational flow, the presence of a boundary layer is excluded from the flow problem. Lastly, potential theory assumes incompressible flow, i.e. the total derivative of the density is zero.
Due to these approximations, vortex lattice methods do not predict separation or turbulence. The absence of separation is a reasonable approximation as long as the angle of attack is limited to the linear region of the lift curve. For very strong wake scenarios, high induced velocities normal to the wing surface can violate these assumptions.
## Summary
In the following scenarios, the solution can lose accuracy:
- High angles of attack
- Modified lift distribution due to large control surface deflection
- Large differences in geometry between simulation and reality
- Extreme gust loads
- 
