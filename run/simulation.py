# First call generator.py to create wake
# Then call encounter.py to create trajectory
# Finally call fwc_wake_encounter.py to calculate wake response 

import subprocess
import os
import warnings

warnings.filterwarnings("ignore")

generator_path = "/home/kruu/git_folder/wake_encounter/wake_generation/generator.py"
encounter_path = "/home/kruu/git_folder/wake_encounter/wake_generation/trailer.py"
fwc_path = "/home/kruu/git_folder/wake_encounter/wake_modeling/fwc_wake_encounter.py"
# output_path = "/home/kruu/git_folder/wake_encounter/data/simulations"
output_path = "/backup/kruu/wake_encounter_simulations"

#TODO: If the ID_wake already exists, don't run that
id_run_wake = 0 
if str(id_run_wake) in os.listdir(os.path.join(output_path, "wakes")):
    print("Wake already computed")
else:
    print("Computing wakes")
    subprocess.run(["python", generator_path, output_path, str(id_run_wake), "2000", "1000", "10", "2", "130", "15", "1013", "10", "100", "136000", "47.3", "--gpa", "0"],
                    stdout=subprocess.DEVNULL)

#TODO: Try to have a lookup table: encounter is only one single point at a fiven time on the wake gate. 

for id_run in range(0,100):
    print(f"Running Simulation {id_run}...")
    subprocess.run(["python", encounter_path, output_path, str(id_run), "--wake_id", "0", "--aircraft_type", "A320", "--crop_distance", "2000"],
                    stdout=subprocess.DEVNULL)
    subprocess.run(["python", fwc_path, output_path, str(id_run), "--wake_id", "0", "--v_field", "True", "--typecode", "A320"],
                   stdout=subprocess.DEVNULL)