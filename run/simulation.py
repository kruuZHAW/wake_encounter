
import argparse
import os
import warnings
from concurrent.futures import ProcessPoolExecutor
import subprocess

warnings.filterwarnings("ignore")

generator_path = "/home/kruu/git_folder/wake_encounter/wake_generation/generator.py"
encounter_path = "/home/kruu/git_folder/wake_encounter/wake_generation/trailer.py"
fwc_path = "/home/kruu/git_folder/wake_encounter/wake_modeling/fwc_wake_encounter.py"
output_path = "/store/kruu/wake_encounter_simulations/t20_phi0_speed80_discrete_theta_right_only_wake_A320"


parser = argparse.ArgumentParser()
parser.add_argument("--n_sim", type=int, default=None, help="Number of simulations")
parser.add_argument("--max-workers", type=int, default=None, help="Number of workers for multiprocessing")
args = parser.parse_args()

max_workers = args.max_workers
n_sim = args.n_sim

#Create directories:
if not os.path.exists(output_path):
    os.makedirs(output_path)
if not os.path.exists(os.path.join(output_path, "wakes")):
    os.makedirs(os.path.join(output_path, "wakes"))
if not os.path.exists(os.path.join(output_path, "encounters")):
    os.makedirs(os.path.join(output_path, "encounters"))


# Generate wakes if not already done
id_run_wake = 0
timestep = 0.1
if str(id_run_wake) not in os.listdir(os.path.join(output_path, "wakes")):
    try: 
        subprocess.run([
            "python", generator_path, output_path, str(id_run_wake),
            "2000", "1000", "10", "2", "130", "15", "1013", "10", "100", 
            "136000", "47.3", str(timestep), "--gpa", "0"
        ],stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, check=True) #stdout=subprocess.DEVNULL)
    except subprocess.CalledProcessError as e:
        print(f"Wake generation failed with exit code {e.returncode}")

# Function to process a single simulation
def process_simulation(id_run):
    print(f"Running Simulation {id_run}...")
    subprocess.run([
        "python", encounter_path, output_path, str(id_run),
        "--wake_id", "0", str(timestep), "--aircraft_type", "A320", "--crop_distance", "1500"
    ], stdout=subprocess.DEVNULL)

    subprocess.run([
        "python", fwc_path, output_path, str(id_run),
        "--wake_id", "0", "--v_field", "True", "--typecode", "A320"
    ], stdout=subprocess.DEVNULL)
    
# Multiprocessing to run simulations in parallel
if __name__ == "__main__":
    start_id = 0
    end_id = n_sim

    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        executor.map(process_simulation, range(start_id, end_id))
