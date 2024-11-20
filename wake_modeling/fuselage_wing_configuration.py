import os
import configobj
from fwc_structure import FWC_Structure
from fwc_aero import FWC_Aero
from fwc_fuselage import FWC_Fuselage
import sharpy.sharpy_main
import h5py
import numpy as np


class Fuselage_Wing_Configuration:
    """
        Fuselage_Wing_Configuration is a template class to create an
        aircraft model of a simple wing-fuselage configuration within
        SHARPy.

    """

    def __init__(self, case_name, case_route, output_route):
        self.case_name = case_name
        self.case_route = case_route
        self.output_route = output_route

        self.structure = None
        self.aero = None
        self.fuselage = None

        self.settings = None

    def init_aeroelastic(self, **kwargs):
        self.clean()
        self.init_structure(**kwargs)
        self.init_aero(**kwargs)
        if not kwargs.get('lifting_only', True):
            self.init_fuselage(**kwargs)

    def init_structure(self, **kwargs):
        self.structure = FWC_Structure(self.case_name, self.case_route, **kwargs)

    def init_aero(self, **kwargs):
        self.aero = FWC_Aero(self.structure, self.case_name, self.case_route,**kwargs)

    def init_fuselage(self, **kwargs):
        self.fuselage = FWC_Fuselage(self.structure, self.case_name, self.case_route,**kwargs)
    def generate(self):
        if not os.path.isdir(self.case_route):
            os.makedirs(self.case_route)
        self.structure.generate()
        self.aero.generate()
        if self.fuselage is not None:
            self.fuselage.generate()

    def create_settings(self, settings):
        file_name = self.case_route + '/' + self.case_name + '.sharpy'
        config = configobj.ConfigObj()
        config.filename = file_name
        for k, v in settings.items():
            config[k] = v
        config.write()
        self.settings = settings

    def clean(self):
        list_files = ['.fem.h5', '.aero.h5', '.nonlifting_body.h5', '.dyn.h5', '.mb.h5', '.sharpy', '.flightcon.txt']
        for file in list_files:
            path_file = self.case_route + '/' + self.case_name + file
            if os.path.isfile(path_file):
                os.remove(path_file)

    def run(self):
        sharpy.sharpy_main.main(['', self.case_route + '/' + self.case_name + '.sharpy'])  
        
    
    def calculate_aero_coefficients(self, data_path, u_inf, rho, s_ref, c_ref):
        """
        Calculates aerodynamic coefficients from SHARPy output data.
    
        Parameters:
        - data_path: str, path to the .h5 file containing aerodynamic forces.
        - u_inf: float, freestream velocity in m/s.
        - rho: float, air density in kg/m^3.
        - s_ref: float, reference area in m^2.
        - c_ref: float, reference chord in m.
    
        Returns:
        - coeffs: dict, containing calculated aerodynamic coefficients.
        """
        
        # Open the .h5 file
        with h5py.File(data_path, 'r') as h5file:
            # Navigate to the aero forces within each timestep
            timestep_info = list(h5file['data/aero/timestep_info'])
            # Initialize lists to store forces for each timestep
            loads = np.zeros((len(timestep_info)-1, 6))
    
      # Iterate through actual timestep keys (like '00001', '00002', etc.)
            for i in range(len(timestep_info)-1):
                counter_string = f"{i:05d}"
                # Access forces directly as a dataset
                forces = np.array(h5file['data/aero/timestep_info/' + counter_string + '/total_steady_body_forces'])

                loads[i,0] = (forces[0])  # Drag force
                loads[i,1] = (forces[1])  # Side force
                loads[i,2] = (forces[2])  # Lift force
                loads[i,3] = (forces[3])  # Rolling moment
                loads[i,4] = (forces[4])  # Pitching moment
                loads[i,5] = (forces[5])  # Yawing moment

        # Compute dynamic pressure
        q_inf = 0.5 * rho * u_inf**2
    
        # Calculate coefficients
        C_D = loads[:,0] / (q_inf * s_ref)     # Drag coefficient
        C_Y = loads[:,1] / (q_inf * s_ref)     # Side-force coefficient
        C_L = loads[:,2] / (q_inf * s_ref)     # Lift coefficient

        C_l = loads[:,3] / (q_inf * s_ref * c_ref)  # Rolling moment coefficient
        C_m = loads[:,4] / (q_inf * s_ref * c_ref)  # Pitching moment coefficient
        C_n = loads[:,5] / (q_inf * s_ref * c_ref)  # Yawing moment coefficient
    
        # Pack results in a dictionary
        coeffs = {
            'C_L': C_L,
            'C_D': C_D,
            'C_Y': C_Y,
            'C_l': C_l,
            'C_m': C_m,
            'C_n': C_n
        }
    
        return coeffs
