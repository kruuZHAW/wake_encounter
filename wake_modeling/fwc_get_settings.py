import numpy as np
import sharpy.utils.algebra as algebra

def define_simulation_settings(flow, model, alpha_deg, u_inf, 
                               dt=1,
                               rho = 1.225, 
                               lifting_only=True, 
                               nonlifting_only=False, 
                               phantom_test=False,
                               horseshoe=False, **kwargs):
    gravity = kwargs.get('gravity',False)
    nonlifting_body_interactions = not lifting_only and not nonlifting_only
    wake_length = kwargs.get('wake_length', 10)
    velocity_field_route = kwargs.get('velocity_field_route',)
    if not velocity_field_route:
        raise ValueError("Velocity field route has to be specified")
    
    # Other parameters
    if horseshoe:
        mstar = 1 
    else:
        mstar = wake_length*model.aero.num_chordwise_panels
    # numerics
    n_step = kwargs.get('n_step', 5)
    structural_relaxation_factor = kwargs.get('structural_relaxation_factor', 0.6)
    tolerance = kwargs.get('tolerance', 1e-6)
    fsi_tolerance = kwargs.get('fsi_tolerance', 1e-6)
    num_cores = kwargs.get('num_cores',4)

    if not lifting_only:
        nonlifting_body_interactions = True
    else:
        nonlifting_body_interactions = False
    settings = {}
    settings['SHARPy'] = {'case': model.case_name,
                        'route': model.case_route,
                        'flow': flow,
                        'write_screen': True,
                        'write_log': 'on',
                        'log_folder': model.output_route,
                        'log_file': model.case_name + '.log'}


    settings['BeamLoader'] = {'unsteady': False,
                                'orientation': algebra.euler2quat(np.array([0.,
                                                                            np.deg2rad(alpha_deg),
                                                                            0.]))}

    settings['BeamLoads']  = {}
    settings['LiftDistribution'] = {'coefficients': True}
    
    settings['NonLinearStatic'] = {'print_info': True,
                                'max_iterations': 150,
                                'num_load_steps': 1,
                                'delta_curved': 1e-1,
                                'min_delta': tolerance,
                                'gravity_on': gravity,
                                'gravity': 9.81}

    settings['StaticUvlm'] = {'print_info': True,
                            'horseshoe': horseshoe,
                            'num_cores': num_cores,
                            'velocity_field_generator': 'SteadyVelocityField',
                            'velocity_field_input': {'u_inf': u_inf,
                                                    'u_inf_direction': [1., 0, 0]},
                            'rho': rho,
                            'nonlifting_body_interactions': nonlifting_body_interactions,
                            'only_nonlifting': nonlifting_only,
                            'phantom_wing_test': phantom_test,
                            'ignore_first_x_nodes_in_force_calculation': kwargs.get('ignore_first_x_nodes_in_force_calculation', 0),
                            }

    settings['StaticCoupled'] = {'print_info': True,
                                'structural_solver': 'NonLinearStatic',
                                'structural_solver_settings': settings['NonLinearStatic'],
                                'aero_solver': 'StaticUvlm',
                                'aero_solver_settings': settings['StaticUvlm'],
                                'max_iter': 100,
                                'n_load_steps': n_step,
                                'tolerance': fsi_tolerance,
                                'relaxation_factor': structural_relaxation_factor,
                                'nonlifting_body_interactions': nonlifting_body_interactions}

    settings['AerogridLoader'] = {'unsteady': False,
                                'aligned_grid': True,
                                'mstar': mstar, #int(20/tstep_factor),
                                'wake_shape_generator': 'StraightWake',
                                'wake_shape_generator_input': {
                                    'u_inf': u_inf,
                                    'u_inf_direction': [1., 0., 0.],
                                    'dt': dt,
                                },
                            }
    if 'WriteVariablesTime' in flow:
        settings['WriteVariablesTime'] = {'cleanup_old_solution': True}
        if kwargs.get('writeCpVariables', False):
            settings['WriteVariablesTime']['nonlifting_nodes_variables'] = ['pressure_coefficients']                                                           
            settings['WriteVariablesTime']['nonlifting_nodes_isurf'] = np.zeros((model.structure.n_node_fuselage,))
            settings['WriteVariablesTime']['nonlifting_nodes_im'] =  np.zeros((model.structure.n_node_fuselage))
            settings['WriteVariablesTime']['nonlifting_nodes_in'] = list(range(model.structure.n_node_fuselage))
        if kwargs.get('writeWingPosVariables', False):
            settings['WriteVariablesTime']['structure_variables'] = ['pos']
            settings['WriteVariablesTime']['structure_nodes'] = list(range(model.structure.n_node-2))

    settings['NonliftingbodygridLoader'] = {}

    settings['AeroForcesCalculator'] = {'coefficients': False}
    settings['AerogridPlot'] = {'plot_nonlifting_surfaces': nonlifting_body_interactions}
    settings['BeamPlot'] = {}
    settings['SaveData'] = {}
    settings['NonLinearDynamicPrescribedStep'] = {'print_info': True,
                                            'max_iterations': 950,
                                            'delta_curved': 1e-1,
                                            'min_delta': 1e-6,
                                            'newmark_damp': 1e-4,
                                            'gravity_on': gravity,
                                            'gravity': 9.81,
                                            'num_steps': kwargs.get('n_tsteps',10),
                                            'dt': dt,
                                            }
    settings['RigidDynamicPrescribedStep'] = {'print_info': True,
                                            'max_iterations': 950,
                                            'delta_curved': 1e-1,
                                            'min_delta': 1e-6,
                                            'newmark_damp': 1e1,
                                            'gravity_on': gravity,
                                            'gravity': 9.81,
                                            'num_steps': kwargs.get('n_tsteps',10),
                                            'dt': dt,
                                            }
                                 
    settings['StepUvlm'] = {'print_info': True,
                                'num_cores': num_cores,
                                'convection_scheme': 0,
                                'velocity_field_generator': 'TurbVelocityField',
                                'velocity_field_input': {'turbulent_field': velocity_field_route,
                                                         'frozen': False,
                                                         'centre_y': False,
                                                         'offset': np.array([0, 0, 100]),
                                                         'store_field': False,
                                                         'periodicity': 'x'
                                                         },
                                'rho': rho,
                                'n_time_steps': kwargs.get('n_tsteps',10),
                                'dt': dt,
                                'cfl1': False,
                                'phantom_wing_test': phantom_test,
                                'nonlifting_body_interactions': not lifting_only,
                                'gamma_dot_filtering': 3,                                
                                'ignore_first_x_nodes_in_force_calculation': kwargs.get('ignore_first_x_nodes_in_force_calculation', 0),}

    dynamic_structural_solver = kwargs.get('structural_solver','RigidDynamicPrescribedStep')
    dynamic_aero_solver = kwargs.get('aero_solver', 'StepUvlm')
                                 
    settings['DynamicCoupled'] = {'structural_solver': dynamic_structural_solver,
                                    'structural_solver_settings': settings[dynamic_structural_solver],
                                    'aero_solver': 'StepUvlm',
                                    'aero_solver_settings': settings[dynamic_aero_solver],
                                    'fsi_substeps': kwargs.get('fsi_substeps', 20),
                                    'fsi_tolerance': fsi_tolerance,
                                    'relaxation_factor': kwargs.get('relaxation_factor',0.1),
                                    'minimum_steps': 1,
                                    'relaxation_steps': 150,
                                    'final_relaxation_factor': kwargs.get('final_relaxation_factor', 0.05),
                                    'n_time_steps': kwargs.get('n_tsteps',10),
                                    'dt': dt,
                                    'nonlifting_body_interactions': not lifting_only,
                                    'include_unsteady_force_contribution': kwargs.get('unsteady_force_distribution', True), 
                                    'postprocessors': ['BeamLoads'],
                                    'postprocessors_settings': {
                                                                'BeamLoads': {'csv_output': False},
                                                                },
                                }
    settings['DynamicUVLM'] = {'print_info': True,
                                'structural_solver': dynamic_structural_solver,
                                'structural_solver_settings': settings[dynamic_structural_solver],
                                'aero_solver': dynamic_aero_solver,
                                'aero_solver_settings': settings[dynamic_aero_solver],
                                'n_time_steps': kwargs.get('n_tsteps',10),
                                'dt': dt,
                                'include_unsteady_force_contribution': False,
                                'postprocessors': ['BeamLoads',
                                                       ],
                                'postprocessors_settings': {
                                                                'BeamLoads': {'csv_output': False},
                                                                },
                                }

    return settings
