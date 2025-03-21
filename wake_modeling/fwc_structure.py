import numpy as np
import h5py as h5


class FWC_Structure:

    """
        FWC_Structure contains all attributes to define the beam geometriy
        and discretisation and makes it accessible for SHARPy.
    """
    def __init__(self, case_name, case_route, **kwargs):
        self.n_elem_multiplier = kwargs.get('n_elem_multiplier', 1)

        self.route = case_route
        self.case_name = case_name

        self.n_node_elem = 3
        self.n_surfaces = 5
        self.n_material = 4 # wing + fuselage + htail + vtail

        self.half_wingspan              = kwargs.get('half_wingspan', 2.)
        self.half_htailspan             = kwargs.get('half_htailspan', 0.5)
        self.vtailspan                  = kwargs.get('vtailspan', 0.5)
        self.fuselage_length            = kwargs.get('fuselage_length', 10)
        self.offset_nose_wing_beam      = kwargs.get('offset_nose_wing', 1.)
        self.vertical_wing_position     = kwargs.get('vertical_wing_position', 0.)
        self.htail_offset               = kwargs.get('htail_offset', self.fuselage_length-self.offset_nose_wing_beam)
        self.vertical_htail_position    = kwargs.get('vertical_htail_position',2.)
        self.vtail_offset               = kwargs.get('vtail_offset', self.fuselage_length-self.offset_nose_wing_beam)
        self.vertical_vtail_position    = kwargs.get('vertical_vtail_position', 0.)
        self.root_chord_wing            = kwargs.get('root_chord_wing', 5.)
        
        self.global_offset = 0.25*self.root_chord_wing
        
        self.n_elem_per_wing    = kwargs.get('n_elem_per_wing', 10)
        self.n_elem_fuselage    = kwargs.get('n_elem_fuselage', 10)
        self.n_elem_per_htail   = kwargs.get('n_elem_htail',5)
        self.n_elem_vtail       = kwargs.get('n_elem_vtail',5)

        self.sigma = kwargs.get('sigma', 10)
        self.sigma_fuselage = kwargs.get('sigma_fuselage', 100.)

        self.sweep_wing     = np.deg2rad(kwargs.get('sweep_wing', 10.))
        self.sweep_htail    = np.deg2rad(kwargs.get('sweep_htail', 30.))
        self.sweep_vtail    = np.deg2rad(kwargs.get('sweep_vtail', 40.))
        self.dihed_wing     = np.deg2rad(kwargs.get('dihed_wing', 0.))
        self.dihed_htail    = np.deg2rad(kwargs.get('dihed_htail', 0.))

        self.enforce_uniform_fuselage_discretisation = kwargs.get(
            'enforce_uniform_fuselage_discretisation', False
        )
        self.fuselage_discretisation = kwargs.get('fuselage_discretisation', 'uniform')

        self.thrust = 0.

    def generate(self):
        """
            Function to set up all necessary parameter inputs to define the geometry and discretisation
            of the beam. Finally, informations are saved to an .h5 file serving as an input
            file for SHARPy.
        """
        self.set_element_and_nodes()
        self.initialize_parameters()

        self.set_stiffness_and_mass_properties()
        self.set_beam_properties_right_wing()
        self.mirror_wing_beam()

        self.app_forces[0] = [0, self.thrust, 0, 0, 0, 0]

        self.set_beam_properties_fuselage()
        self.set_beam_properties_right_htail()
        self.mirror_htail_beam()
        self.set_beam_properties_vtail()
        self.write_input_file()

    def set_element_and_nodes(self):
        """
            Based on the specified number of elements of the wing and fuselage, the number of nodes of each
            component and the total number of elements and nodes are defined here.
        """

        self.n_node_fuselage = self.n_elem_fuselage*(self.n_node_elem - 1)
        self.n_node_right_wing = self.n_elem_per_wing*(self.n_node_elem - 1) + 1
        # the left wing beam has one node less than the right one, since they share the center node
        self.n_node_left_wing = self.n_node_right_wing - 1
        # since the fuselage is already generated, the center node already exists
        self.n_node_right_htail = self.n_elem_per_htail*(self.n_node_elem - 1)
        self.n_node_left_htail = self.n_node_right_htail
        # the vtail should take one of the fuselage's nodes, so this is one less than for the wing
        self.n_node_vtail = self.n_elem_vtail*(self.n_node_elem - 1)

        self.n_node_wing_total = self.n_node_right_wing + self.n_node_left_wing
        self.n_node_htail_total = self.n_node_right_htail + self.n_node_left_htail
        self.n_node_tail_total =  self.n_node_htail_total + self.n_node_vtail
        self.n_node_fuselage_wing = self.n_node_wing_total + self.n_node_fuselage
        self.n_node_fuselage_wing_htail = self.n_node_fuselage_wing + self.n_node_htail_total
        self.n_node = self.n_node_fuselage + self.n_node_wing_total + self.n_node_tail_total

        self.n_elem_fuselage_wing = self.n_elem_fuselage + 2 * self.n_elem_per_wing
        self.n_elem_fuselage_wing_htail = self.n_elem_fuselage_wing + 2 * self.n_elem_per_htail
        self.n_elem = self.n_elem_fuselage_wing_htail + self.n_elem_vtail

        # vertical wing offset requires a connection element between wing and fuselage beam
        if not self.vertical_wing_position == 0:
            self.n_elem += 1
            self.n_node += 1

    def initialize_parameters(self):
        self.x = np.zeros((self.n_node, ))
        self.y = np.zeros((self.n_node, ))
        self.z = np.zeros((self.n_node, ))

        self.frame_of_reference_delta = np.zeros((self.n_elem, self.n_node_elem, 3))
        self.conn = np.zeros((self.n_elem, self.n_node_elem), dtype=int)

        self.beam_number = np.zeros((self.n_elem, ), dtype=int)
        self.boundary_conditions = np.zeros((self.n_node, ), dtype=int)
        self.structural_twist = np.zeros((self.n_elem, self.n_node_elem))

        self.app_forces = np.zeros((self.n_node, 6))


        self.stiffness = np.zeros((self.n_material, 6, 6))
        self.mass = np.zeros((self.n_material, 6, 6))
        self.elem_stiffness = np.zeros((self.n_elem, ), dtype=int)
        self.mass = np.zeros((self.n_material, 6, 6))
        self.elem_mass = np.zeros((self.n_elem, ), dtype=int)

    def set_stiffness_and_mass_properties(self):
        # Define aeroelastic properties
        ea = 1e10
        ga = 1e10
        gj = 1e10
        eiy = 1e10
        eiz = 1e10
        m_bar_main = 0.75
        j_bar_main = 0.075

        m_bar_fuselage = 0.3*1.5
        j_bar_fuselage = 0.08

        base_stiffness_main = self.sigma*np.diag([ea, ga, ga, gj, eiy, eiz])
        base_stiffness_fuselage = base_stiffness_main.copy()*self.sigma_fuselage
        base_stiffness_fuselage[4, 4] = base_stiffness_fuselage[5, 5]

        self.stiffness[0, ...] = base_stiffness_main
        self.stiffness[1, ...] = base_stiffness_fuselage
        self.stiffness[2, ...] = base_stiffness_main
        self.stiffness[3, ...] = base_stiffness_main

        self.mass[0, ...] = self.generate_mass_matrix(m_bar_main, j_bar_main)
        self.mass[1, ...] = self.generate_mass_matrix(m_bar_fuselage, j_bar_fuselage)
        self.mass[2, ...] = self.generate_mass_matrix(m_bar_main, j_bar_main)
        self.mass[3, ...] = self.generate_mass_matrix(m_bar_main, j_bar_main)

    def generate_mass_matrix(self, m_bar, j_bar):
        return np.diag([m_bar, m_bar, m_bar,
                j_bar, 0.5*j_bar, 0.5*j_bar])

    def set_beam_properties_right_wing(self):
        """
            Defines all necessary parameters to define the beam including node coordinate,
            elements with their associated nodes, frame of reference delta, boundary conditions
            such as reference node and free tips, stiffness and mass property ID for each element,
            and twist.
        """
        self.x[:self.n_node_right_wing] = np.linspace(0, self.half_wingspan*np.tan(self.sweep_wing), self.n_node_right_wing) + self.global_offset
        self.y[:self.n_node_right_wing] = np.linspace(0, self.half_wingspan, self.n_node_right_wing)
        self.z[:self.n_node_right_wing] += np.linspace(0, self.half_wingspan*np.tan(self.dihed_wing), self.n_node_right_wing)
        for ielem in range(self.n_elem_per_wing):
            self.conn[ielem, :] = ((np.ones((3, )) * ielem * (self.n_node_elem - 1)) +
                                [0, 2, 1])
            for ilocalnode in range(self.n_node_elem):
                self.frame_of_reference_delta[ielem, ilocalnode, :] = [-1.0, 0.0, 0.0]

        self.boundary_conditions[0] = 1
        self.boundary_conditions[self.n_node_right_wing-1] = -1 # free tip

    def mirror_wing_beam(self):
        """
            Mirrors the parameters from the beam representing the right free-flying wing
            for the left one.
        """

        self.x[self.n_node_right_wing:self.n_node_wing_total]  = self.x[1:self.n_node_right_wing]
        self.y[self.n_node_right_wing:self.n_node_wing_total]  = - self.y[1:self.n_node_right_wing]
        self.z[self.n_node_right_wing:self.n_node_wing_total]  = self.z[1:self.n_node_right_wing]
        self.frame_of_reference_delta[self.n_elem_per_wing:2*self.n_elem_per_wing, :, :] = self.frame_of_reference_delta[:self.n_elem_per_wing, :, :] * (-1)
        self.elem_stiffness[self.n_elem_per_wing:2*self.n_elem_per_wing] = self.elem_stiffness[:self.n_elem_per_wing]
        self.elem_mass[self.n_elem_per_wing:2*self.n_elem_per_wing] = self.elem_mass[:self.n_elem_per_wing]

        self.beam_number[self.n_elem_per_wing:2*self.n_elem_per_wing] = 1
        self.boundary_conditions[self.n_node_wing_total-1] = -1 # free tip
        self.conn[self.n_elem_per_wing:2*self.n_elem_per_wing, :] = self.conn[:self.n_elem_per_wing, :] + self.n_node_right_wing - 1
        self.conn[self.n_elem_per_wing, 0] = 0

    def set_x_coordinate_fuselage(self):
        if self.vertical_wing_position == 0:
            n_nodes_fuselage = self.n_node_fuselage + 1
        else:
            n_nodes_fuselage = self.n_node_fuselage
        if self.fuselage_discretisation == 'uniform':
            x_coord_fuselage = np.linspace(0, self.fuselage_length, n_nodes_fuselage) - self.offset_nose_wing_beam
        elif self.fuselage_discretisation == '1-cosine':
            x_coord_fuselage = np.linspace(0, 1, n_nodes_fuselage)
            x_coord_fuselage =  0.5*(1.0 - np.cos(x_coord_fuselage*np.pi))
            x_coord_fuselage *= self.fuselage_length
            x_coord_fuselage -= self.offset_nose_wing_beam
        else:
            raise "ERROR Specified fuselage discretisation '{}' unknown".format(self.fuselage_discretisation)
        self.idx_junction = self.find_index_of_closest_entry(x_coord_fuselage, self.x[0])

        self.idx_junction_global = self.idx_junction + self.n_node_wing_total
        if self.vertical_wing_position == 0:
            if self.enforce_uniform_fuselage_discretisation:
                self.x[:self.n_node_wing_total] += x_coord_fuselage[self.idx_junction]
            x_coord_fuselage = np.delete(x_coord_fuselage, self.idx_junction)
        self.x[self.n_node_wing_total:self.n_node_fuselage_wing] = x_coord_fuselage

    def adjust_fuselage_connectivities(self):
        idx_in_conn = np.where(self.conn ==  self.idx_junction_global)
        self.conn[idx_in_conn[0][0]+1:, :] -= 1

        if idx_in_conn[0][0] == 2:
            # if middle node, correct end node of element
            self.conn[idx_in_conn[0][0], 1] -= 1
        for i_match in range(np.shape(idx_in_conn)[1]):
            #several matches possible if junction node is not middle node
            self.conn[idx_in_conn[0][i_match], idx_in_conn[1][i_match]] = 0

    def add_additional_element_for_low_wing(self):
        self.x[-1] = self.x[0]
        self.y[-1] = self.y[0]
        self.z[-1] = self.vertical_wing_position / 2
        self.conn[-1, 0] = 0
        self.conn[-1, 1] = self.idx_junction_global
        self.conn[-1, 2] = self.n_node - 1
        self.elem_stiffness[-1] = 1
        self.elem_mass[-1] = 1
        self.beam_number[-1] = 3


    def set_beam_properties_fuselage(self):
        self.set_x_coordinate_fuselage()
        self.beam_number[self.n_elem_per_wing*2:self.n_elem_fuselage_wing] = 2
        for ielem in range(self.n_elem_per_wing * 2,self.n_elem_fuselage_wing):
            self.conn[ielem, :] = ((np.ones((3, )) * (ielem-self.n_elem_per_wing * 2) * (self.n_node_elem - 1)) +
                                [0, 2, 1])  + self.n_node_wing_total
            for ilocalnode in range(self.n_node_elem):
                self.frame_of_reference_delta[ielem, ilocalnode, :] = [0.0, 1.0, 0.0]
        self.elem_stiffness[self.n_elem_per_wing*2:] = 1
        self.elem_mass[self.n_elem_per_wing*2:] = 1

        if self.vertical_wing_position == 0:
            self.adjust_fuselage_connectivities()
        else:
            self.add_additional_element_for_low_wing()
        self.boundary_conditions[self.n_node_wing_total] = -1 # fuselage nose

    def set_beam_properties_right_htail(self):
       """
           Defines all necessary parameters to define the beam including node coordinate,
           elements with their associated nodes, frame of reference delta, boundary conditions
           such as reference node and free tips, stiffness and mass property ID for each element,
           and twist.
       """
       # the first line is added, however there is no sweep/dihedral yet
       self.x[(self.n_node_fuselage_wing - 1):(self.n_node_fuselage_wing  + self.n_node_right_htail)] = self.htail_offset + np.linspace(0, self.half_htailspan*np.tan(self.sweep_htail), (self.n_node_right_htail + 1)) + self.global_offset
       self.y[(self.n_node_fuselage_wing - 1):(self.n_node_fuselage_wing  + self.n_node_right_htail)] = np.linspace(0, self.half_htailspan, (self.n_node_right_htail + 1))
       self.z[(self.n_node_fuselage_wing - 1):(self.n_node_fuselage_wing  + self.n_node_right_htail)] += (-self.vertical_wing_position + self.vertical_htail_position) + np.linspace(0, self.half_htailspan*np.tan(self.dihed_htail), (self.n_node_right_htail + 1))
       for ielem in range(self.n_elem_fuselage_wing, self.n_elem_fuselage_wing + self.n_elem_per_htail):
           self.conn[ielem, :] = ((np.ones((3, )) * (ielem-self.n_elem_fuselage_wing) * (self.n_node_elem - 1)) +
                               [0, 2, 1]) + self.n_node_fuselage_wing - 1
           for ilocalnode in range(self.n_node_elem):
               self.frame_of_reference_delta[ielem, ilocalnode, :] = [-1.0, 0.0, 0.0]
               
       self.beam_number[self.n_elem_fuselage_wing:(self.n_elem_fuselage_wing + self.n_elem_per_htail)] = 3
       self.boundary_conditions[(self.n_node_fuselage_wing - 1 + self.n_node_right_htail)] = -1 # free tip

    def mirror_htail_beam(self):
        """
            Mirrors the parameters from the beam representing the right free-flying htail
            for the left one.
        """
        self.x[(self.n_node_fuselage_wing + self.n_node_right_htail):self.n_node_fuselage_wing_htail]  = self.x[self.n_node_fuselage_wing:(self.n_node_fuselage_wing + self.n_node_right_htail)]
        self.y[(self.n_node_fuselage_wing + self.n_node_right_htail):self.n_node_fuselage_wing_htail]  = - self.y[self.n_node_fuselage_wing:(self.n_node_fuselage_wing + self.n_node_right_htail)]
        self.z[(self.n_node_fuselage_wing + self.n_node_right_htail):self.n_node_fuselage_wing_htail]  = self.z[self.n_node_fuselage_wing:(self.n_node_fuselage_wing + self.n_node_right_htail)]
        self.frame_of_reference_delta[self.n_elem_fuselage_wing + self.n_elem_per_htail:self.n_elem_fuselage_wing_htail, :, :] = self.frame_of_reference_delta[self.n_elem_fuselage_wing:(self.n_elem_fuselage_wing + self.n_elem_per_htail), :, :] * (-1)
        self.elem_stiffness[self.n_elem_fuselage_wing + self.n_elem_per_htail:self.n_elem_fuselage_wing_htail] = self.elem_stiffness[:self.n_elem_per_htail]
        self.elem_mass[self.n_elem_fuselage_wing + self.n_elem_per_htail:self.n_elem_fuselage_wing_htail] = self.elem_mass[:self.n_elem_per_htail]

        self.beam_number[self.n_elem_fuselage_wing + self.n_elem_per_htail:self.n_elem_fuselage_wing_htail] = 4
        self.boundary_conditions[self.n_node_fuselage_wing_htail-1] = -1 # free tip
        self.conn[self.n_elem_fuselage_wing + self.n_elem_per_htail:self.n_elem_fuselage_wing_htail, :] = self.conn[self.n_elem_fuselage_wing:self.n_elem_fuselage_wing + self.n_elem_per_htail, :] + self.n_node_right_htail
        self.conn[self.n_elem_fuselage_wing + self.n_elem_per_htail, 0] = self.conn[self.n_elem_fuselage_wing,0]
        
    def set_beam_properties_vtail(self):
       """
           Defines all necessary parameters to define the beam including node coordinate,
           elements with their associated nodes, frame of reference delta, boundary conditions
           such as reference node and free tips, stiffness and mass property ID for each element,
           and twist.
       """
       # the first line is added, however there is no sweep/dihedral yet
       self.y[self.n_node_fuselage_wing_htail:self.n_node] = 0
       self.z[self.n_node_fuselage_wing_htail:self.n_node] = (-self.vertical_wing_position + self.vertical_vtail_position + self.vertical_htail_position) + np.linspace(self.vtailspan/self.n_node_vtail, self.vtailspan, (self.n_node_vtail))
       self.x[self.n_node_fuselage_wing_htail:self.n_node] = self.vtail_offset + np.linspace(self.z[self.n_node_fuselage_wing_htail]*np.tan(self.sweep_vtail),self.vtailspan*np.tan(self.sweep_vtail),self.n_node_vtail) + self.global_offset

       for ielem in range(self.n_elem_fuselage_wing_htail, self.n_elem):
           self.conn[ielem, :] = ((np.ones((3, )) * (ielem-self.n_elem_fuselage_wing_htail) * (self.n_node_elem - 1)) +
                               [0, 2, 1]) + self.n_node_fuselage_wing_htail - 1
           for ilocalnode in range(self.n_node_elem):
               self.frame_of_reference_delta[ielem, ilocalnode, :] = [-1.0, 0.0, 0.0]
               
       self.beam_number[self.n_elem_fuselage_wing_htail:(self.n_elem_fuselage_wing_htail + self.n_elem_vtail)] = 5
       self.boundary_conditions[(self.n_node_fuselage_wing_htail - 1 + self.n_node_vtail)] = -1 # free tip
       self.conn[self.n_elem_fuselage_wing_htail, 0] = self.conn[self.n_elem_fuselage_wing,0]


    def find_index_of_closest_entry(self, array_values, target_value):
        return np.argmin(np.abs(array_values - target_value))

    def set_thrust(self, value):
        self.thrust = value

    def write_input_file(self):
        """
            Writes previously defined parameters to an .h5 file which serves later as an
            input file for SHARPy.
        """
        with h5.File(self.route + '/' + self.case_name + '.fem.h5', 'a') as h5file:
            h5file.create_dataset('coordinates', data=np.column_stack((self.x, self.y, self.z)))
            h5file.create_dataset('connectivities', data=self.conn)
            h5file.create_dataset('num_node_elem', data=self.n_node_elem)
            h5file.create_dataset('num_node', data=self.n_node)
            h5file.create_dataset('num_elem', data=self.n_elem)
            h5file.create_dataset('stiffness_db', data=self.stiffness)
            h5file.create_dataset('elem_stiffness', data=self.elem_stiffness)
            h5file.create_dataset('mass_db', data=self.mass)
            h5file.create_dataset('elem_mass', data=self.elem_mass)
            h5file.create_dataset('frame_of_reference_delta', data=self.frame_of_reference_delta)
            h5file.create_dataset('boundary_conditions', data=self.boundary_conditions)
            h5file.create_dataset('beam_number', data=self.beam_number)

            h5file.create_dataset('structural_twist', data=self.structural_twist)
            h5file.create_dataset('app_forces', data=self.app_forces)
