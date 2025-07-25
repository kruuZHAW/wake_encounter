## import declaration
import pandas as pd
import numpy as np
import scipy as sp
import os

def calcGrid(encounter: np.array, dt):
    # this function needs to be modified to calculate a grid
    # step 1: retrieve orientation matrix for each timestep of encounter
    # step 2: calculate all grid points, i.e. pointing away normal to LE (x),
    #           LE +- 5 meters (y), pointing normal to the wing up/down (z)
    # step 2a: initialize grid as a 4d matrix (t,x,y,z)
    # step 2b: multiply each timestep with R
    # step 2c: add the encounter coordinates for each timestep
    # in that way, only the grid needs to be exported
    # the grid can then be given to vortxl, and in that way, the vortxl function
    # can be vectorised
    
    r_t_dir = encounter[1:] - encounter[:-1]
    r_t_dir = r_t_dir[:,1:]

    norms = np.linalg.norm(r_t_dir, axis=1, keepdims=True)
    
    V_inf = np.mean(norms)/dt
    
    dir_t = r_t_dir/norms

    yaw = np.arctan2(dir_t[:,1], dir_t[:,0])
    pitch = np.arcsin(dir_t[:,2] / norms[1,:])

    # Compute cos and sin for pitch and yaw
    cos_pitch = np.cos(pitch)
    sin_pitch = np.sin(pitch)
    cos_yaw = np.cos(yaw)
    sin_yaw = np.sin(yaw)

    # Construct the R_z rotation matrices (yaw rotation around the z-axis)
    R_z = np.zeros((dir_t.shape[0], 3, 3))
    R_z[:, 0, 0] = cos_yaw
    R_z[:, 0, 1] = -sin_yaw
    R_z[:, 1, 0] = sin_yaw
    R_z[:, 1, 1] = cos_yaw
    R_z[:, 2, 2] = 1

    # Construct the R_y rotation matrices (pitch rotation around the y-axis)
    R_y = np.zeros((dir_t.shape[0], 3, 3))
    R_y[:, 0, 0] = cos_pitch
    R_y[:, 0, 2] = -sin_pitch
    R_y[:, 1, 1] = 1
    R_y[:, 2, 0] = sin_pitch
    R_y[:, 2, 2] = cos_pitch

   # Construct the R_y rotation matrices (roll rotation around the x-axis to transform from NWU to NED coordinate system)
    R_x = np.zeros((dir_t.shape[0], 3, 3))
    R_x[:, 0, 0] = 1
    R_x[:, 1, 1] = 1
    R_x[:, 2, 2] = -1

    # Combine the R_z and R_y matrices by matrix multiplication
    R = np.einsum('tij,tjk,tkm->tim', R_x, R_y, R_z)
    
    # Define the grid ranges and step size
    x = np.linspace(0, 100, 21)  # x from -100 to 0 with step 5
    y = np.linspace(-50, 50, 21)  # y from -20 to 20 with step 5
    z = np.linspace(-50, 50, 21)  # z from -10 to 10 with step 5
        
    t_dim = np.shape(encounter)[0]-1
    x_dim = len(x)
    y_dim = len(y)
    z_dim = len(z)
        
    # Create a meshgrid of points
    X, Y, Z = np.meshgrid(x, y, z, indexing='ij')

    # Flatten the grid points to (N, 3) where N = x_dim * y_dim * z_dim
    grid_points = np.vstack([X.ravel(), Y.ravel(), Z.ravel()]).T  # Shape (N, 3)

    # Multiply each rotation matrix (t, 3, 3) with grid points (N, 3)
    rotated_grids = np.einsum('tij,nj->tni', R, grid_points)

    # Add encounter coordinates to each timestep's grid points to make the grid follow the trajectory
    encounter_coords = encounter[1:, 1:4]  # Assuming encounter columns are time, x, y, z
    final_grids = rotated_grids + encounter_coords[:, np.newaxis, :]  # Broadcasting (t_dim, 1, 3) + (t_dim, N, 3)
    
    grids = final_grids.reshape(t_dim, x_dim, y_dim, z_dim, 3)
    return grids,R,V_inf

def vortxl(wake: np.array, grid: np.array, R: np.array):
    t_dim = np.shape(grid)[0]
    x_dim = np.shape(grid)[1]
    y_dim = np.shape(grid)[2]
    z_dim = np.shape(grid)[3]
    cutoff = 0.01 # cutoff radius in m to prevent infinite velocities 
    reg_norm = 5
    
    r_v_l = wake[:,[1,2,3]]
    r_v_r = wake[:,[1,5,6]]
    
    gamma_l = wake[:,4]
    gamma_r = -wake[:,7]
    
    # calculate r1 and r2 for both vortices
    r1_l = r_v_l[1:, np.newaxis, np.newaxis, np.newaxis, :] - grid[:, :, :, :, :] - np.array([50, 0, 0])
    r2_l = r1_l + np.array([200, 0, 0])
    r1_r = r_v_r[1:, np.newaxis, np.newaxis, np.newaxis, :] - grid[:, :, :, :, :] - np.array([50, 0, 0])
    r2_r = r1_r + np.array([200, 0, 0])
     
    ## calculate the vector cross product for each time step
    cross_product_l = np.cross(r1_l, r2_l, axis=-1)
    cross_product_r = np.cross(r1_r, r2_r, axis=-1)
    c_p_l_norm = np.linalg.norm(cross_product_l, axis=-1)
    c_p_r_norm = np.linalg.norm(cross_product_r, axis=-1)
       
    ## calculate distance between vortex filament and encounter for each time step
    r1_l_norm = np.linalg.norm(r1_l, axis=-1)
    r1_r_norm = np.linalg.norm(r1_r, axis=-1)
    r2_l_norm = np.linalg.norm(r2_l, axis=-1)
    r2_r_norm = np.linalg.norm(r2_r, axis=-1)
    
    # calculate wake trajectory
    r0_l = r2_l - r1_l
    r0_r = r2_r - r1_r
    
    # dot product 
    dot_product_l1 = np.einsum('...i,...i->...', r0_l, r1_l)
    dot_product_l2 = np.einsum('...i,...i->...', r0_l, r2_l)
    dot_product_r1 = np.einsum('...i,...i->...', r0_r, r1_r)
    dot_product_r2 = np.einsum('...i,...i->...', r0_r, r2_r)
    
    # check for zeros
    
    indZeros_r1_l = np.where(r1_l_norm == 0)
    indZeros_r2_l = np.where(r2_l_norm == 0)    
    indZeros_r1_r = np.where(r1_r_norm == 0)
    indZeros_r2_r = np.where(r2_r_norm == 0)
    
    indZeros_c_p_l = np.where(c_p_l_norm == 0)
    indZeros_c_p_r = np.where(c_p_r_norm == 0)
    
    r2_l_norm[indZeros_r2_l] = reg_norm
    r1_l_norm[indZeros_r1_l] = reg_norm
    r2_r_norm[indZeros_r2_r] = reg_norm
    r1_r_norm[indZeros_r1_r] = reg_norm
    c_p_l_norm[indZeros_c_p_l] = 1
    c_p_r_norm[indZeros_c_p_r] = 1
    
    ## apply Biot-Savart Law for a finite vortex filament
    K_l = (gamma_l[1:, np.newaxis, np.newaxis, np.newaxis] / (4 * np.pi * c_p_l_norm**2)) * \
          (dot_product_l1 / r1_l_norm - dot_product_l2 / r2_l_norm)
    K_r = (gamma_r[1:, np.newaxis, np.newaxis, np.newaxis] / (4 * np.pi * c_p_r_norm**2)) * \
          (dot_product_r1 / r1_r_norm - dot_product_r2 / r2_r_norm)
    
    u = K_l * cross_product_l[..., 0] + K_r * cross_product_r[..., 0]
    v = K_l * cross_product_l[..., 1] + K_r * cross_product_r[..., 1]
    w = K_l * cross_product_l[..., 2] + K_r * cross_product_r[..., 2]

    V_ind = np.stack([u, v, w], axis=-1) #
    
    # Reshape V_ind for efficient einsum operation
    V_ind_reshaped = V_ind.reshape(t_dim, -1, 3).transpose(0,2,1)  # Shape: (time, x*y*z, u,v,w)
    
    # Compute the inverse of the rotation matrix R to transform the velocities back to the aircraft coordinate system
    R_inv = R.transpose(0, 2, 1)
    
    # Perform the matrix-vector multiplication using einsum
    V_ind_rotated = np.einsum('tij,tjn->tin', R_inv, V_ind_reshaped)  # Shape: (t, u,v,w , x*y*z)
    # Transpose back and reshape to original dimensions
    V_ind = V_ind_rotated.transpose(0, 2, 1).reshape(t_dim, x_dim, y_dim, z_dim, 3)  # Shape: (time, x, y, z, u,v,w)
    
    # Set velocity to zero where the cutoff criterion is met
    mask = (c_p_l_norm < cutoff) | (c_p_r_norm < cutoff) | \
           (r1_l_norm < cutoff) | (r1_r_norm < cutoff) | \
           (r2_l_norm < cutoff) | (r2_r_norm < cutoff)
    
    V_ind[mask] = 0
    
    # Separate the velocities into U, V, W
    U = V_ind[..., 0]
    V = V_ind[..., 1]
    W = V_ind[..., 2]
    
    U = np.moveaxis(U, 0, -1)
    V = np.moveaxis(V, 0, -1)
    W = np.moveaxis(W, 0, -1)

    U = np.transpose(U, (0, 2, 1, 3))
    V = np.transpose(V, (0, 2, 1, 3))
    W = np.transpose(W, (0, 2, 1, 3))
    
    return U,V,W

def create_xdmf(U, V, W, save_path, nt, dt):
    """
    Creates XDMF and binary files for the velocity field data and saves them to the specified path.
    
    Parameters:
    - U, V, W: 4D numpy arrays containing velocity components.
    - save_path: str, the directory path where files will be saved.
    - nt: number of timesteps
    - dt: timestep
    """
    
    # Parameters
    dx, dy, dz = 5, 5, 5  # Spatial grid step sizes
    x_range, y_range, z_range = 100, 100, 100  # Spatial ranges
    x_origin, y_origin, z_origin = 0, -50, -50  # Spatial coordinate origin
    nx, ny, nz = int(x_range/dx+1), int(y_range/dy+1), int(z_range/dz+1)  # Number of spatial grid points
    
    # Generate spatial and time grids
    x = np.linspace(0, x_range, nx)
    y = np.linspace(0, y_range, ny)
    z = np.linspace(0, z_range, nz)
    [X, Y, Z] = np.meshgrid(x, y, z, indexing='ij')
    
    # Create binary files for each time step
    filename_prefix = 'velocity_field'
    for i in range(nt):
        filename_x = os.path.join(save_path, f'{filename_prefix}_t{i:02d}_ux.bin')
        filename_y = os.path.join(save_path, f'{filename_prefix}_t{i:02d}_uy.bin')
        filename_z = os.path.join(save_path, f'{filename_prefix}_t{i:02d}_uz.bin')
    
        U_slice = U[:, :, :, i]
        W_slice = W[:, :, :, i]
        V_slice = -V[:, :, :, i]  # Note the sign change for V
    
        U_slice.astype('float32').tofile(filename_x)
        W_slice.astype('float32').tofile(filename_y)
        V_slice.astype('float32').tofile(filename_z)
    
    # Create XDMF file
    xdmf_file = os.path.join(save_path, 'velocity_field.xdmf')
    with open(xdmf_file, 'w') as fid:
        fid.write('<?xml version="1.0" ?>\n')
        fid.write('<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd" []>\n')
        fid.write('<Xdmf Version="2.0">\n')
        fid.write('  <Domain>\n')
        fid.write(f'    <Topology TopologyType="3DRectMesh" Dimensions="{ny} {nz} {nx}"/>\n')
        fid.write('    <Geometry GeometryType="ORIGIN_DXDYDZ">\n')
        fid.write(f'      <DataItem Name="Origin" Dimensions="3" Format="XML"> 0.0 {y_origin} {z_origin} </DataItem>\n')
        fid.write(f'      <DataItem Name="Spacing" Dimensions="3" Format="XML"> {-dy} {dz} {dx} </DataItem>\n')
        fid.write('    </Geometry>\n')
        fid.write('    <Grid Name="GeneralGrid" GridType="Tree">\n')
        fid.write(f'      <Time>\n')
        fid.write(f'        <DataItem Dimensions="2" Format="XML"> 0.0 {dt} </DataItem>\n')
        fid.write('      </Time>\n')
    
        for i in range(nt):
            filename = f'{filename_prefix}_t{i:02d}'
            fid.write(f'      <Grid Name="SubGrid{i+1}" GridType="Uniform">\n')
            fid.write('        <Attribute Name="ux" Center="Node">\n')
            fid.write(f'           <DataItem Format="Binary" Precision="4">{filename}_ux.bin</DataItem>\n')
            fid.write('        </Attribute>\n')
            fid.write('        <Attribute Name="uy" Center="Node">\n')
            fid.write(f'           <DataItem Format="Binary" Precision="4">{filename}_uy.bin</DataItem>\n')
            fid.write('        </Attribute>\n')
            fid.write('        <Attribute Name="uz" Center="Node">\n')
            fid.write(f'           <DataItem Format="Binary" Precision="4">{filename}_uz.bin</DataItem>\n')
            fid.write('        </Attribute>\n')
            fid.write('      </Grid>\n')
    
        fid.write('    </Grid>\n')
        fid.write('  </Domain>\n')
        fid.write('</Xdmf>\n')
    return

def main(wake_path, trajectory_path, save_path, flag):
   encounter_df = pd.read_parquet(trajectory_path)
   encounter_df_light = encounter_df.drop(columns=['hit_gate','dist_left_wake','dist_right_wake'])
   encounter = np.array(encounter_df_light.reset_index())

   # encounter looks like this
   # t, x, y, z
   
   wake_df = pd.read_parquet(wake_path)
   wake_df = wake_df.drop(columns=['z_2s_lo','z_2s_hi','z_3s_lo','z_3s_hi','y_2s_l','y_2s_r',
                         'y_3s_l','y_3s_r','g_2s_lo','g_2s_hi','g_3s_lo','g_3s_hi'])
   indices_to_drop = wake_df.index.difference(encounter_df.index)
   wake_df_light = wake_df.drop(indices_to_drop)
   wake = np.array(wake_df_light.reset_index())
   index = [0,7,2,3,1,5,6,4]
   wake = wake[:,index]
   
   # rearranged wake so that it looks like this
   # t,x,y_l,z_l,gamma_l,y_r,z_r,gamma_r

   dt = np.round(encounter_df.index[1] - encounter_df.index[0], decimals=2)
   n_tsteps = len(encounter_df) - 1
   time = np.round(dt * n_tsteps, decimals=2)
   grid,R,V_inf = calcGrid(encounter, dt)
   print(f"dt: {dt}")
   print(f"Simulation time: {time}")
   print(f"Created grid: {grid.shape}")
   
   if flag:
       U,V,W = vortxl(wake,grid,R)
       
       U = np.concatenate([np.zeros_like(V[:,:,:,0])[..., np.newaxis], U], axis=3)
       U = U + V_inf
       
       V = np.concatenate([np.zeros_like(V[:,:,:,0])[..., np.newaxis], V], axis=3)

       W = np.concatenate([np.zeros_like(W[:,:,:,0])[..., np.newaxis], W], axis=3)

       create_xdmf(U, V, W, save_path, n_tsteps, dt)
   
   return V_inf,time, n_tsteps, dt
