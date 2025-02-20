import os
import numpy as np
import pandas as pd
import plotly.graph_objects as go
import plotly.offline as py_offline
import sys


def calculate_wake_locations(wake_df: pd.DataFrame, viz) -> pd.DataFrame:
    data_l = wake_df[["yl", "zl", "gam_l"]].rename(columns={"yl": "y", "zl": "z", "gam_l": "gam"})
    data_l["t"] = wake_df.index
    data_l["side"] = "left"

    data_r = wake_df[["yr", "zr", "gam_r"]].rename(columns={"yr": "y", "zr": "z", "gam_r": "gam"})
    data_r["t"] = wake_df.index
    data_r["side"] = "right"

    data_wakes = pd.concat([data_l, data_r])
    data_wakes["size"] = 30

    for index, row in wake_df.iterrows():
        t = index
        gam = 1
        size = 1
        x_vals_ellipse_2s, y_vals_ellipse_2s = viz.get_ellipse_points(row["y_2s_l"], row["y_2s_r"], row["z_2s_lo"], row["z_2s_hi"], num_points=100)
        x_vals_ellipse_3s, y_vals_ellipse_3s = viz.get_ellipse_points(row["y_3s_l"], row["y_3s_r"], row["z_3s_lo"], row["z_3s_hi"], num_points=100)
        df_2s = pd.DataFrame({"t": t, "gam": gam, "y": x_vals_ellipse_2s, "z": y_vals_ellipse_2s, "side": "uncert_2s", "size": size})
        df_3s = pd.DataFrame({"t": t, "gam": gam, "y": x_vals_ellipse_3s, "z": y_vals_ellipse_3s, "side": "uncert_3s", "size": size})
        data_wakes = pd.concat([data_wakes, df_2s, df_3s])
    
    return data_wakes

def generate_wake_plot(encounter_path: str, wake_path: str) -> go.Figure:
    sys.path.append('/home/kruu/git_folder/wake_encounter/P2P_base')
    import utils.viz as viz
    
    encounter_df = pd.read_parquet(encounter_path)
    wake_df = pd.read_parquet(wake_path)
    
    data_wakes = calculate_wake_locations(wake_df, viz)
    py_offline.init_notebook_mode()

    t_target = encounter_df.query("hit_gate == True").index.values[0]
    x_target = encounter_df.query("hit_gate == True").X.values[0]
    y_target = encounter_df.query("hit_gate == True").Y.values[0]
    z_target = encounter_df.query("hit_gate == True").Z.values[0]
    
    wake_color = 'rgba(255, 165, 0, 0.6)'
    plane_color = 'rgba(128, 128, 128, 0.3)'

    target_marker = {'size': 8, 'opacity': 0.9, 'color': 'red', 'symbol': 'x'}
    target_trace = go.Scatter3d(x=[x_target], y=[y_target], z=[z_target], mode='markers', marker=target_marker, name="Crossing point")
    
    x_plane = [wake_df.x.iloc[0]]*100
    y_plane = np.linspace(data_wakes.y.min(), data_wakes.y.max(), 100) 
    z_plane = np.linspace(data_wakes.z.min(), data_wakes.z.max(), 100) 
    plane_trace = go.Surface(x=x_plane, y=y_plane, z=np.tile(z_plane, (len(x_plane), 1)), opacity=0.5, colorscale=[[0, plane_color], [1, plane_color]], showscale=False)

    radius = 10
    x = np.linspace(-400, 400, 100)
    cyl_angle = np.linspace(0, 2 * np.pi, 50)

    left_wake_center = [data_wakes.query("side == 'left'").loc[t_target].y, data_wakes.query("side == 'left'").loc[t_target].z]
    left_wake_y = left_wake_center[0] + radius * np.cos(cyl_angle)
    left_wake_z = left_wake_center[1] + radius * np.sin(cyl_angle)
    left_wake_Yc, left_wake_Xc = np.meshgrid(left_wake_y, x)
    left_wake_Zc, _ = np.meshgrid(left_wake_z, x)
    left_wake_cylinder = go.Surface(x=left_wake_Xc, y=left_wake_Yc, z=left_wake_Zc, opacity=0.5, colorscale=[[0, wake_color], [1, wake_color]], showscale=False)

    right_wake_center = [data_wakes.query("side == 'right'").loc[t_target].y, data_wakes.query("side == 'right'").loc[t_target].z]
    right_wake_y = right_wake_center[0] + radius * np.cos(cyl_angle)
    right_wake_z = right_wake_center[1] + radius * np.sin(cyl_angle)
    right_wake_Yc, right_wake_Xc = np.meshgrid(right_wake_y, x)
    right_wake_Zc, _ = np.meshgrid(right_wake_z, x)
    right_wake_cylinder = go.Surface(x=right_wake_Xc, y=right_wake_Yc, z=right_wake_Zc, opacity=0.5, colorscale=[[0, wake_color], [1, wake_color]], showscale=False)

    trajectory_trace = go.Scatter3d(x=encounter_df['X'], y=encounter_df['Y'], z=encounter_df['Z'], mode='lines', line=dict(color='crimson', width=6), name='Aircraft Trajectory')
    layout = go.Layout(title='Wake Encounter Scenario', scene=dict(xaxis=dict(title='X (m)', range=[-300, 300]), yaxis=dict(title='Y (m)', range=[data_wakes.y.min() - 100, data_wakes.y.max() + 100]), zaxis=dict(title='Altitude (m)', range=[data_wakes.z.min() - 100, data_wakes.z.max() + 100]), bgcolor='rgba(230, 230, 230, 0.8)'), margin={'l': 0, 'r': 0, 'b': 0, 't': 40})

    figure = go.Figure(layout=layout)
    figure.add_traces([plane_trace, target_trace, trajectory_trace, left_wake_cylinder, right_wake_cylinder])
    
    return figure

def params_encounter(source_path: str) -> pd.DataFrame:
    directories = sorted([d for d in os.listdir(os.path.join(source_path, "encounters")) if d.isdigit()])
    dataframes = []

    for dir_name in directories:
        file_path = os.path.join(source_path, "encounters", dir_name, "param.parquet")
        
        if os.path.exists(file_path):
            df = pd.read_parquet(file_path)
            df['simulation_index'] = int(dir_name)
            dataframes.append(df)

    return pd.concat(dataframes, ignore_index=True)

def params_wake(source_path: str, id: int) -> pd.DataFrame:
    wake_params = pd.read_parquet(os.path.join(source_path, "wakes", str(id), "param.parquet"))
    wake_df = pd.read_parquet(os.path.join(source_path, "wakes", str(id), "wakes_df.parquet"))
    return wake_params, wake_df

def read_encounter_results(source_path: str) -> pd.DataFrame:
    directories = sorted([d for d in os.listdir(os.path.join(source_path, "encounters")) if d.isdigit()])
    results = {}

    for dir_name in directories:
        results_file_path = os.path.join(source_path, "encounters", dir_name, "results.parquet")
        encounter_file_path = os.path.join(source_path, "encounters", dir_name, "encounter_df.parquet")
        
        if os.path.exists(results_file_path):
            df = pd.read_parquet(results_file_path)
            enc = pd.read_parquet(encounter_file_path)
            df.index = df.index + enc.index[0] + 1 # +1 because the results are calculated for t=1 and not t=0
            results[int(dir_name)] = df
    return results