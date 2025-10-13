#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 28 10:19:21 2025

@author: g
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import glob
import time
import re
import matplotlib.animation as animation
import pyvista as pv
from IO import extract_number
import pdb

# Load CSV file
def load_data(file_path):
    data = pd.read_csv(file_path, header=None)
    theta, phi = data[0], data[1]  # Convert to radians
    base_height, flow_height = data[2], data[3]
    vel_x, vel_y = data[4], data[5]
    return theta, phi, base_height, flow_height, vel_x, vel_y

def get_color_limits_vel(folder_path):
    file_paths = sorted(glob.glob(f"{folder_path}/field_*.csv"))
    min_val, max_val = float('inf'), float('-inf')
    for file_path in file_paths:
        _, _, _, flow_height, vel_x, vel_y = load_data(file_path)
        min_val = min(min_val, vel_x.min())
        max_val = max(max_val, vel_x.max())
    return min_val, max_val


def get_color_limits_height(folder_path,epsilon):
    file_paths = sorted(glob.glob(f"{folder_path}/field_*.csv"))
    min_val, max_val = float('inf'), float('-inf')
    i=0
    for file_path in file_paths:
        _, _, base, flow_height, vel_x, _ = load_data(file_path)
        flow_height =(flow_height +base)*epsilon[i]*250
        i=i+1
        min_val = min(min_val, flow_height.min())
        max_val = max(max_val, flow_height.max())
    return min_val, max_val

# Plot height as contour
def plot_height(theta, phi, flow_height, time_step,vmin,vmax):
    plt.rcParams.update({'font.size' : 14})
    custom_levels = np.arange(-1.0, 1.0, 0.15)
    sc = plt.tricontourf(phi*180/np.pi, theta*180/np.pi, flow_height, levels = custom_levels, cmap='viridis')#,vmin=-4,vmax=4)
    plt.colorbar(sc, label='Height')
    plt.xlabel('Phi (degree)')
    plt.ylabel('Theta (degree)')
    #plt.title(f'Flow Height Contour (Time Step {time_step})')
    plt.draw()
    plt.pause(0.5)
    #plt.axis('equal')
    plt.minorticks_on()
    plt.tick_params(direction='in',right=True, top=True, left=True, bottom=True)
    plt.tick_params(labelsize=14)
    plt.tick_params(labelbottom=True, labeltop=False, labelright=False, labelleft=True)
    plt.title(f"{time_step+1} impacts")
    #plt.tick_params(direction='in',which='minor', length=5, bottom=True, top=True, left=True, right=True)
    #plt.tick_params(direction='in',which='major', length=10, bottom=True, top=True, left=True, right=True)
    #plt.legend(loc='best')
    #plt.legend(fontsize=14) 
    plt.savefig(f'/home/g/Asteroids/output/Spherical_2D_3.5/img/{time_step}.png',dpi = 300)

# Plot velocity as quiver
def plot_velocity(theta, phi, vel_x, vel_y, time_step):
    plt.figure(figsize=(8, 6))
    plt.quiver(theta, phi, vel_x, vel_y, scale=10, color='red')
    plt.xlabel('Theta (radians)')
    plt.ylabel('Phi (radians)')
    plt.title(f'Velocity Field (Time Step {time_step})')
    plt.show()
    
# Initialize figure
def init_plot():
    fig = plt.figure(figsize=(8, 6))
    ax = fig.add_axes([0.1, 0.1, 0.70, 0.8])  # Manually set plot area (left, bottom, width, height)
    
    cbar_ax = fig.add_axes([0.85, 0.15, 0.03, 0.7])  # Colorbar axis
    return fig, ax, cbar_ax

def update(frame, file_paths, ax, cbar_ax, vmin, vmax):
    ax.clear()
    file_path = file_paths[frame]
    theta, phi, _, flow_height, vel_x, vel_y = load_data(file_path)
    sc = ax.tricontourf(theta*180/np.pi, phi*180/np.pi, flow_height, cmap='viridis', vmin=vmin, vmax=vmax)
    #ax.quiver(theta, phi, vel_x, vel_y, scale=10, color='red')
    ax.set_xlabel('Theta (degree)')
    ax.set_ylabel('Phi (degree)')
    ax.set_title(f'Flow Visualization (Time Step {frame})')
    cbar_ax.clear()
    plt.colorbar(sc, cax=cbar_ax, label='Height')
    return sc

# Main function to execute visualization and save animation
def main(folder_path, output_file="animation.mp4"):

    
    vmin, vmax = get_color_limits_height(folder_path)  # Get global color limits
    file_paths = sorted(glob.glob(f"{folder_path}/field_*.csv"),key=extract_number)
    fig, ax, cbar_ax = init_plot()
    ani = animation.FuncAnimation(fig, update, frames=len(file_paths), fargs=(file_paths, ax, cbar_ax, vmin, vmax), interval=500)
    ani.save(output_file, writer='ffmpeg', fps=2)
    plt.show()


# Main function to execute visualization
def plot_height_main(folder_path,epsilon,Gamma):  
    
    vmin, vmax = get_color_limits_height(folder_path)  # Get global color limits
    file_paths = sorted(glob.glob(f"{folder_path}/field_*.csv"),key=extract_number)
   # print(file_paths)
    plt.figure(figsize=(8, 6))
    plt.ion()  # Turn on interactive mode
    for i, file_path in enumerate(file_paths):
        print(f"Processing: {file_path}")
        theta, phi, base_height, flow_height, vel_x, vel_y = load_data(file_path)
        
        plt.clf()

        plot_height(theta, phi, (flow_height*epsilon)*250, i,vmin,vmax)
        #plot_velocity(theta, phi, vel_x, vel_y, i)
        time.sleep(1)  # Add delay to observe each step
        print(vmin,vmax)


# Main function to execute visualization
def plot_surface_main(folder_path):
    
    #pdb.set_trace()
    file1 = "/home/g/Asteroids/output/Spherical_2D_3.5/run1/data/dia.txt"
    epsilon=np.loadtxt(file1,dtype=float)[:,2]
    Gamma = np.loadtxt(file1,dtype=float)[:,3]
    vmin, vmax = (0,0) #get_color_limits_height(folder_path,epsilon)  # Get global color limits
    file_paths = sorted(glob.glob(f"{folder_path}/field_*.csv"),key=extract_number)
   # print(file_paths)
    plt.figure(figsize=(8, 6))
    plt.ion()  # Turn on interactive mode
    for i, file_path in enumerate(file_paths):
        print(f"Processing: {file_path}")

        theta, phi, base_height, flow_height, vel_x, vel_y = load_data(file_path)
        
        plt.clf()

        plot_height(theta, phi, (Gamma[i]*base_height+epsilon[i]*flow_height[i])*250, i,vmin,vmax)
        #plot_velocity(theta, phi, vel_x, vel_y, i)
        time.sleep(1)  # Add delay to observe each step
        print(vmin,vmax)

        
def plot_vel_main(folder_path):
    
    vmin, vmax = get_color_limits_vel(folder_path)  # Get global color limits
    file_paths = sorted(glob.glob(f"{folder_path}/field_*.csv"),key=extract_number)
    plt.figure(figsize=(8, 6))
    plt.ion()  # Turn on interactive mode
    for i, file_path in enumerate(file_paths):
        print(f"Processing: {file_path}")
        theta, phi, base_height, flow_height, vel_x, vel_y = load_data(file_path)
        
        plt.clf()
        plot_height(theta, phi, vel_x, i,vmin,vmax)

        #plot_velocity(theta, phi, vel_x, vel_y, i)
        time.sleep(1)  # Add delay to observe each step

        
        
# Example usage
# file1 = "/home/g/Asteroids/output/Crater/run1/data/dia.txt"
# epsilon=np.loadtxt(file1,dtype=float)[:,2]
# Gamma = np.loadtxt(file1,dtype=float)[:,3]
# for i in range(1,2):
#   #  plot_height_main(f"/home/g/Asteroids/output/Crater/run1/data/landslides_{i}",epsilon[i],Gamma[i])  # Uncomment and replace with your folder path
#     plt.close()

# plot_surface_main("/home/g/Asteroids/output/crater_3/run1/data")  # Uncomment and replace with your folder path




def generate_asteroid_shape(reference_radius, lats_rad, lons_rad, deviations):
    """
    Generates mesh vertices and faces from spherical coordinates and deviations.

    Args:
        reference_radius (float): The radius of the base sphere.
        lats_rad (np.ndarray): 2D array of latitudes in radians (shape: [n_lat, n_lon]).
        lons_rad (np.ndarray): 2D array of longitudes in radians (shape: [n_lat, n_lon]).
        deviations (np.ndarray): 2D array of height deviations from the reference
                                 radius (shape: [n_lat, n_lon]).

    Returns:
        tuple: (vertices, faces) suitable for creating a PyVista mesh,
               or (None, None) on error.
               vertices (np.ndarray): Array of vertex coordinates (N, 3).
               faces (np.ndarray): Array defining triangular faces.
    """
    if not (lats_rad.shape == lons_rad.shape == deviations.shape):
        print("Error: Input array shapes (latitude, longitude, deviations) must match.")
        return None, None

    n_lat, n_lon = lats_rad.shape
    print(f"Generating shape for grid: {n_lat} latitudes, {n_lon} longitudes.")

    # --- 1. Calculate actual radius for each point ---
    actual_radius = reference_radius + deviations
    # Ensure radius doesn't go below zero (or a very small positive number)
    actual_radius = np.maximum(actual_radius, 1e-6)

    # --- 2. Convert spherical coordinates (with actual radius) to Cartesian ---
    x = actual_radius * np.sin(lats_rad) * np.cos(lons_rad)
    y = actual_radius * np.sin(lats_rad) * np.sin(lons_rad)
    z = actual_radius * np.cos(lats_rad)

    # --- 3. Create a PyVista StructuredGrid ---
    # PyVista expects X, Y, Z arrays defining the grid structure.
    # Note: The dimensions should match the original grid.
    grid = pv.StructuredGrid(x, y, z)

    # --- 4. Extract the outer surface mesh (automatically generates faces) ---
    # This is the easiest way to get a topologically correct surface
    # It handles poles and longitude wrap-around implicitly.
    try:
        surface = grid.extract_surface(nonlinear_subdivision=2) # Subdivision can smooth near poles
        # Recalculate normals for better shading
        surface.compute_normals(cell_normals=False, point_normals=True, inplace=True)
    except Exception as e:
        print(f"Error extracting surface from grid: {e}")
        # Fallback if extract_surface fails (less robust)
        print("Attempting simple triangulation...")
        try:
             surface = grid.delaunay_3d().extract_surface() # Alternative approach
             surface.compute_normals(cell_normals=False, point_normals=True, inplace=True)
        except Exception as e2:
             print(f"Fallback triangulation also failed: {e2}")
             return None, None


    print(f"Generated surface mesh with {surface.n_points} vertices and {surface.n_cells} faces.")
    return surface.points, surface.faces # Return vertices and faces

# --- Main Execution ---
def visualization_3D(R, lats_rad, lons_rad, dev):



    # --- Generate Vertices and Faces ---
    vertices, faces = generate_asteroid_shape(R, lats_rad, lons_rad, dev)

    if vertices is not None and faces is not None:
        # --- Create PyVista Mesh Object ---
        asteroid_polydata = pv.PolyData(vertices, faces=faces)

        # Add the deviation data itself as scalars to the mesh for coloring
        # Need to reshape the deviation data to match the number of vertices
        # Note: This simple flatten assumes vertex order matches grid flatten order
        try:
            asteroid_polydata['Deviation'] = dev.flatten(order='C') # Try Fortran order first (often matches StructuredGrid)
            if len(asteroid_polydata['Deviation']) != asteroid_polydata.n_points:
                 print("Warning: Deviation array length mismatch, trying 'F' order flatten.")
                 asteroid_polydata['Deviation'] = dev.flatten(order='F') # Try C order
            if len(asteroid_polydata['Deviation']) != asteroid_polydata.n_points:
                 print("Warning: Could not map deviation data to vertices for coloring.")
                 if 'Deviation' in asteroid_polydata.point_data: del asteroid_polydata['Deviation'] # Remove if incorrect
        except Exception as e:
             print(f"Error assigning deviation scalars: {e}")


        # --- Plotting ---
        print("Plotting the reconstructed shape...")
        plotter = pv.Plotter(window_size=[800, 800])

        # Add the mesh. Color by deviation if available, otherwise use a solid color.
        if 'Deviation' in asteroid_polydata.point_data:
            plotter.add_mesh(asteroid_polydata, scalars='Deviation', cmap='viridis',
                             scalar_bar_args={'title': 'Deviation from Ref. Radius'},
                             show_edges=False) # show_edges=True can help visualize grid
        else:
            plotter.add_mesh(asteroid_polydata, color='tan', show_edges=False)

        plotter.show_axes()
        plotter.add_text("Reconstructed Asteroid Shape", position='upper_edge', font_size=10)
        print("Showing plotter window...")
        plotter.show()
        print("Plotter closed.")

        # --- Optional: Save the generated mesh ---
        # output_filename = "reconstructed_asteroid.obj" # Or .ply, .stl, .vtk
        # print(f"Saving mesh to {output_filename}...")
        # asteroid_polydata.save(output_filename)
        # print("Mesh saved.")

    else:
        print("Failed to generate asteroid shape.")