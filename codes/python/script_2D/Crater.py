#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Mar 29 06:00:08 2025

@author: g
"""

import numpy as np
from IO import Output_File
import matplotlib.pyplot as plt
# Import color normalization tools
import matplotlib.colors as mcolors
from Initialize import grid_2D
# Optional: for 3D plotting if you want to visualize
# from mpl_toolkits.mplot3d import Axes3D # Already implicitly imported by projection='3d'

# --- Helper Functions (keep as before) ---

def spherical_to_cartesian(rho, theta, phi):
    """
    Converts spherical coordinates (radius, inclination, azimuth)
    to Cartesian coordinates (x, y, z).
    Assumes theta is inclination (polar angle, from +Z) in [0, pi]
    Assumes phi is azimuth (from +X towards +Y) in [0, 2*pi]
    rho, theta, phi can be arrays.
    """
    x = rho * np.sin(theta) * np.cos(phi)
    y = rho * np.sin(theta) * np.sin(phi)
    z = rho * np.cos(theta)
    return x, y, z

def calculate_angular_distance(theta1, phi1, theta2, phi2):
    """
    Calculates the central angle (angular distance) between two points
    on a sphere given their spherical coordinates (theta, phi).
    Uses the formula derived from the dot product of Cartesian vectors.
    Assumes angles are in radians.
    """
    # Convert to unit Cartesian vectors
    x1 = np.sin(theta1) * np.cos(phi1)
    y1 = np.sin(theta1) * np.sin(phi1)
    z1 = np.cos(theta1)

    x2 = np.sin(theta2) * np.cos(phi2)
    y2 = np.sin(theta2) * np.sin(phi2)
    z2 = np.cos(theta2)

    # Calculate dot product - can be done element-wise if theta1/phi1 are arrays
    dot_product = x1*x2 + y1*y2 + z1*z2

    # Clamp dot product to [-1, 1] due to potential floating point errors
    dot_product = np.clip(dot_product, -1.0, 1.0)

    # Calculate angle
    angle = np.arccos(dot_product)
    return angle

# --- Crater Function (keep as before) ---

def add_gaussian_crater_to_grid(
    rho_grid, theta_grid, phi_grid,
    crater_theta, crater_phi,
    crater_radius_angle, max_depth,
    rim_height=0.0, rim_width_factor=1.5
    ):
    """
    Modifies rho_grid to add a Gaussian-shaped crater depression and optional rim.

    Args:
        rho_grid (np.ndarray): 2D array of current radial distances. Modified in-place or returned.
        theta_grid (np.ndarray): 2D array of theta coordinates (inclination).
        phi_grid (np.ndarray): 2D array of phi coordinates (azimuth).
        crater_theta (float): Inclination angle (radians) of the crater center.
        crater_phi (float): Azimuthal angle (radians) of the crater center.
        crater_radius_angle (float): Angular radius (radians) defining the crater width (related to Gaussian sigma).
        max_depth (float): Maximum depth of the crater depression (positive value, subtracted from radius).
        rim_height (float): Maximum height of the raised rim (positive value, added to radius). Defaults to 0.0.
        rim_width_factor (float): How much wider the rim effect extends compared to the
                                  depression radius (e.g., 1.5 means rim extends to
                                  1.5 * crater_radius_angle). Defaults to 1.5.

    Returns:
        np.ndarray: The modified rho_grid.
    """
    if not (rho_grid.shape == theta_grid.shape == phi_grid.shape):
        raise ValueError("Input grid shapes must match")

    #print(f"Adding crater at theta={np.degrees(crater_theta):.1f} deg, phi={np.degrees(crater_phi):.1f} deg")

    # Calculate angular distance from EACH grid point to the crater center
    angular_distances = calculate_angular_distance(theta_grid, phi_grid, crater_theta, crater_phi)

    # --- Calculate Depression ---
    sigma_depression = crater_radius_angle 
    depression = max_depth * np.exp(-angular_distances**2 / (2 * sigma_depression**2))


    # --- >>> Calculate the 2-sigma mask <<< ---
    threshold_angle = 2.0 * sigma_depression
    mask_within_2sigma = angular_distances <= threshold_angle
    num_marked = np.sum(mask_within_2sigma)
 #   print(f"  Marking {num_marked} points within 2-sigma angular distance ({np.degrees(threshold_angle):.2f} deg).")
    # --- End of mask calculation ---


    # --- Calculate Rim (Optional) ---
    rim_displacement = np.zeros_like(rho_grid)
    if rim_height > 0:
        rim_outer_angle = crater_radius_angle * rim_width_factor
        sigma_rim = (rim_outer_angle - crater_radius_angle) / 2.5
        if sigma_rim < 1e-6:
             print("Warning: Rim width factor too small, rim sigma near zero.")
             sigma_rim = 1e-6

        peak_angle = crater_radius_angle
        rim_raw = rim_height * np.exp(-(angular_distances - peak_angle)**2 / (2 * sigma_rim**2))

        taper_mask = angular_distances > peak_angle
        taper_factor = np.clip(1.0 - (angular_distances[taper_mask] - peak_angle) / (rim_outer_angle - peak_angle + 1e-9), 0.0, 1.0)
        rim_displacement = rim_raw
        rim_displacement[taper_mask] *= taper_factor
        rim_displacement[angular_distances > rim_outer_angle] = 0.0

    modified_rho_grid = rho_grid - depression + rim_displacement
    # modified_rho_grid = np.maximum(modified_rho_grid, 0.0) # Optional floor
  #  print("Crater calculation complete.")
    return modified_rho_grid, mask_within_2sigma

def crater_radius_to_angle(crater_radius, sphere_radius):
    """
    Converts a physical crater radius (measured along the surface)
    to an angular radius in radians.

    Args:
        crater_radius (float): Physical radius of the crater (e.g., in km).
        sphere_radius (float): Radius of the base sphere (in the SAME units as crater_radius).

    Returns:
        float: Crater radius in radians.
    """
    if sphere_radius <= 0:
        raise ValueError("Sphere radius must be positive.")
    if crater_radius < 0:
        raise ValueError("Crater radius cannot be negative.")
    # Ensure calculation is done in float
    return float(crater_radius) / float(sphere_radius)


def Crater(parameters,target,impactor):
    
    # 1. Define the Spherical Grid
    initial_sphere_radius = target.d/2 
    epsilon =float(parameters['epsilon'])
    mydir=Output_File (parameters,"output",["base.txt"]) 
    
    base=np.loadtxt(mydir,delimiter=",",dtype=float)
    n_theta = int(parameters['X Resolution'])
    n_phi = int(parameters['Y Resolution'])
    
    theta, phi = grid_2D(parameters)
    phi_grid, theta_grid = np.meshgrid(phi, theta)
    
    
    expected_rows = n_theta * n_phi
    if base.shape[0] != expected_rows:
        raise ValueError(f"Input data has {base.shape[0]} rows, but expected {expected_rows} based on n_theta={n_theta}, n_phi={n_phi}")
    if base.shape[1] != 7:
        raise ValueError(f"Input data must have 7 columns (theta, phi, rho, height), but found {base.shape[1]}")



    # 3. <<< --- REPLACEMENT for np.full_like --- >>>
    # Extract rho values (3rd column, index 2)
    rho_values_flat = base[:, 2]

    # Reshape rho values into the grid, assuming C-style (row-major) order
    # This order corresponds to iterating through theta (rows) then phi (columns)
    # when the input data was created/flattened.
    rho_grid_initial = rho_values_flat.reshape((n_theta, n_phi))
    rho_grid_initial = initial_sphere_radius*(1+epsilon*rho_grid_initial)
  #  print(f"Reshaped initial rho_grid shape: {rho_grid_initial.shape}")
    

    # --- Define Crater Parameters (using physical radius) ---
    crater_theta = impactor.Theta
    crater_phi = impactor.Phi
    # --- > Define physical radius FIRST < ---
    crater_physical_radius = (impactor.M/target.dens)**(1/3)*0.62*(target.grav*impactor.d/2/impactor.vel**2)**(-0.17)# 15*impactor.d # In the same 'units' as initial_sphere_radius
    # --- > CONVERT physical radius to angular radius < ---
    crater_radius_angle = crater_radius_to_angle(crater_physical_radius, initial_sphere_radius)
    # --- > Use the calculated angle below < ---
    print(f"Crater : Physical Radius={crater_physical_radius} -> Angular Radius={np.degrees(crater_radius_angle):.2f} degrees")

    crater_depth = 0.3*crater_physical_radius 
    crater_rim_h = 0.0 # No rim
    crater_rim_w = 1.5

    # --- Add Craters (using the calculated angular radius) ---
    rho_grid_mod, mask = add_gaussian_crater_to_grid(
        rho_grid_initial,
        theta_grid, phi_grid,
        crater_theta, crater_phi, crater_radius_angle, crater_depth, # Pass the angle
        rim_height=crater_rim_h, rim_width_factor=crater_rim_w
    )


    # 1. Flatten the final rho grid and the corresponding theta/phi grids
    # Ensure the flattening order is consistent (default 'C' or row-major)
    theta_flat = theta_grid.ravel()
    phi_flat = phi_grid.ravel()
    rho_mod_flat = (rho_grid_mod.ravel()/initial_sphere_radius -1 )
    mask_flat_int = mask.ravel().astype(int)
    # 2. Stack them as columns [theta, phi, rho_modified]
    # This creates an array of shape (n_theta * n_phi, 3)
    output_spherical_data = np.column_stack((theta_flat, phi_flat, rho_mod_flat, mask_flat_int))

   # print(f"Output data shape: {output_spherical_data.shape}")



    # --- Visualization (keep as before, with colorbar fix) ---
    # print("\nPreparing visualization...")
    # X, Y, Z = spherical_to_cartesian(rho_grid_mod, theta_grid, phi_grid)
    # height_deviation = rho_grid_mod - initial_sphere_radius
    # max_abs_dev = np.max(np.abs(height_deviation))
    # if max_abs_dev < 1e-9: max_abs_dev = 1e-9
    # norm = mcolors.Normalize(vmin=-max_abs_dev, vmax=max_abs_dev)
    # cmap = plt.get_cmap('coolwarm')
    # face_colors = cmap(norm(height_deviation))

    # fig = plt.figure(figsize=(10, 8))
    # ax = fig.add_subplot(111, projection='3d')
    # surf = ax.plot_surface(X, Y, Z, facecolors=face_colors, rstride=1, cstride=1, linewidth=0, antialiased=False, shade=True)
    # mappable = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    # mappable.set_array(height_deviation)
    # cbar = fig.colorbar(mappable, ax=ax, shrink=0.6, aspect=20, pad=0.1) # Corrected line
    # cbar.set_label('Height Deviation from Base Sphere')

    # max_range = np.array([X.max()-X.min(), Y.max()-Y.min(), Z.max()-Z.min()]).max() / 2.0
    # mid_x = (X.max()+X.min()) * 0.5
    # mid_y = (Y.max()+Y.min()) * 0.5
    # mid_z = (Z.max()+Z.min()) * 0.5
    # ax.set_xlim(mid_x - max_range, mid_x + max_range)
    # ax.set_ylim(mid_y - max_range, mid_y + max_range)
    # ax.set_zlim(mid_z - max_range, mid_z + max_range)
    # ax.set_xlabel("X axis")
    # ax.set_ylabel("Y axis")
    # ax.set_zlabel("Z axis")
    # ax.set_title("Sphere with Craters (Colored by Height Deviation)")

    # print("Showing plot...")
    # plt.show()
    # print("Plot window closed.")
    
    ############################################################################################
    
    
    return output_spherical_data, crater_depth