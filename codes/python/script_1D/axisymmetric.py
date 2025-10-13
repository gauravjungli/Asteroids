#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 18 11:13:39 2023

@author: kumargaurav
"""

import numpy as np
import matplotlib.pyplot as plt
from IO import Output_File, load_asteroid_data

# --- Helper function to plot the (r, theta) profile ---
def plot_polar_profile(theta, r_values, title="Shape Profile"):
    """
    Plots the shape defined by r(theta) by converting to Cartesian coordinates.
    Assumes theta is in radians.
    """
    x = r_values * np.sin(theta)
    y = r_values * np.cos(theta)

    plt.figure(figsize=(8, 6))
    plt.plot(x, y, label=f'Profile: {title}')
    
    
    plt.gca().set_aspect('equal', adjustable='box')
    
    plt.title(title)
    plt.xlabel("x (Axis of Revolution)")
    plt.ylabel("y (Profile Radius)")
    plt.grid(True, linestyle=':', alpha=0.7)
    plt.legend()
    plt.xlim(left=0) 
    plt.show()

# --- Shape Generation Functions (same as before) ---

def get_oblate_spheroid_profile(theta, polar_radius_a, equatorial_radius_b):
    if polar_radius_a <= 0 or equatorial_radius_b <= 0:
        raise ValueError("Radii must be positive.")
    if polar_radius_a >= equatorial_radius_b:
        print("Warning: For oblate spheroid, polar_radius_a should be < equatorial_radius_b.")
    
    num = polar_radius_a * equatorial_radius_b
    den = np.sqrt((equatorial_radius_b * np.cos(theta))**2 + \
                  (polar_radius_a * np.sin(theta))**2)
    return num / (den + 1e-12) # Added small epsilon for strict positivity

def get_prolate_spheroid_profile(theta, polar_radius_a, equatorial_radius_b):
    if polar_radius_a <= 0 or equatorial_radius_b <= 0:
        raise ValueError("Radii must be positive.")
    if polar_radius_a <= equatorial_radius_b:
        print("Warning: For prolate spheroid, polar_radius_a should be > equatorial_radius_b.")
    
    return get_oblate_spheroid_profile(theta, polar_radius_a, equatorial_radius_b)

def get_sinusoidal_profile(theta, base_radius, amplitude, frequency, phase=0):
    if base_radius <= abs(amplitude):
        print(f"Warning: For sinusoidal profile, base_radius ({base_radius}) "
              f"should be > abs(amplitude) ({abs(amplitude)}) to ensure r > 0 strictly.")
    
    r = base_radius + amplitude * np.sin(frequency * theta + phase)
    return np.maximum(r, 1e-9) # Ensure r is practically positive

def get_gaussian_profile(theta, base_radius, amplitude, peak_angle, width_sigma):
    if base_radius < 0:
        print("Warning: base_radius for Gaussian profile should ideally be non-negative.")
    if base_radius == 0 and amplitude <= 0:
        print("Warning: If base_radius is 0, amplitude must be positive for r > 0.")
    
    exponent = -((theta - peak_angle)**2) / (2 * width_sigma**2)
    r = base_radius + amplitude * np.exp(exponent)
    return np.maximum(r, 1e-9) # Ensure r is practically positive

# --- Function to calculate derivatives and save data ---
def save_profile_data(filename, theta_vals, r_vals):
    """
    Calculates dr/dtheta, d^2r/dtheta^2 and saves the data to a CSV file.
    Columns: theta, r, 1, dr_dtheta, d2r_dtheta2
    """
    if len(theta_vals) < 3:
        print("Warning: Not enough points to accurately calculate second derivatives.")
        dr_dtheta = np.full_like(theta_vals, np.nan)
        d2r_dtheta2 = np.full_like(theta_vals, np.nan)
    else:
        # Calculate first derivative: dr/dtheta
        dr_dtheta = np.gradient(r_vals, theta_vals, edge_order=2)
        
        # Calculate second derivative: d^2r/dtheta^2
        d2r_dtheta2 = np.gradient(dr_dtheta, theta_vals, edge_order=2)

    # Create the column of ones
    ones_column = np.ones_like(theta_vals)

    # Stack the data columns: theta, r, 1, dr/dtheta, d^2r/dtheta^2
    data_to_save = np.stack((theta_vals, r_vals, ones_column, dr_dtheta, d2r_dtheta2), axis=-1)


    # Save to file
    np.savetxt(filename, data_to_save, delimiter=",", comments='')
    print(f"Data saved to {filename}")



# --- Main part: Generate shapes and optionally save data ---
def axisymmetric(base,a,c,R):
    # Define the theta range (0 to pi)
    
    theta_vals = base[:,0]

    # 1. Oblate Spheroid Profile
    r_oblate = get_oblate_spheroid_profile(theta_vals, polar_radius_a=c, equatorial_radius_b=a)
 

    # # 2. Prolate Spheroid Profile
    # r_prolate = get_prolate_spheroid_profile(theta_input_vals, polar_radius_a=2.5, equatorial_radius_b=1.0)
    # plot_polar_profile(theta_input_vals, r_prolate, "Prolate Spheroid Profile (a=2.5, b=1.0)")
    
    # # Sphere Profile
    # r_sphere = get_oblate_spheroid_profile(theta_input_vals, polar_radius_a=1.5, equatorial_radius_b=1.5)
    # plot_polar_profile(theta_input_vals, r_sphere, "Sphere Profile (a=1.5, b=1.5)")

    # # 3. Sinusoidal Profile
    # r_sinusoidal = get_sinusoidal_profile(theta_vals, 
    #                                       base_radius=R, 
    #                                       amplitude=0.3*R, 
    #                                       frequency=6,
    #                                       phase=np.pi/2)
    # plot_polar_profile(theta_vals, r_sinusoidal)


    # # 4. Gaussian-like Profile
    # r_gaussian1 = get_gaussian_profile(theta_input_vals,
    #                                    base_radius=0.2,
    #                                    amplitude=1.0,
    #                                    peak_angle=np.pi/2,
    #                                    width_sigma=0.5)
    # plot_polar_profile(theta_input_vals, r_gaussian1, "Gaussian Bump Profile (centered)")

    # r_gaussian_dip = get_gaussian_profile(theta_input_vals,
    #                                       base_radius=1.5,
    #                                       amplitude=-0.7,
    #                                       peak_angle=np.pi/2,
    #                                       width_sigma=0.4)
    # plot_polar_profile(theta_input_vals, r_gaussian_dip, "Gaussian Dip Profile")

    # --- Choose a profile and save its data to base.txt ---
    
    # MODIFY THE LINE BELOW TO CHOOSE WHICH 'r_values' TO SAVE
    # For example, to save the prolate spheroid data:
    r_vals = r_oblate

    # Or, to save the Gaussian bump data:
    # chosen_r_values_for_saving = r_gaussian1
    # profile_description = "Gaussian Bump Profile (centered)"
    # Or, to save the sinusoidal data:
    # chosen_r_values_for_saving = r_sinusoidal1
    # profile_description = "Sinusoidal Profile (freq=6)"
    if len(theta_vals) < 3:
        print("Warning: Not enough points to accurately calculate second derivatives.")
        dr_dtheta = np.full_like(theta_vals, np.nan)
        d2r_dtheta2 = np.full_like(theta_vals, np.nan)
    else:
        # Calculate first derivative: dr/dtheta
        dr_dtheta = np.gradient(r_vals, theta_vals, edge_order=2)
        
        # Calculate second derivative: d^2r/dtheta^2
        d2r_dtheta2 = np.gradient(dr_dtheta, theta_vals, edge_order=2)
    base[:,1] = r_vals/R
    base[:,3] = dr_dtheta/R
    base[:,4] = d2r_dtheta2/R
    
    

    
