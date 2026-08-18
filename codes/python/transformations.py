#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 11 12:04:53 2026

@author: g
"""
import numpy as np
    
def spherical_to_cartesian(theta, phi, r=1):
    x = r * np.sin(theta) * np.cos(phi)
    y = r * np.sin(theta) * np.sin(phi)
    z = r * np.cos(theta)
    return np.stack([x, y, z], axis=-1)

def cartesian_to_spherical(cartesian_coords):
    x, y, z = cartesian_coords[..., 0], cartesian_coords[..., 1], cartesian_coords[..., 2]
    r = np.sqrt(x**2 + y**2 + z**2)
    theta = np.arccos(np.clip(z / r, -1.0, 1.0))  # Ensuring theta is in [0, pi]
    phi = np.arctan2(y, x)
    return np.stack([theta, phi], axis=-1)

def rotation_matrix(theta_0, phi_0):
    Rz = np.array([[np.cos(phi_0), -np.sin(phi_0), 0],
                   [np.sin(phi_0), np.cos(phi_0), 0],
                   [0, 0, 1]])
    
    Ry = np.array([[np.cos(theta_0), 0, np.sin(theta_0)],
                   [0, 1, 0],
                   [-np.sin(theta_0), 0, np.cos(theta_0)]])
    
    return Ry @ Rz  # Apply rotation about z-axis first, then y-axis

def rotate_spherical_coordinates(theta, phi, theta_0, phi_0, r=1):
    cartesian_coords = spherical_to_cartesian(theta, phi, r)
    R = rotation_matrix(theta_0, phi_0)
    rotated_cartesian = np.einsum('ij,...j->...i', R, cartesian_coords)
    return cartesian_to_spherical(rotated_cartesian)


def transform_spherical_coords(theta, phi, theta_0, phi_0):
    """
    Transforms spherical coordinates (theta, phi) to a new system (theta', phi')
    where the direction (theta_0, phi_0) becomes the new North Pole (Z' axis).

    Assumes standard physics spherical coordinates:
    - theta: polar angle from Z+ (radians, 0 to pi)
    - phi: azimuthal angle from X+ (radians, 0 to 2*pi)

    Args:
        theta (float or np.ndarray): Original polar angle(s) in radians.
        phi (float or np.ndarray): Original azimuthal angle(s) in radians.
        theta_0 (float): Polar angle of the new North Pole in the old system (radians).
        phi_0 (float): Azimuthal angle of the new North Pole in the old system (radians).

    Returns:
        tuple[float or np.ndarray, float or np.ndarray]:
            - theta_prime: Polar angle(s) in the new coordinate system (radians).
            - phi_prime: Azimuthal angle(s) in the new coordinate system (radians, 0 to 2*pi).
    """
    # Ensure inputs are numpy arrays for broadcasting
    theta = np.asarray(theta)
    phi = np.asarray(phi)

    # === Step 1: Define the New Z' axis ===
    sin_t0 = np.sin(theta_0)
    cos_t0 = np.cos(theta_0)
    sin_p0 = np.sin(phi_0)
    cos_p0 = np.cos(phi_0)

    Z_prime_vec = np.array([sin_t0 * cos_p0, sin_t0 * sin_p0, cos_t0])

    # === Step 2: Define the New X' and Y' axes ===
    # Choose X' to point "downhill" from Z' towards the original Z axis, projected.
    # Or equivalently, Z = cos(theta_0)*Z' + sin(theta_0)*X' (if phi_0=0)
    # So X' is proportional to Z - cos(theta_0)*Z'
    # Unit vector for original Z
    Z_vec = np.array([0.0, 0.0, 1.0])
    # Component of Z along Z'
    Z_proj_Zprime = cos_t0 * Z_prime_vec
    # Vector pointing from Z' towards Z (perpendicular to Z')
    X_prime_dir_temp = Z_vec - Z_proj_Zprime

    # Handle edge case where new pole is old pole (theta_0 = 0 or pi)
    # Use small epsilon to avoid division by zero if theta_0 is exactly 0 or pi
    epsilon = 1e-12
    norm_X_prime_temp = np.linalg.norm(X_prime_dir_temp)

    if norm_X_prime_temp < epsilon:
         # If Z' aligns with Z, the standard X,Y axes (rotated by phi_0) can be used.
         # X' should correspond to theta'=pi/2, phi'=0 in new system.
         # A simple choice is rotation of original X by phi_0: [cos(phi_0), sin(phi_0), 0]
         # Let's refine this using Y' first for this case.
         # Y' = Z' x X'. If Z' = Z, Y' should be related to original Y rotated by phi_0.
         # Let's use the standard cross product derived Y' vector below, it should simplify.
         Y_prime_vec = np.array([sin_p0, -cos_p0, 0.0])
         # Then X' = Y' x Z' (note order for right-handed system: X=YxZ)
         X_prime_vec = np.cross(Y_prime_vec, Z_prime_vec)
    else:
        # Normalize the derived X' direction
        X_prime_vec = X_prime_dir_temp / norm_X_prime_temp
        # Calculate Y' = Z' x X' to complete the right-handed system
        Y_prime_vec = np.cross(Z_prime_vec, X_prime_vec)


    # === Step 3: Convert Original Point(s) to Cartesian ===
    # Assume unit sphere (rho=1) as only direction matters for angles
    sin_t = np.sin(theta)
    cos_t = np.cos(theta)
    sin_p = np.sin(phi)
    cos_p = np.cos(phi)

    # Need to handle theta possibly being an array
    if theta.shape: # If theta is an array
        P_vec_x = sin_t * cos_p
        P_vec_y = sin_t * sin_p
        P_vec_z = cos_t
        # Stack into (N, 3) array if theta/phi are 1D, or handle higher dims if needed
        # For dot product using np.dot, easier if P_vec is (3, N) or similar
        P_vec = np.stack((P_vec_x, P_vec_y, P_vec_z), axis=0) # Shape (3, ...)
    else: # If theta is a scalar
         P_vec = np.array([sin_t * cos_p, sin_t * sin_p, cos_t]) # Shape (3,)


    # === Step 4: Project onto New Axes ===
    # Use np.dot which handles matrix/vector multiplication correctly
    # If P_vec is (3,) -> result is scalar
    # If P_vec is (3, N) -> result is (N,)
    x_prime = np.dot(X_prime_vec, P_vec)
    y_prime = np.dot(Y_prime_vec, P_vec)
    z_prime = np.dot(Z_prime_vec, P_vec)

    # === Step 5: Convert New Cartesian to New Spherical ===
    # Calculate new polar angle theta'
    # Clip argument to arccos to avoid domain errors due to floating point inaccuracies
    rho_prime = np.sqrt(x_prime**2 + y_prime**2 + z_prime**2) # Should be ~1.0
    cos_theta_prime_arg = np.clip(z_prime / rho_prime, -1.0, 1.0)
    theta_prime = np.arccos(cos_theta_prime_arg)

    # Calculate new azimuthal angle phi'
    phi_prime = np.arctan2(y_prime, x_prime)

    # Adjust phi_prime range from [-pi, pi] to [0, 2*pi)
    phi_prime = np.mod(phi_prime, 2 * np.pi)

    # Handle poles explicitly: phi' is degenerate at the poles (theta'=0 or pi)
    # arctan2(0, 0) returns 0, which is a fine convention.
    # Check if any input points were exactly poles where phi is degenerate
    # if theta.shape:
    #     near_pole = (theta < epsilon) | (theta > np.pi - epsilon)
    #     # Or check if output points are near new poles
    #     near_new_pole = (theta_prime < epsilon) | (theta_prime > np.pi - epsilon)
    #     # Can set phi_prime to 0 for these cases if needed, but arctan2(0,0) handles it.
    # else:
    #     if (theta < epsilon) or (theta > np.pi - epsilon):
    #         pass # phi_prime already 0 from arctan2(0,0) if input is pole

    return theta_prime, phi_prime
