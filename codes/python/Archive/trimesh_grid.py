#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Mar 30 05:59:04 2025

@author: g
"""

import numpy as np
import trimesh

# --- Helper Function for Coordinate Conversion ---
def spherical_to_cartesian(rho, theta, phi):
    """
    Converts spherical coordinates (ISO convention: radius, inclination, azimuth)
    to Cartesian coordinates (x, y, z).
    Assumes theta is inclination (polar angle, from +Z) in [0, pi]
    Assumes phi is azimuth (from +X towards +Y) in [0, 2*pi]
    """
    x = rho * np.sin(theta) * np.cos(phi)
    y = rho * np.sin(theta) * np.sin(phi)
    z = rho * np.cos(theta)
    return np.array([x, y, z])

# --- Crater Function ---
def add_gaussian_crater(mesh, center_theta, center_phi, crater_radius_angle, max_depth, rim_height=0.0, rim_width_factor=1.5):
    """
    Adds a Gaussian-shaped crater depression and an optional raised rim to a mesh.

    Args:
        mesh (trimesh.Trimesh): The mesh object to modify (should be sphere-like).
        center_theta (float): Inclination angle (polar angle, from +Z) of the crater center [0, pi].
        center_phi (float): Azimuthal angle (from +X) of the crater center [0, 2*pi].
        crater_radius_angle (float): The angular radius of the crater depression (e.g., in radians).
                                     Controls the width of the Gaussian sigma.
        max_depth (float): Maximum depth of the crater depression (positive value).
        rim_height (float): Maximum height of the raised rim (positive value). Defaults to 0.0 (no rim).
        rim_width_factor (float): How much wider the rim effect extends compared to the
                                  depression radius (e.g., 1.5 means rim extends to
                                  1.5 * crater_radius_angle). Defaults to 1.5.

    Returns:
        trimesh.Trimesh: The modified mesh object.
    """
    if not isinstance(mesh, trimesh.Trimesh):
        raise TypeError("Input mesh must be a trimesh.Trimesh object")
    if mesh.is_empty:
        raise ValueError("Input mesh is empty")

    print(f"Adding crater at theta={np.degrees(center_theta):.1f} deg, phi={np.degrees(center_phi):.1f} deg")

    # Assuming the mesh represents a sphere centered at the origin
    # Calculate the approximate radius of the sphere from vertex distances
    sphere_radius = np.mean(np.linalg.norm(mesh.vertices, axis=1))
    if sphere_radius < 1e-6:
        raise ValueError("Could not determine sphere radius or radius is near zero.")

    # Convert crater center spherical coordinates to Cartesian
    impact_center_cartesian = spherical_to_cartesian(sphere_radius, center_theta, center_phi)
    impact_center_normalized = impact_center_cartesian / np.linalg.norm(impact_center_cartesian)

    # Define Gaussian sigma based on angular radius (controls crater width)
    # A smaller sigma makes a sharper crater for the same radius_angle
    sigma = crater_radius_angle / 2.0  # Adjust this factor as needed

    # Define angular radius for the rim effect
    rim_radius_angle = crater_radius_angle * rim_width_factor

    modified_vertices = np.copy(mesh.vertices) # Work on a copy

    # Iterate through each vertex
    for i, vertex in enumerate(mesh.vertices):
        vertex_normalized = vertex / np.linalg.norm(vertex)

        # Calculate the angle between the vertex and the impact center
        # Using dot product: dot(A, B) = |A||B|cos(angle)
        dot_product = np.dot(vertex_normalized, impact_center_normalized)
        dot_product = np.clip(dot_product, -1.0, 1.0) # Clamp for numerical stability
        angle_to_center = np.arccos(dot_product)

        # Calculate displacement only if the vertex is reasonably close to the crater center
        # Check against the larger rim radius if a rim exists, otherwise the crater radius
        check_radius = rim_radius_angle if rim_height > 0 else crater_radius_angle * 1.5 # Extend check slightly beyond crater
        if angle_to_center < check_radius: # Only modify vertices within influence zone

            # --- Calculate Depression ---
            # Gaussian function for depth based on angle
            # depth = max_depth * exp(-angle^2 / (2 * sigma^2))
            depression_depth = max_depth * np.exp(-angle_to_center**2 / (2 * sigma**2))

            # --- Calculate Rim (Optional) ---
            rim_displacement = 0.0
            if rim_height > 0 and angle_to_center > sigma : # Only apply rim outside the central depression zone sigma
                 # Use another Gaussian or similar function for the rim shape, peaking near crater_radius_angle
                 # Example: A Gaussian centered around crater_radius_angle
                 rim_sigma = (rim_radius_angle - crater_radius_angle) / 2.5 # Width of the rim bump
                 if rim_sigma > 1e-6: # Avoid division by zero if rim_width_factor is too small
                     rim_displacement = rim_height * np.exp(-(angle_to_center - crater_radius_angle)**2 / (2 * rim_sigma**2))
                     # Optional: Taper rim effect smoothly to zero at rim_radius_angle
                     taper = np.clip(1.0 - (angle_to_center - crater_radius_angle) / (rim_radius_angle - crater_radius_angle + 1e-6) , 0.0, 1.0)
                     rim_displacement *= taper


            # Total displacement: Rim pushes outwards, Depression pushes inwards
            # Displacement is along the surface normal (which is the vertex direction for a sphere at origin)
            total_displacement = rim_displacement - depression_depth
            displacement_vector = vertex_normalized * total_displacement

            # Apply the displacement
            modified_vertices[i] += displacement_vector

    # Create a new mesh with the modified vertices
    # Using process=False speeds things up as we assume topology doesn't change drastically
    new_mesh = trimesh.Trimesh(vertices=modified_vertices, faces=mesh.faces, process=False)
    print("Crater addition complete.")
    return new_mesh

# --- Main Example ---
if __name__ == "__main__":
    # 1. Create a base sphere mesh
    # Icosphere provides more uniform triangulation than UV sphere
    sphere = trimesh.primitives.Sphere(radius=1.0, subdivisions=5) # Higher subdivisions = more detail
    print(f"Created sphere with {len(sphere.vertices)} vertices.")

    # --- Define Craters ---
    # Angles should be in RADIANS for numpy functions
    # Format: (theta_rad, phi_rad, angular_radius_rad, depth, rim_height, rim_width_factor)
    craters_to_add = [
        # Large crater near North Pole
        (np.radians(20), np.radians(45), np.radians(15), 0.15, 0.01, 1.5),
        # Smaller crater in mid-latitudes
        (np.radians(70), np.radians(180), np.radians(8), 0.08, 0.005, 1.4),
        # Tiny crater near South Pole with no rim
        (np.radians(160), np.radians(270), np.radians(3), 0.03, 0.0, 1.5),
         # Crater right on the equator
        (np.radians(90), np.radians(0), np.radians(10), 0.1, 0.008, 1.6),
    ]

    # 3. Add craters iteratively
    cratere_mesh = sphere # Start with the base sphere
    for i, params in enumerate(craters_to_add):
        print(f"\n--- Adding Crater {i+1} ---")
        theta, phi, radius_angle, depth, rim_h, rim_w = params
        cratere_mesh = add_gaussian_crater(
            cratere_mesh,
            center_theta=theta,
            center_phi=phi,
            crater_radius_angle=radius_angle,
            max_depth=depth,
            rim_height=rim_h,
            rim_width_factor=rim_w
        )
        # Optional: Recalculate normals if needed for lighting/shading downstream
        # cratere_mesh.fix_normals()


    # 4. Visualize the result
    print("\nDisplaying final mesh...")
    # Ensure vertices have distinct colors from faces for better visualization
    # cratere_mesh.visual.face_colors = [100, 100, 100, 255] # Grey faces
    cratere_mesh.show()

    # Optional: Save the mesh
    # try:
    #     cratere_mesh.export("sphere_with_craters.stl")
    #     print("Mesh saved to sphere_with_craters.stl")
    # except Exception as e:
    #     print(f"Could not save mesh: {e}")