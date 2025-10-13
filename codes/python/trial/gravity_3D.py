#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Apr 12 10:37:30 2025

@author: g
"""

import numpy as np
import trimesh
import logging
import warnings
import csv
import matplotlib.pyplot as plt
import matplotlib.cm as cm

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

# --- Constants ---
G = 6.67430e-11  # Gravitational constant (N m^2 / kg^2)

# --- Helper Functions (calculate_Le, calculate_omega_f - unchanged from previous) ---
# ... (Keep the exact calculate_Le and calculate_omega_f functions from the previous answer) ...
def calculate_Le(r, p_i, p_j):
    """Calculates the L_e term (Eq. 11)."""
    x_i_vec = r - p_i
    x_j_vec = r - p_j
    x_i = np.linalg.norm(x_i_vec)
    x_j = np.linalg.norm(x_j_vec)
    e_ij_vec = p_j - p_i
    e_ij = np.linalg.norm(e_ij_vec)
    denominator = x_i + x_j - e_ij
    numerator = x_i + x_j + e_ij
    if x_i < 1e-10 or x_j < 1e-10: return 0.0
    if np.linalg.norm(np.cross(x_i_vec, e_ij_vec)) < 1e-10 * e_ij: return 0.0
    if abs(numerator) < 1e-10 or abs(denominator) < 1e-10: return 0.0
    if denominator <= 1e-10 or numerator <= 0: return 0.0
    return np.log(numerator / denominator)

def calculate_omega_f(r, p_i, p_j, p_k):
    """Calculates the solid angle ω_f term (Eq. 12)."""
    x_i_vec = r - p_i
    x_j_vec = r - p_j
    x_k_vec = r - p_k
    x_i = np.linalg.norm(x_i_vec)
    x_j = np.linalg.norm(x_j_vec)
    x_k = np.linalg.norm(x_k_vec)
    if x_i < 1e-10 or x_j < 1e-10 or x_k < 1e-10: return 0.0
    triple_product = np.dot(x_i_vec, np.cross(x_j_vec, x_k_vec))
    denominator = (x_i * x_j * x_k +
                   np.dot(x_i_vec, x_j_vec) * x_k +
                   np.dot(x_i_vec, x_k_vec) * x_j +
                   np.dot(x_j_vec, x_k_vec) * x_i)
    if abs(triple_product) < 1e-10 * (x_i*x_j*x_k): return 0.0
    if abs(denominator) < 1e-12: # Added check for very small denominator
         # logging.warning(f"omega_f denominator is very small ({denominator}) for r={r}, face ({p_i}, {p_j}, {p_k}). Potential instability. Returning 0.")
         return 0.0 # Avoid potential division issues with atan2 if denominator is tiny
    omega = 2.0 * np.arctan2(triple_product, denominator)
    if np.isnan(omega): return 0.0
    return omega
# --- Main Calculation Function (modified mesh loading slightly) ---

def calculate_polyhedral_gravity(mesh, density, evaluation_points):
    """
    Calculates gravitational acceleration using the polyhedral method.
    (Core logic mostly unchanged from previous version)

    Args:
        mesh (trimesh.Trimesh): Pre-loaded and validated trimesh object.
        density (float): Uniform density of the object (kg/m^3).
        evaluation_points (np.ndarray): Array of points (Nx3) where gravity
                                        should be calculated (meters).

    Returns:
        np.ndarray: Array of gravitational acceleration vectors (Nx3) in m/s^2.
    """
    logging.info("Extracting mesh topology...")
    edges = mesh.edges_unique
    edge_face_pairs = mesh.edges_face
    vertices = mesh.vertices
    faces = mesh.faces
    face_normals = mesh.face_normals # Outward-pointing normals

    logging.info("Starting gravity calculation loop...")
    gravity_vectors = np.zeros_like(evaluation_points, dtype=float)

    # Pre-calculate edge vectors for efficiency if needed (optional)
    # edge_vectors = vertices[edges[:, 1]] - vertices[edges[:, 0]]

    for idx, r in enumerate(evaluation_points):
        if idx > 0 and idx % 500 == 0: # Log progress less frequently for large meshes
             logging.info(f"Calculating gravity for point {idx}/{len(evaluation_points)}...")

        sum_edge_term = np.zeros(3, dtype=float)
        sum_face_term = np.zeros(3, dtype=float)

        # --- Loop over Edges ---
        for i_edge, edge_verts_idx in enumerate(edges):
            p_A = vertices[edge_verts_idx[0]]
            p_B = vertices[edge_verts_idx[1]]
            face_idx_1, face_idx_2 = edge_face_pairs[i_edge]

            if face_idx_1 == -1 or face_idx_2 == -1: continue # Skip boundary edges

            n_f1 = face_normals[face_idx_1]
            n_f2 = face_normals[face_idx_2]
            edge_vec = p_B - p_A # Vector along the edge

            # Calculate n_e for each face (perpendicular to edge in face plane)
            # Note: Orientation relative to edge direction might matter based on strict interpretation
            # Standard methods often use dyadic products of face normal & edge vector directly.
            # Let's stick to Eq.8's formula:
            n_e_f1 = np.cross(edge_vec, n_f1)
            n_e_f2 = np.cross(edge_vec, n_f2)
            norm_n_e_f1 = np.linalg.norm(n_e_f1)
            norm_n_e_f2 = np.linalg.norm(n_e_f2)
            if norm_n_e_f1 > 1e-9: n_e_f1 /= norm_n_e_f1
            else: continue
            if norm_n_e_f2 > 1e-9: n_e_f2 /= norm_n_e_f2
            else: continue

            E_e = np.outer(n_f1, n_e_f1) + np.outer(n_f2, n_e_f2)
            L_e = calculate_Le(r, p_A, p_B)
            x_A_vec = r - p_A
            if np.isfinite(L_e):
                 sum_edge_term += (E_e @ x_A_vec) * L_e
            # else: logging.debug(...) # Reduce verbosity of warnings

        # --- Loop over Faces ---
        for i_face, face_verts_idx in enumerate(faces):
            p_i = vertices[face_verts_idx[0]]
            p_j = vertices[face_verts_idx[1]]
            p_k = vertices[face_verts_idx[2]]
            n_f = face_normals[i_face]
            F_f = np.outer(n_f, n_f)
            omega_f = calculate_omega_f(r, p_i, p_j, p_k)
            x_i_vec = r - p_i
            if np.isfinite(omega_f):
                sum_face_term += (F_f @ x_i_vec) * omega_f
            # else: logging.debug(...) # Reduce verbosity

        # --- Final Gravity Vector ---
        if not np.all(np.isfinite(sum_edge_term)) or not np.all(np.isfinite(sum_face_term)):
             logging.warning(f"Non-finite sum term encountered for point {idx} ({r}). Setting gravity to zero.")
             gravity_vectors[idx] = np.zeros(3)
        else:
             gravity_vectors[idx] = -G * density * (sum_edge_term - sum_face_term)

    logging.info(f"Gravity calculation finished for {len(evaluation_points)} points.")
    return gravity_vectors

# --- Saving Functions ---

def save_gravity_data_csv(filename, points, vectors):
    """Saves evaluation points and gravity vectors to a CSV file."""
    if points.shape[0] != vectors.shape[0]:
        logging.error("Mismatch between number of points and vectors. Cannot save CSV.")
        return
    try:
        header = ['point_x', 'point_y', 'point_z', 'gravity_x', 'gravity_y', 'gravity_z']
        data_to_save = np.hstack((points, vectors))
        with open(filename, 'w', newline='') as f:
            writer = csv.writer(f)
            writer.writerow(header)
            writer.writerows(data_to_save)
        logging.info(f"Gravity data saved to CSV: {filename}")
    except Exception as e:
        logging.error(f"Failed to save data to CSV {filename}: {e}")

def save_gravity_data_npz(filename, points, vectors):
    """Saves evaluation points and gravity vectors to a NumPy NPZ file."""
    if points.shape[0] != vectors.shape[0]:
        logging.error("Mismatch between number of points and vectors. Cannot save NPZ.")
        return
    try:
        np.savez_compressed(filename, evaluation_points=points, gravity_vectors=vectors)
        logging.info(f"Gravity data saved to NPZ: {filename}")
    except Exception as e:
        logging.error(f"Failed to save data to NPZ {filename}: {e}")


# --- Plotting Function ---

def plot_normal_gravity(mesh, points, gravity_vectors, surface_normals):
    """
    Plots the mesh colored by the component of gravity normal to the surface.

    Args:
        mesh (trimesh.Trimesh): The mesh object.
        points (np.ndarray): Points where gravity was evaluated (should correspond to normals).
        gravity_vectors (np.ndarray): Calculated gravity vectors at the points.
        surface_normals (np.ndarray): Surface normal vectors at the points.
    """
    if not all([points.shape[0] == gravity_vectors.shape[0],
                points.shape[0] == surface_normals.shape[0]]):
        logging.error("Mismatch in array lengths for plotting. Aborting plot.")
        return

    logging.info("Calculating normal component of gravity...")
    # Project gravity vector onto the normal vector: g_normal = dot(g, n)
    # A negative value means gravity points opposite to the normal (i.e., inwards)
    g_normal_component = np.sum(gravity_vectors * surface_normals, axis=1)

    # We typically want to visualize the *inward pull*, so plot -g_normal_component
    # Higher positive values will mean stronger inward gravity.
    values_to_plot = -g_normal_component
    min_val = np.min(values_to_plot)
    max_val = np.max(values_to_plot)
    logging.info(f"Normal gravity component range (inward positive): Min={min_val:.3e}, Max={max_val:.3e} m/s^2")

    # Normalize values to 0-1 range for colormapping
    norm = plt.Normalize(vmin=min_val, vmax=max_val)
    colormap = cm.viridis # Or choose another colormap like 'plasma', 'inferno', 'magma'

    # Map normalized values to RGBA colors
    colors = colormap(norm(values_to_plot))

    # --- Assign colors to the mesh ---
    # We assume points correspond directly to faces (if using centroids)
    # or vertices (if using vertices). Check the points generation logic.
    # If points are face centroids:
    mesh.visual.face_colors = colors
    # If points are vertices:
    # mesh.visual.vertex_colors = colors

    logging.info("Displaying mesh with normal gravity colormap...")

    # Create a scene and add the colored mesh
    scene = trimesh.Scene(mesh)

    # --- Add a color bar ---
    # Create a dummy scalar mappable for the colorbar
    sm = plt.cm.ScalarMappable(cmap=colormap, norm=norm)
    sm.set_array([]) # You need to set an array for the mappable

    # Getting the colorbar figure (can be shown separately or saved)
    fig, ax = plt.subplots(figsize=(1, 6)) # Small figure for the colorbar
    cbar = plt.colorbar(sm, cax=ax)
    cbar.set_label('Inward Normal Gravity Component (m/s^2)')
    fig.suptitle('Color Scale', fontsize=10)
    fig.tight_layout(rect=[0, 0.03, 1, 0.95]) # Adjust layout

    # Show the main scene (this blocks until closed)
    scene.show()

    # Show the color bar figure (doesn't block)
    plt.show() # Show the separate colorbar window


# --- Main Execution Block ---
if __name__ == "__main__":
    obj_file = '/home/g/Asteroids/input/Bennu.obj' # <--- REPLACE WITH YOUR OBJ FILE PATH
    uniform_density = 1250.0   # kg/m^3 (Example: typical rock density)
    output_csv_file = 'gravity_results.csv'
    output_npz_file = 'gravity_results.npz'

    # --- Load and Prepare Mesh ---
    logging.info(f"Loading mesh from {obj_file}...")
    mesh = None
    try:
        mesh = trimesh.load(obj_file, force='mesh', process=True)
        if not isinstance(mesh, trimesh.Trimesh):
             logging.error(f"Loaded object is not a Trimesh instance (type: {type(mesh)}). Cannot proceed.")
             mesh = None
        else:
             logging.info("Mesh loaded.")
             logging.info(f"  Vertices: {len(mesh.vertices)}")
             logging.info(f"  Faces: {len(mesh.faces)}")
             # Validate and potentially fix
             if not mesh.is_watertight:
                 logging.warning("Mesh not watertight. Attempting to fill holes...")
                 mesh.fill_holes()
                 if not mesh.is_watertight:
                     logging.warning("Mesh still not watertight after attempting fix.")
             if not mesh.is_orientable:
                 logging.warning("Mesh not orientable. Attempting to fix windings.")
                 mesh.fix_normals()
             mesh.fix_normals() # Ensure normals are consistent after fixes
             logging.info("Mesh validation/fixing complete.")

    except Exception as e:
        logging.error(f"Failed to load or process mesh: {e}")
        mesh = None

    if mesh is None:
        exit("Exiting due to mesh loading failure.")

    # --- Define Evaluation Points ---
    # **RECOMMENDED**: Use face centroids slightly offset along the normal
    # This avoids singularities at vertices/edges and provides normals easily.
    offset_distance = 1e-7 # meters (very small offset)
    evaluation_points = mesh.triangles_center + mesh.face_normals * offset_distance
    corresponding_normals = mesh.face_normals # Normals corresponding to faces
    logging.info(f"Calculating gravity at {len(evaluation_points)} face centroids (offset by {offset_distance}m).")

    # Alternative: Use vertices (might be numerically unstable)
    # evaluation_points = mesh.vertices
    # corresponding_normals = mesh.vertex_normals # Normals corresponding to vertices
    # logging.info(f"Calculating gravity at {len(evaluation_points)} mesh vertices.")


    # --- Run Calculation ---
    gravity_results = None
    with warnings.catch_warnings():
         warnings.simplefilter("ignore", category=RuntimeWarning) # Suppress runtime warnings during calc
         gravity_results = calculate_polyhedral_gravity(mesh, uniform_density, evaluation_points)


    # --- Save and Plot Results ---
    if gravity_results is not None:
        # Save results
        save_gravity_data_csv(output_csv_file, evaluation_points, gravity_results)
        save_gravity_data_npz(output_npz_file, evaluation_points, gravity_results)

        # Plot results
        plot_normal_gravity(mesh, evaluation_points, gravity_results, corresponding_normals)

        logging.info("Processing complete.")
    else:
        logging.error("Gravity calculation failed. No results to save or plot.")