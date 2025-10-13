#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 10 07:21:50 2025

@author: g
"""

import numpy as np
import trimesh
from scipy.optimize import minimize
import matplotlib.pyplot as plt
import math
import os
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from scipy.interpolate import griddata
import warnings

def interpolate_nan_deviations(lats_rad, lons_rad, deviations, method='linear', padding_lon=5):
    """
    Interpolates NaN values in a 2D deviation grid defined on a sphere.

    Handles longitude wrap-around by padding the data.

    Args:
        lats_rad (np.ndarray): 2D array of latitudes (radians).
        lons_rad (np.ndarray): 2D array of longitudes (radians).
        deviations (np.ndarray): 2D array of deviations, may contain NaNs.
        method (str): Interpolation method for griddata ('linear', 'nearest', 'cubic').
                      'linear' is generally recommended.
        padding_lon (int): Number of longitude columns to pad on each side
                           to handle wrap-around during interpolation.

    Returns:
        np.ndarray: A copy of the deviations array with NaNs filled by interpolation,
                    or the original array if no NaNs were present.
    """
    if not np.isnan(deviations).any():
        print("No NaN values found in deviations grid. Returning original.")
        return deviations.copy() # Return a copy

    n_lat, n_lon = deviations.shape
    if padding_lon >= n_lon // 2:
        print(f"Warning: Longitude padding ({padding_lon}) is large relative to grid width ({n_lon}). Adjusting.")
        padding_lon = max(1, n_lon // 4) # Use a smaller padding if too large

    print(f"Interpolating {np.isnan(deviations).sum()} NaN values using '{method}' method...")

    # --- 1. Create Padded Arrays to Handle Longitude Wrap-around ---

    # Pad deviations array
    padded_deviations = np.pad(deviations, ((0, 0), (padding_lon, padding_lon)), mode='wrap')

    # Create corresponding padded longitude array
    lon_step = lons_rad[0, 1] - lons_rad[0, 0] # Assumes regular lon spacing
    left_pad_lons = lons_rad[:, -padding_lon:] - 2 * np.pi
    right_pad_lons = lons_rad[:, :padding_lon] + 2 * np.pi
    padded_lons_rad = np.concatenate((left_pad_lons, lons_rad, right_pad_lons), axis=1)

    # Latitude array doesn't need padding in the same way, but needs matching shape
    padded_lats_rad = np.pad(lats_rad, ((0, 0), (padding_lon, padding_lon)), mode='edge') # Pad lats using edge values

    # --- 2. Identify Valid and Invalid Points in the PADDED Grid ---
    valid_mask_padded = ~np.isnan(padded_deviations)
    invalid_mask_original = np.isnan(deviations) # Identify NaNs in the ORIGINAL grid

    # Coordinates of known points (from padded grid)
    points_known = np.stack((padded_lats_rad[valid_mask_padded],
                             padded_lons_rad[valid_mask_padded]), axis=-1)

    # Values at known points
    values_known = padded_deviations[valid_mask_padded]

    # Coordinates where we need to interpolate (from original grid)
    points_to_interpolate = np.stack((lats_rad[invalid_mask_original],
                                      lons_rad[invalid_mask_original]), axis=-1)

    if points_to_interpolate.shape[0] == 0:
         print("Internal check: No points need interpolation after masking? Returning original.")
         return deviations.copy()
    if points_known.shape[0] < 3 and method != 'nearest':
         print(f"Warning: Not enough valid data points ({points_known.shape[0]}) for '{method}' interpolation. Trying 'nearest'.")
         method = 'nearest'
    if points_known.shape[0] == 0:
        print("Error: No valid data points found to interpolate from.")
        return deviations.copy() # Cannot interpolate

    # --- 3. Perform Interpolation ---
    print(f"Interpolating {points_to_interpolate.shape[0]} points from {points_known.shape[0]} valid points...")
    start_time = time.time() # Requires `import time`
    interpolated_values = griddata(points_known, values_known, points_to_interpolate, method=method)
    print(f"Interpolation finished in {time.time() - start_time:.2f} seconds.")


    # --- 4. Fill NaNs in a Copy of the Original Array ---
    filled_deviations = deviations.copy()
    nan_indices = np.where(invalid_mask_original) # Get indices of NaNs

    # Check if interpolation produced NaNs (can happen with 'linear'/'cubic' if outside convex hull)
    num_interp_nan = np.isnan(interpolated_values).sum()
    if num_interp_nan > 0:
        print(f"Warning: Interpolation resulted in {num_interp_nan} NaN values.")
        # Option 1: Leave them as NaN
        # Option 2: Try 'nearest' for the remaining NaNs (more robust fallback)
        if method != 'nearest':
            print("Attempting 'nearest' neighbor interpolation for remaining NaNs...")
            points_nan_interp = points_to_interpolate[np.isnan(interpolated_values)]
            nearest_values = griddata(points_known, values_known, points_nan_interp, method='nearest')
            # Fill only where the original interpolation failed
            original_nan_indices_to_fill = (nan_indices[0][np.isnan(interpolated_values)],
                                            nan_indices[1][np.isnan(interpolated_values)])
            filled_deviations[original_nan_indices_to_fill] = nearest_values
            # Fill the ones that succeeded with the primary method
            valid_interp_mask = ~np.isnan(interpolated_values)
            original_nan_indices_succeeded = (nan_indices[0][valid_interp_mask],
                                              nan_indices[1][valid_interp_mask])
            filled_deviations[original_nan_indices_succeeded] = interpolated_values[valid_interp_mask]
        else:
             # If even nearest failed, something is fundamentally wrong, leave NaN
             filled_deviations[nan_indices] = interpolated_values # Fill anyway, might contain NaNs
    else:
        # All points interpolated successfully
        filled_deviations[nan_indices] = interpolated_values

    # Final check
    final_nan_count = np.isnan(filled_deviations).sum()
    if final_nan_count > 0:
         print(f"Warning: {final_nan_count} NaN values remain after interpolation attempts.")
    else:
         print("All NaN values successfully interpolated.")


    return filled_deviations


def intersect_chunk(mesh_data, ray_origins_chunk, ray_directions_chunk, original_indices):
    """Worker function to intersect a chunk of rays."""
    # Important: Pass necessary mesh data or reload mesh if needed
    # Re-creating the mesh object might be necessary depending on how processes work
    # For simplicity, assume mesh object can be pickled or reconstruct it
    # mesh = trimesh.Trimesh(**mesh_data) # Example reconstruction
    mesh = trimesh.load('/home/g/Asteroids/input/small_Bennu.obj', force='mesh') # Or reload

    locations, index_ray_chunk, index_tri = mesh.ray.intersects_location(
        ray_origins=ray_origins_chunk,
        ray_directions=ray_directions_chunk,
        multiple_hits=False
    )
    # Map local chunk indices back to original global indices
    index_ray_global = original_indices[index_ray_chunk]
    return locations, index_ray_global

def calculate_surface_deviations_parallel(mesh, center, radius, grid_points_xyz, lats_rad, lons_rad, num_workers=24):
    n_points = grid_points_xyz.shape[0]
    deviations_flat = np.full(n_points, np.nan)
    ray_origins = np.tile(center, (n_points, 1))
    ray_directions = grid_points_xyz - center
    norms = np.linalg.norm(ray_directions, axis=1)
    valid_norms = norms > 1e-9
    ray_directions[valid_norms] /= norms[valid_norms, np.newaxis]

    # Determine chunk size (e.g., 10000 rays per chunk)
    chunk_size = 1000
    num_chunks = int(np.ceil(n_points / chunk_size))

    all_locations = []
    all_index_ray = []

    start_time = time.time()
    # Use ProcessPoolExecutor for CPU-bound tasks
    with ProcessPoolExecutor(max_workers=num_workers) as executor:
        futures = []
        for i in range(num_chunks):
            start_idx = i * chunk_size
            end_idx = min((i + 1) * chunk_size, n_points)
            
            # Get chunks
            origins_chunk = ray_origins[start_idx:end_idx]
            directions_chunk = ray_directions[start_idx:end_idx]
            original_indices_chunk = np.arange(start_idx, end_idx) # Track original indices

            # Prepare mesh data (only if mesh object isn't easily picklable)
            # mesh_data = {'vertices': mesh.vertices, 'faces': mesh.faces} # Example

            # Submit task
            # futures.append(executor.submit(intersect_chunk, mesh_data, origins_chunk, directions_chunk, original_indices_chunk))
            # If mesh object is picklable (often is), you might pass it directly BUT test carefully
            # Be cautious passing large objects - reloading might be safer/faster if pickling is slow
            futures.append(executor.submit(intersect_chunk, None, origins_chunk, directions_chunk, original_indices_chunk)) # Assuming reload in worker

        # Collect results as they complete
        for future in as_completed(futures):
            try:
                locations_chunk, index_ray_global_chunk = future.result()
                if len(locations_chunk) > 0:
                    all_locations.append(locations_chunk)
                    all_index_ray.append(index_ray_global_chunk)
            except Exception as exc:
                print(f'A chunk generated an exception: {exc}')

    print(f"Parallel ray casting took: {time.time() - start_time:.2f}s")

    if not all_locations:
         print("Warning: No rays hit the mesh surface in parallel execution.")
         return deviations_flat.reshape(lats_rad.shape) # Return NaN array

    # Combine results from all chunks
    locations = np.vstack(all_locations)
    index_ray = np.concatenate(all_index_ray)
    print(f"Parallel processing: {len(locations)} rays hit the surface out of {n_points}.")

    if len(locations) > 0:
        surface_distances = np.linalg.norm(locations - center, axis=1)
        hit_deviations = surface_distances - radius
        deviations_flat[index_ray] = hit_deviations

    return deviations_flat.reshape(lats_rad.shape)


# --- Helper Functions ---

def load_asteroid_mesh(filepath):
    """Loads an asteroid mesh from OBJ or GLB/glTF file."""
    try:
        # force='mesh' ensures we get a single mesh geometry
        mesh = trimesh.load(filepath, force='mesh')
        print(f"Successfully loaded mesh from: {filepath}")
        print(f"Mesh has {len(mesh.vertices)} vertices and {len(mesh.faces)} faces.")
        # If the mesh has multiple disconnected parts, combine them
        if isinstance(mesh, trimesh.Scene):
             # Combine all meshes in the scene into a single mesh
             mesh = trimesh.util.concatenate(mesh.dump())
        elif isinstance(mesh, list):
             # If load returns a list of meshes
             mesh = trimesh.util.concatenate(mesh)

        # Ensure vertices are float64 for stability in calculations
        mesh.vertices = mesh.vertices.astype(np.float64)

        # Sometimes models are not centered, pre-center for potentially better fitting
        # mesh.apply_translation(-mesh.centroid) # Optional: Center on centroid first
        return mesh
    except Exception as e:
        print(f"Error loading mesh from {filepath}: {e}")
        return None

def fit_sphere_least_squares(points):
    """
    Fits a sphere to a set of 3D points using least squares minimization.
    Returns the center (x, y, z) and radius of the best-fit sphere.
    """
    if points.shape[0] < 4:
        raise ValueError("Need at least 4 points to fit a sphere")

    # Initial guess: Centroid and average distance to centroid
    center_init = np.mean(points, axis=0)
    distances_init = np.linalg.norm(points - center_init, axis=1)
    radius_init = np.mean(distances_init)
    
    # --- Objective function ---
    # Minimize the sum of squared differences between point distances
    # to the center and the sphere radius.
    def sphere_error(params, points):
        cx, cy, cz, r = params
        center = np.array([cx, cy, cz])
        # Ensure radius is positive during optimization
        r_positive = np.abs(r) 
        distances = np.linalg.norm(points - center, axis=1)
        return np.sum((distances - r_positive)**2)

    # --- Optimization ---
    initial_guess = np.append(center_init, radius_init)
    result = minimize(
        sphere_error,
        initial_guess,
        args=(points,),
        method='L-BFGS-B', # Algorithm supporting bounds if needed
        options={'ftol': 1e-6, 'maxiter': 1000} 
    )

    if not result.success:
        print(f"Warning: Sphere fitting optimization did not converge successfully. Status: {result.status}, Message: {result.message}")
        # Fallback or raise error? Using result anyway for now.
        # raise RuntimeError("Sphere fitting optimization failed.")


    center_fit = result.x[:3]
    radius_fit = np.abs(result.x[3]) # Ensure radius is positive

    print(f"Best-fit sphere center: {center_fit}")
    print(f"Best-fit sphere radius: {radius_fit}")
    return center_fit, radius_fit

def create_spherical_grid(center, radius, n_lat, n_lon):
    """
    Generates points on a sphere surface based on latitude/longitude.

    Args:
        center (np.array): Sphere center [x, y, z].
        radius (float): Sphere radius.
        n_lat (int): Number of latitude lines (excluding poles potentially).
        n_lon (int): Number of longitude lines.

    Returns:
        tuple: (lats_rad, lons_rad, grid_points_xyz)
               lats_rad: 2D array of latitudes in radians.
               lons_rad: 2D array of longitudes in radians.
               grid_points_xyz: Nx3 array of Cartesian coordinates on the sphere.
    """
    # Create latitude and longitude arrays (radians)
    # Latitude from -pi/2 (South Pole) to +pi/2 (North Pole)
    # Linspace includes endpoints, adjust number if poles shouldn't be duplicated across longitude
    lat = np.linspace(0, np.pi, n_lat)
    # Longitude from -pi (-180 deg) to +pi (+180 deg)
    # Endpoint=False to avoid duplicating the -180/180 meridian
    lon = np.linspace(0, 2*np.pi, n_lon, endpoint=False)

    # Create a meshgrid
    # Note: ordering='ij' makes lat the first index, lon the second
    lats_rad, lons_rad = np.meshgrid(lat, lon,indexing='ij')

    # Convert spherical coordinates (lat, lon, radius) to Cartesian (x, y, z)
    x = center[0] + radius * np.sin(lats_rad) * np.cos(lons_rad)
    y = center[1] + radius * np.sin(lats_rad) * np.sin(lons_rad)
    z = center[2] + radius * np.cos(lats_rad)

    # Stack coordinates into an (N, 3) array where N = n_lat * n_lon
    grid_points_xyz = np.vstack([x.ravel(), y.ravel(), z.ravel()]).T

    return lats_rad, lons_rad, grid_points_xyz


def calculate_surface_deviations(mesh, center, radius, grid_points_xyz, lats_rad, lons_rad):
    """
    Calculates the distance from the center to the mesh surface along rays
    defined by the spherical grid points, and finds the deviation from the
    best-fit sphere radius.

    Args:
        mesh (trimesh.Trimesh): The asteroid mesh.
        center (np.array): The sphere center.
        radius (float): The sphere radius.
        grid_points_xyz (np.array): Nx3 array of grid points on the sphere surface.
        lats_rad (np.array): 2D array of latitudes (radians).
        lons_rad (np.array): 2D array of longitudes (radians).


    Returns:
        np.array: A 2D array (matching lats/lons shape) containing the deviation
                  (surface_distance - sphere_radius) for each grid point.
                  NaN indicates the ray did not intersect the mesh.
    """
    n_points = grid_points_xyz.shape[0]
    deviations_flat = np.full(n_points, np.nan) # Initialize with NaN

    # --- Ray casting ---
    # Ray origins are all the sphere center
    ray_origins = np.tile(center, (n_points, 1))
    # Ray directions point from the center towards each grid point
    ray_directions = grid_points_xyz - center
    # Normalize directions (important!)
    norms = np.linalg.norm(ray_directions, axis=1)
    # Avoid division by zero for potential center points (shouldn't happen with grid)
    valid_norms = norms > 1e-9 
    ray_directions[valid_norms] /= norms[valid_norms, np.newaxis]
    
    # Use trimesh's ray intersection
    # `locations`: intersection points in 3D space
    # `index_ray`: index of the ray that hit something (maps back to our grid_points_xyz)
    # `index_tri`: index of the triangle face that was hit
    locations, index_ray, index_tri = mesh.ray.intersects_location(
        ray_origins=ray_origins,
        ray_directions=ray_directions,
        multiple_hits=False # We only want the first hit from the center
    )

    print(f"Ray casting performed. {len(locations)} rays hit the surface out of {n_points}.")

    if len(locations) > 0:
        # Calculate the distance from the center to each intersection point
        surface_distances = np.linalg.norm(locations - center, axis=1)

        # Calculate the deviation: surface_distance - sphere_radius
        hit_deviations = surface_distances - radius

        # Place the calculated deviations into the correct spot in the flat array
        # using the index_ray mapping
        deviations_flat[index_ray] = hit_deviations
    else:
        print("Warning: No rays hit the mesh surface. Check mesh orientation, center, or ray directions.")

    # Reshape the flat deviations array back into the 2D grid shape
    deviations_2d = deviations_flat.reshape(lats_rad.shape)

    return deviations_2d
#%%

def plot_deviation_map(lats_rad_ignored, lons_rad_ignored, deviations, title="Asteroid Surface Deviation"):
    """Plots the deviation data as a 2D map using pcolormesh.
    Note: lats_rad_ignored and lons_rad_ignored are no longer directly used
          for plotting but kept for API consistency if needed elsewhere.
    """

    if deviations is None or np.all(np.isnan(deviations)):
        print("Warning: Deviation data is empty or all NaN. Skipping plot.")
        return

    # Get dimensions FROM THE DATA ARRAY (deviations)
    # This ensures consistency regardless of initial N_LATITUDE/N_LONGITUDE settings,
    # assuming deviations correctly reflects the (latitude, longitude) structure.
    n_lat, n_lon = deviations.shape
    print(f"Plotting deviations with shape (n_lat={n_lat}, n_lon={n_lon})") # Debug print

    # --- Calculate corner coordinates required by pcolormesh ---

    # Calculate the step size for lat and lon based on the number of cells
    # Handle edge case of only 1 division
    dlat = (np.pi / (n_lat - 1)) if n_lat > 1 else np.pi
    dlon = (2 * np.pi / n_lon) if n_lon > 0 else (2 * np.pi)

    # Create 1D boundary arrays (size N+1)
    lat_bounds = np.linspace(0, np.pi , n_lat + 1) # Correct number of boundaries for n_lat cells
    lon_bounds = np.linspace(0, 2*np.pi, n_lon + 1) # Correct number of boundaries for n_lon cells

    # Create 2D meshgrid for the corners. `indexing='ij'` makes the first dimension
    # correspond to the first input array (`lat_bounds`), matching `deviations`.
    lat_corners_rad, lon_corners_rad = np.meshgrid(lat_bounds, lon_bounds, indexing='ij')
    # Expected shapes: lat_corners_rad (n_lat+1, n_lon+1), lon_corners_rad (n_lat+1, n_lon+1)

    # Convert corner coordinates to degrees for plotting
    lons_deg_corners = np.rad2deg(lon_corners_rad)
    lats_deg_corners = np.rad2deg(lat_corners_rad)

    # --- Create the plot ---
    fig, ax = plt.subplots(figsize=(10, 5))

    # pcolormesh expects X, Y, C
    # X: Longitude corners (shape n_lat+1, n_lon+1)
    # Y: Latitude corners (shape n_lat+1, n_lon+1)
    # C: Deviation data (shape n_lat, n_lon)
    # This now matches the requirement that C is one smaller in each dim than X, Y.
    im = ax.pcolormesh(lons_deg_corners, lats_deg_corners, deviations,
                       shading='flat', # 'flat' needs dimensions (M+1, N+1) for X,Y and (M,N) for C
                       cmap='coolwarm',
                       vmin=np.nanmin(deviations), vmax=np.nanmax(deviations))

    ax.set_xlabel("Longitude (degrees)")
    ax.set_ylabel("Latitude (degrees)")
    ax.set_title(title)

    ax.set_xlim(0, 360)
    ax.set_ylim(0, 180)
    ax.set_xticks(np.linspace(0, 360, 7))
    ax.set_yticks(np.linspace(0, 180, 7))

    cbar = fig.colorbar(im, ax=ax)
    cbar.set_label("Deviation (Surface Distance - Sphere Radius) [units of mesh]")

    plt.tight_layout()
    plt.show()


def save_asteroid_data(filename, reference_radius, lats_rad, lons_rad, deviations):
    """
    Saves the reference radius, latitude, longitude, and deviation grids
    to a compressed NumPy (.npz) file.

    Args:
        filename (str): The path to the file to save (e.g., 'asteroid_data.npz').
        reference_radius (float): The reference radius used.
        lats_rad (np.ndarray): 2D array of latitudes in radians.
        lons_rad (np.ndarray): 2D array of longitudes in radians.
        deviations (np.ndarray): 2D array of height deviations.

    Returns:
        bool: True if saving was successful, False otherwise.
    """
    try:
        # Use savez_compressed for smaller file size
        np.savez_compressed(
            filename,
            radius=np.array(reference_radius), # Store radius as a 0-D numpy array
            lats_rad=lats_rad,
            lons_rad=lons_rad,
            deviations=deviations
        )
        print(f"Successfully saved asteroid data to: {filename}")
        return True
    except Exception as e:
        print(f"Error saving asteroid data to {filename}: {e}")
        return False
    
#%%

# --- Main Execution ---
if __name__ == "__main__":

    asteroid_filepath = '/home/g/Asteroids/input/Itokawa.obj' # Example: Replace with your file

    N_LATITUDE = 1000  # Number of latitude samples
    N_LONGITUDE =1000  # Number of longitude samples

    # --- Processing ---
    mesh = load_asteroid_mesh(asteroid_filepath)
    
    rotation = False #make it true for Itokawa
    if rotation:

        angle_rad = np.pi / 2
        rotation_axis = [0,1,0]
        # Create the 4x4 transformation matrix
        rotation_matrix = trimesh.transformations.rotation_matrix(angle=angle_rad, direction=rotation_axis)
       
       # 3. Apply the transformation to the mesh
        mesh.apply_transform(rotation_matrix)

    if mesh:
        # 1. Fit the best-fit sphere
        # Using vertices is common, could use face midpoints or sample points too
        points_for_fitting = mesh.vertices
        center, radius = fit_sphere_least_squares(points_for_fitting)

        # 2. Create the spherical grid based on the fitted sphere
        lats_rad, lons_rad, grid_points = create_spherical_grid(center, radius, N_LATITUDE, N_LONGITUDE)

        # 3. Calculate deviations from the sphere to the actual surface
        start_time =time.time()
        deviations_raw = calculate_surface_deviations(mesh, center, radius, grid_points, lats_rad, lons_rad)
        deviations = interpolate_nan_deviations(lats_rad, lons_rad, deviations_raw, method='linear')
        print("The time taken for creating rays: ", time.time()-start_time)
        
        # 4. Report basic statistics
        valid_deviations = deviations[~np.isnan(deviations)]
        if len(valid_deviations) > 0:
            print("\nDeviation Statistics:")
            print(f"  Min Deviation: {np.nanmin(deviations):.4f}")
            print(f"  Max Deviation: {np.nanmax(deviations):.4f}")
            print(f"  Mean Deviation: {np.nanmean(deviations):.4f}")
            print(f"  Std Dev: {np.nanstd(deviations):.4f}")
            print(f"  Points with data: {len(valid_deviations)} / {deviations.size}")
        else:
            print("\nNo valid deviation data could be calculated.")

        # 5. Visualize the deviation map 
        plot_deviation_map(lats_rad, lons_rad, deviations,
                           title=f"Asteroid Surface Deviation from Best-Fit Sphere\n(File: {asteroid_filepath.split('/')[-1]})")

        # Optional: Save deviations to a file (e.g., NumPy array or CSV)
        # np.save('asteroid_deviations.npy', deviations)
        # Or save lats, lons, deviations together
        # np.savez('asteroid_topo_data.npz', lats=lats_rad, lons=lons_rad, deviations=deviations)
        print("\nProcessing complete.")
        save_asteroid_data("/home/g/Asteroids/input/Spherical.npz", radius, lats_rad, lons_rad, deviations) 

    else:
        print("Could not load mesh. Exiting.")