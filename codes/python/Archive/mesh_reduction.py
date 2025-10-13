#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 10 08:27:06 2025

@author: g
"""

import trimesh
import os
import time
import pyfqmr

def simplify_mesh(input_filepath, output_filepath, reduction_factor=10):
    """
    Loads a mesh, simplifies it by reducing the number of faces,
    and saves the simplified mesh.

    Args:
        input_filepath (str): Path to the input mesh file (OBJ, GLB, etc.).
        output_filepath (str): Path to save the simplified mesh file.
        reduction_factor (float): The target factor by which to reduce the
                                  number of faces (e.g., 10 means 1/10th).
    """
    print(f"Loading mesh from: {input_filepath}...")
    try:
        # Start timer
        start_time = time.time()
        mesh = trimesh.load(input_filepath, force='mesh')
        load_time = time.time() - start_time
        print(f"Mesh loaded in {load_time:.2f} seconds.")

        # Handle potential scene or multi-body meshes
        if isinstance(mesh, trimesh.Scene):
             print("Input was a scene, concatenating geometries...")
             mesh = trimesh.util.concatenate(mesh.dump())
        elif isinstance(mesh, list):
             print("Input was a list of meshes, concatenating...")
             mesh = trimesh.util.concatenate(mesh)

        if not isinstance(mesh, trimesh.Trimesh):
            print(f"Error: Loaded object is not a Trimesh instance after processing. Type: {type(mesh)}")
            return

        original_vertices = len(mesh.vertices)
        original_faces = len(mesh.faces)
        print(f"Original mesh: {original_vertices} vertices, {original_faces} faces.")

        if original_faces == 0:
            print("Error: Mesh has no faces to simplify.")
            return

        # Calculate the target number of faces
        target_faces = max(10, int(original_faces / reduction_factor)) # Ensure at least a few faces
        print(f"Targeting reduction factor: {reduction_factor}x")
        print(f"Target face count: {target_faces} (approx)")

        # --- Perform the Simplification ---
        # Trimesh uses simplify_quadratic_decimation which is generally good
        # It might require the 'pyfqmr' package for best performance/results
        # If pyfqmr is not installed, it might fall back to a slower version.
        # Try installing it: pip install pyfqmr
        print("Starting mesh simplification (this may take time for large meshes)...")
        start_time = time.time()
        # The simplification function directly targets a face count
        simplified_mesh = mesh.simplify_quadric_decimation(face_count = target_faces)
        simplify_time = time.time() - start_time
        print(f"Simplification finished in {simplify_time:.2f} seconds.")

        final_vertices = len(simplified_mesh.vertices)
        final_faces = len(simplified_mesh.faces)
        actual_reduction_factor_v = original_vertices / final_vertices if final_vertices > 0 else float('inf')
        actual_reduction_factor_f = original_faces / final_faces if final_faces > 0 else float('inf')

        print(f"Simplified mesh: {final_vertices} vertices, {final_faces} faces.")
        print(f"Actual vertex reduction: ~{actual_reduction_factor_v:.2f}x")
        print(f"Actual face reduction: ~{actual_reduction_factor_f:.2f}x")


        # --- Save the Simplified Mesh ---
        print(f"Saving simplified mesh to: {output_filepath}...")
        # Make sure the output directory exists
        os.makedirs(os.path.dirname(output_filepath), exist_ok=True)
        
        start_time = time.time()
        simplified_mesh.export(output_filepath)
        save_time = time.time() - start_time
        print(f"Simplified mesh saved in {save_time:.2f} seconds.")

    except ImportError as e:
         print(f"An Import error occurred during simplification: {e}")
    except Exception as e:
        print(f"An error occurred during simplification: {e}")
        import traceback
        traceback.print_exc()


# --- Main Execution ---
if __name__ == "__main__":
    # !!! --- Configuration --- !!!
    # !!! REPLACE with the path to YOUR large OBJ file !!!
    input_obj_file = '/home/g/Asteroids/input/Bennu.obj'
    # Define where to save the simplified file (can be OBJ, GLB, PLY, STL etc.)
    output_simplified_file = '/home/g/Asteroids/input/small_Bennu.obj' # Or e.g., 'simplified_asteroid.glb'

    # Define the desired reduction factor (10 for 1/10th)
    factor = 10.0
    # --- End Configuration ---

    if not os.path.exists(input_obj_file):
        print(f"Error: Input file not found at '{input_obj_file}'")
    else:
        simplify_mesh(input_obj_file, output_simplified_file, reduction_factor=factor)
        print("\nSimplification process complete.")