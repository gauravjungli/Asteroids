#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 17 12:33:42 2025

@author: g
"""

import pyvista as pv
import numpy as np
from PIL import Image, ImageDraw

# --- 1. Create the 3D Shape ---
# Example: A section of a sphere
radius = 1.0
lat_min, lat_max = np.radians(-60), np.radians(60)
lon_min, lon_max = np.radians(-120), np.radians(120)
n_lat, n_lon = 50, 100  # Grid resolution

# Create latitude and longitude arrays
lat = np.linspace(lat_min, lat_max, n_lat)
lon = np.linspace(lon_min, lon_max, n_lon)
lon_grid, lat_grid = np.meshgrid(lon, lat)

# Convert spherical coordinates to Cartesian coordinates
x = radius * np.cos(lat_grid) * np.cos(lon_grid)
y = radius * np.cos(lat_grid) * np.sin(lon_grid)
z = radius * np.sin(lat_grid)

# Create a PyVista StructuredGrid
mesh = pv.StructuredGrid(x, y, z)

# --- 2. Generate Texture Coordinates (UV) ---
# Method 1: Calculate manually (more general if not a perfect sphere section)
# u_coords = (lon_grid - lon_min) / (lon_max - lon_min)
# v_coords = (lat_grid - lat_min) / (lat_max - lat_min)
# mesh.active_texture_coordinates = np.stack((u_coords.ravel(), v_coords.ravel()), axis=1)

# Method 2: Use PyVista's built-in sphere mapping (often better for spheres)
# This assumes your mesh is indeed part of a sphere centered at the origin
mesh.texture_map_to_sphere(inplace=True)
# Note: If using this method, the texture image should ideally cover the
# full 0-360 deg longitude and -90 to +90 deg latitude range for correct mapping.
# Our texture creation below will be simpler and match the *range* of the mesh instead.
# For a full sphere texture map, you'd create a texture representing the whole globe.

# --- 3. Create the Graticule Texture Image ---


# --- (Previous code: Imports, Shape Creation, UV Coords) ---
# ... (Keep the code for shape creation and texture coordinates as before) ...

# --- 3. Create the Graticule Texture Image ---

def create_graticule_texture(
    width=500,
    height=250,
    n_lat_lines=7,
    n_lon_lines=13,
    line_color="gray",
    bg_color="white",
):
    """Creates a PIL Image with latitude and longitude lines."""
    img = Image.new("RGB", (width, height), color=bg_color)
    draw = ImageDraw.Draw(img)

    # Draw longitude lines (vertical)
    lon_spacing = width / (n_lon_lines - 1)
    for i in range(n_lon_lines):
        x = int(i * lon_spacing)
        x = max(0, min(width - 1, x)) # Clamp to bounds
        draw.line([(x, 0), (x, height - 1)], fill=line_color, width=1)

    # Draw latitude lines (horizontal)
    lat_spacing = height / (n_lat_lines - 1)
    for i in range(n_lat_lines):
        y = int(i * lat_spacing)
        y = max(0, min(height - 1, y)) # Clamp to bounds
        draw.line([(0, y), (width - 1, y)], fill=line_color, width=1)

    return img # Return the PIL Image

# Create the texture image using PIL
graticule_pil_image = create_graticule_texture(
    width=600, height=300, n_lat_lines=10, n_lon_lines=15, line_color="blue"
)

# *** CORRECTED STEP: Convert PIL Image to NumPy array ***
graticule_numpy_array = np.array(graticule_pil_image)

# Create the PyVista texture from the NumPy array
graticule_texture = pv.Texture(graticule_numpy_array)

# --- 4. Apply Texture and Visualize ---
plotter = pv.Plotter(window_size=[800, 800])

# Add the mesh to the plotter and assign the texture
actor = plotter.add_mesh(mesh, texture=graticule_texture, smooth_shading=True)

# Optional: Add axes, set background, camera position etc.
plotter.show_axes()
plotter.background_color = 'black'
plotter.camera_position = 'xy'
plotter.camera.zoom(1.5)

print("Displaying the shape with latitude and longitude lines (graticule texture).")
plotter.show()