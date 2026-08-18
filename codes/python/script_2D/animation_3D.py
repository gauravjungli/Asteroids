#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Apr 12 15:50:14 2025

@author: g
"""

import numpy as np
import pyvista as pv
import math # For pi
from IO import load_asteroid_data, Output_File, extract_number
import os
import time
from scipy.interpolate import griddata
import pdb
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors # For custom colormap
from PIL import Image, ImageDraw

os.environ['IMAGEIO_FFMPEG_EXE'] = '/home/g/miniconda3/envs/G/bin/ffmpeg'

def regrid(lat_low_2d,lon_low_2d,height_low_2d,lat_high_2d,lon_high_2d):
    
    points_low_res = np.vstack((lat_low_2d.flatten(),lon_low_2d.flatten())).T
    points_high_res = np.vstack((lat_high_2d.flatten(), lon_high_2d.flatten())).T
    nlat_high = lat_high_2d.shape[0]
    nlon_high =  lat_high_2d.shape[1]
    values_low_res = height_low_2d.flatten()

    # --- 7. Perform Interpolation ---
    interpolation_method = 'cubic' # Options: 'nearest', 'linear', 'cubic'


    start_time = time.time()
    interpolated_values_flat = griddata(
        points_low_res,      # Coordinates of the known low-res data points (N, 2)
        values_low_res,      # Values at the known low-res data points (N,)
        points_high_res,     # Coordinates where we want to interpolate (M, 2)
        method=interpolation_method,
        fill_value=0    # Use NaN for points outside the low-res domain
    )
    end_time = time.time()
    print(f"Interpolation finished in {end_time - start_time:.2f} seconds.")


    # --- 8. Reshape Output ---
    height_high_2d = interpolated_values_flat.reshape((nlat_high, nlon_high))


    # --- Optional: Verification/Visualization ---
    # Check how many non-NaN values were generated
    num_valid_points = np.count_nonzero(~np.isnan(height_high_2d))
    print(f"\nNumber of valid (non-NaN) points in interpolated high-res grid: {num_valid_points}")
    if num_valid_points == 0:
        print("WARNING: Interpolation resulted in all NaNs. Check coordinate systems, ranges, and data.")
    elif num_valid_points < points_high_res.shape[0]:
         print("INFO: Fewer valid points than high-res points.")


    # Optional: Save the resulting high-resolution grid
    # np.save('high_res_interpolated_data.npy', data_high_res_interpolated)
    # np.savetxt('high_res_interpolated_data.csv', data_high_res_interpolated, delimiter=',') # Can be very large!

    # Optional: Simple plot using matplotlib
    # comment this line if you want to see the plots of regridding
    return height_high_2d

    try:

        plt.figure(figsize=(10, 5))
        # Use pcolormesh for better handling of grid cell boundaries
        plt.pcolormesh(lon_low_2d, lat_low_2d, height_low_2d, shading='auto', cmap='viridis')
        # Add low-res grid points for context
       # plt.scatter(lon_low_2d, lat_low_2d, s=10, c='red', label='Low-res points')
        plt.colorbar(label='Interpolated Fluid Value')
        plt.xlabel('Longitude')
        plt.ylabel('Latitude')
        plt.title(f'Interpolated Data ({interpolation_method}) on Low-Resolution Grid')
        # Limit plot view to the approximate low-res area for clarity
        plt.xlim(0, 2*np.pi)
        plt.ylim(0, np.pi)
        plt.show()

        plt.figure(figsize=(10, 5))
         # Use pcolormesh for better handling of grid cell boundaries
        plt.pcolormesh(lon_high_2d, lat_high_2d, height_high_2d, shading='auto', cmap='viridis')
         # Add low-res grid points for context
        # plt.scatter(lon_low_2d, lat_low_2d, s=10, c='red', label='Low-res points')
        plt.colorbar(label='Interpolated Fluid Value')
        plt.xlabel('Longitude')
        plt.ylabel('Latitude')
        plt.title(f'Interpolated Data ({interpolation_method}) on High-Resolution Grid')
         # Limit plot view to the approximate low-res area for clarity
        plt.xlim(0, 2*np.pi)
        plt.ylim(0, np.pi)
        plt.show()

    except Exception as e:
        print(f"\nError during plotting: {e}")


    return height_high_2d

def create_3D_mesh(radius,lats,lons,devs):
    # Calculate actual radius at each lat/lon point
    actual_radius = radius +devs  
    
    
    # Calculate Cartesian coordinates using the *actual* radius
    x = actual_radius * np.sin(lats) * np.cos(lons)
    y = actual_radius * np.sin(lats) * np.sin(lons)
    z = actual_radius * np.cos(lats)
    
    print(f"Coordinate shapes (X, Y, Z): {x.shape}, {y.shape}, {z.shape}")
    
    # Assign the coordinates (needs to be flattened in Fortran 'F' order or reshaped correctly)
    # Or simply stack them column-wise for PyVista (N x 3 shape)
    points_flat = np.vstack((x.ravel(order='C'),   # Using 'C' order matching numpy default
                        y.ravel(order='C'),   # if data[i, j] is lat[i], lon[j]
                        z.ravel(order='C'))).T
    
    # Also return the non-flattened 3D coordinates array for easier indexing
    # Output shape: (n_lat, n_lon, 3)
    coords_xyz = np.stack((x, y, z), axis=-1)

    return points_flat, coords_xyz


#Used for plotting latitude and longitude lines on the sphere


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

def setup_old_plotter(radius,grid,initial_coords_xyz):
    # Set up the plotter
    # Use off_screen=True if you *only* want to save the animation file
    # Use off_screen=False (or omit) if you want an interactive window first
    vmin = 0
    vmax = 1
    clim = [vmin, vmax]
    
    # --- Define Custom Colormap ---
    # Simple map: grey at the start (vmin=0), yellow at the end (vmax)
    nodes = [0.0, 0.5, 1.0]
# Corresponding colors at the nodes
    colors = ["blue","green","yellow"]

# Create the colormap using from_list with (node, color) pairs
    grey_red_yellow_cmap = mcolors.LinearSegmentedColormap.from_list(
    "GreyRedYellowMap", list(zip(nodes, colors)))
    
    grey_yellow_cmap = mcolors.LinearSegmentedColormap.from_list(
        "grey_yellow_map", ["grey", "blue"]
    )
    plotter = pv.Plotter(off_screen=True, window_size=[1000, 800])
    camera_distance_factor = 3.0    # How many radii away the camera should be (adjust)
    min_camera_distance = 2.0       # Minimum distance to prevent camera getting too close
    
    # Add the mesh to the plotter
    # 'scalars' tells PyVista which point data array to use for coloring
    # 'clim' sets the color limits
    # 'scalar_bar_args' customizes the color bar
    
    graticule_pil_image = create_graticule_texture(
       width=300, height=200, n_lat_lines=10, n_lon_lines=15, line_color="blue")

   ## *** CORRECTED STEP: Convert PIL Image to NumPy array ***
    graticule_numpy_array = np.array(graticule_pil_image)

    ##Create the PyVista texture from the NumPy array
    graticule_texture = pv.Texture(graticule_numpy_array)
    
# --- Define Scalar Bar Customization ---
    scalar_bar_args = {
    'title': " Height (m)",  # Set your desired label text here
    'title_font_size': 30,           # Font size for the title/label
    'label_font_size': 30,           # Font size for the numeric tick labels
    'position_x': 0.88,              # Horizontal position of bottom-left corner (0=left, 1=right)
    'position_y': 0.30,              # Vertical position of bottom-left corner (0=bottom, 1=top)
    #'width': 0.9,                    # Relative width of the color bar (increase for horizontal)
   # 'height': 0.08,                   # Relative height of the color bar (decrease for horizontal)
    'vertical': True,               # Set to False for a horizontal bar
    'n_labels': 5,                   # Approximate number of numeric labels
    'fmt': "%.1f",                   # Format string for numeric labels (e.g., 1 decimal place)
    'color': 'black',                # Color of the text labels and title
    'shadow': True,                  # Add a shadow to the text for better visibility
    'interactive': False,            # Set to True if you want to drag the bar (usually False for animations)
    }
    
    plotter.add_mesh(grid,
                     scalars='fluid_height',
                     cmap=grey_red_yellow_cmap,
                     clim=clim,
                     #texture=graticule_texture, #it is just for a sphere
                     #smooth_shading=True,
                     scalar_bar_args = scalar_bar_args,
                     name='surface') # Give the actor a name for easy updates
    
    plotter.camera.zoom(1.5)
    plotter.background_color = 'white' # Set background color
    
    
    # plotter.camera_position = 'iso' # Or a custom (x, y, z)
    # plotter.camera.focal_point = (0, 0, 0)
    # plotter.camera.distance = radius * 4.0
    # === 3. Update Camera ===
    # Find the 2D indices (lat_idx, lon_idx) of the max height in the current frame
    # Using argmax on the 2D array is slightly easier than on flattened
    #max_indices_2d = np.unravel_index(np.argmax(height_new), height_new.shape)
    lat_idx, lon_idx = (int(50),int(100)) # Shape is (n_lat, n_lon), so result is (lat_idx, lon_idx)

    # Get the 3D coordinate of the peak on the *current* surface using the 3D coord array
    peak_coord_3d = initial_coords_xyz[lat_idx, lon_idx, :] # Shape (3,)

    # Set the camera's focal point (where it looks) to the peak location
    plotter.camera.focal_point = peak_coord_3d

    # Calculate the direction vector from the origin to the peak
    # This determines the line along which the camera will be placed
    vector_to_peak = peak_coord_3d - [0, 0, 0] # Vector from origin
    norm_vector_to_peak = np.linalg.norm(vector_to_peak)

    if norm_vector_to_peak > 1e-9: # Avoid division by zero if peak is exactly at origin
        unit_vector_to_peak = vector_to_peak / norm_vector_to_peak
    else:
        unit_vector_to_peak = np.array([0, 0, 1.0]) # Default view direction if peak at origin

    # Calculate desired camera distance based on object's current max radius
    #current_max_radius = radius + np.max(devs)
    desired_distance = max(min_camera_distance, camera_distance_factor * radius)


    # Adjust the camera clipping range based on the new position/view
    plotter.reset_camera_clipping_range()
    return plotter

def setup_plotter(radius, grid, initial_coords_xyz):
    vmin = 0
    vmax = 1
    clim = [vmin, vmax]
    
    # --- Define Custom Colormap ---
    nodes = [0.0, 0.5, 1.0]
    colors = ["grey", "green", "yellow"] # Note: kept your colors, but renamed variable to match
    grey_green_yellow_cmap = mcolors.LinearSegmentedColormap.from_list(
        "GreyGreenYellowMap", list(zip(nodes, colors))
    )
    
    plotter = pv.Plotter(off_screen=False, window_size=[1000, 800])
    camera_distance_factor = 3.0    
    min_camera_distance = 2.0       
    
    # --- Define Scalar Bar Customization ---
    scalar_bar_args = {
        'title': " Height (m)",  
        'title_font_size': 30,            
        'label_font_size': 30,            
        'position_x': 0.88,               
        'position_y': 0.30,               
        'vertical': True,               
        'n_labels': 5,                   
        'fmt': "%.1f",                   
        'color': 'black',                
        'shadow': True,                  
        'interactive': False,            
    }
    
    # 1. Add the main fluid surface (Removed the texture argument)
    plotter.add_mesh(grid,
                     scalars='fluid_height',
                     cmap=grey_green_yellow_cmap,
                     clim=clim,
                     scalar_bar_args=scalar_bar_args,
                     name='surface')
                     
    # === NEW: Add Latitude and Longitude Lines ===
    # Create a sphere just slightly larger than the base radius to prevent z-fighting
    lat_lon_grid = pv.Sphere(
        radius=radius * 1.01, 
        theta_resolution=15, # Number of longitude lines
        phi_resolution=10    # Number of latitude lines
    )
    
    # Add it to the plotter as a wireframe
    plotter.add_mesh(lat_lon_grid, 
                     style='wireframe', 
                     color='blue', 
                     line_width=1.5,
                     name='graticule_lines')

    plotter.camera.zoom(1.5)
    plotter.background_color = 'white' 
    
    # === Update Camera ===
    lat_idx, lon_idx = (int(50), int(100)) 
    peak_coord_3d = initial_coords_xyz[lat_idx, lon_idx, :] 

    plotter.camera.focal_point = peak_coord_3d

    vector_to_peak = peak_coord_3d - [0, 0, 0] 
    norm_vector_to_peak = np.linalg.norm(vector_to_peak)

    if norm_vector_to_peak > 1e-9: 
        unit_vector_to_peak = vector_to_peak / norm_vector_to_peak
    else:
        unit_vector_to_peak = np.array([0, 0, 1.0]) 

    desired_distance = max(min_camera_distance, camera_distance_factor * radius)
    
    # === CORRECTED STEP: Actually apply the camera position ===
    # Set camera position along the vector away from the focal point
    plotter.camera.position = peak_coord_3d + (unit_vector_to_peak * desired_distance)

    plotter.reset_camera_clipping_range()
    
    return plotter

def save_3D_plot(target):
    #pdb.set_trace()  
    
    #file1 = file+'/run1/base.txt'
    file1 = Output_File(target,'output',['base.txt'])
    epsilon = target.epsilon
    radius = target.d/2
    n_lat = target.x_res
    n_lon = target.y_res
   

    if os.path.exists(file1):
       
        w = np.loadtxt(file1,delimiter=",",dtype='float')
        
        theta = w[:,0]
        phi = w[:,1]
        base = epsilon*w[:,2]

        r = w[:,4]
        dr = w[:,5]

        metric = r**2 + dr**2
         
        rad = np.sqrt((2*r**2*np.sqrt(metric)*(base) + metric*((base)**2 + r**2  ))/metric)
        print("maximum radial distance is ", max(rad) )
        z = np.cos(theta)*r + 1/np.sqrt(metric)*(base*dr*np.sin(theta)+np.cos(theta)*r*base)
        theta_new = np.arccos(z/rad)
        
        x = radius*rad*np.sin(theta_new)*np.cos(phi)
        y = radius*rad*np.sin(theta_new)*np.sin(phi)
        z = radius*rad*np.cos(theta_new)
        

    initial_points = np.vstack((x,   y,  z)).T
    
    initial_coords_xyz = np.stack((x.reshape((n_lat,n_lon)), y.reshape((n_lat,n_lon)), z.reshape((n_lat,n_lon))), axis=-1)
    
    lons = phi.reshape((n_lat,n_lon))

    grid = pv.StructuredGrid()
    #grid.texture_map_to_sphere(inplace=True)  #it is just for a sphere
    grid.dimensions =(n_lat,n_lon,1)
    grid.points = initial_points
    initial_height = np.full_like(lons,0)
    initial_height = initial_height.ravel(order='C')
    grid.point_data['fluid_height'] = initial_height
    grid.point_data.active_scalars_name = 'fluid_height' # Explicitly set active scalar
    
    plotter =  setup_old_plotter(radius,grid,initial_coords_xyz)
    
    main_dir = Output_File(target,'output',['data','landslides_1'])
    

    # Get and sort files by number within the subdirectory
    dirFiles = sorted( [f for f in os.listdir(main_dir) if f.lower() != "log.txt"], key=extract_number)
    
        
    print("Starting animation generation (saving to MP4)...")
        # Open a movie file (requires ffmpeg)
    try:
        movie = Output_File(target,'output',['animation.mp4'])
        plotter.open_movie(movie) # Adjust framerate as needed

        frame_dir = Output_File(target,'output',['frames'])
        if not os.path.exists(frame_dir):
            os.makedirs(frame_dir)
        i=0
        for file in dirFiles:
            if i%1==0:
                myfile = os.path.join(main_dir,file)
                w=np.loadtxt(myfile,delimiter=",",dtype=float)
                
                theta = w[:,0]
                phi = w[:,1]
                base = epsilon*w[:,2]
                height =  epsilon*w[:,3]
                r = w[:,4]
                dr = w[:,5]
    
                metric = r**2 + dr**2
                 
                rad = np.sqrt((2*r**2*np.sqrt(metric)*(base) + metric*((base)**2 + r**2  ))/metric)
                print(file )
                z = np.cos(theta)*r + 1/np.sqrt(metric)*(base*dr*np.sin(theta)+np.cos(theta)*r*base)
                theta_new = np.arccos(z/rad)
                
                x = radius*rad*np.sin(theta_new)*np.cos(phi)
                y = radius*rad*np.sin(theta_new)*np.sin(phi)
                z = radius*rad*np.cos(theta_new)
                
                current_points = np.vstack((x,   y,  z)).T
                
                
                vmax = np.nanmax(height*radius)
                clim = [0, vmax]
                
                max_idx = np.argmax(height)
                current_peak_pos = current_points[max_idx]
    
                # Point the camera at the new peak location
                plotter.camera.focal_point =(0,0,0)
                
                # Update scalars on the mesh
                actor = plotter.actors['surface']
                mesh_in_plotter = actor.mapper.dataset
                mesh_in_plotter.points = current_points
                plotter.camera.elevation = -25 #-i*0.1
                plotter.camera.azimuth = 135 #- min(100,i*1.0)
                plotter.camera.distance = radius * 2.0
               # plotter.camera.zoom(1.25)
                plotter.render()
                
                mesh_in_plotter.point_data['fluid_height'] = height*radius
                mesh_in_plotter.point_data.active_scalars_name = 'fluid_height'
                actor.mapper.scalar_range = clim
           
            
            
                print('saved frames:',i)
                filename = os.path.join(frame_dir, f"frame_{i:04d}.png")
                plotter.screenshot(filename, transparent_background=True)
            i=i+1
            # Write the current view as a frame to the movie file
        #    plotter.write_frame()

# ... rest of the loop ...
    
        # Finalize and close the movie file and plotter
        plotter.close()
        print("Animation saved successfully to Bennu.mp4")
    
    except ImportError:
        print("\nERROR: Failed to save animation.")
        print("Saving MP4 requires 'imageio' and 'imageio-ffmpeg'.")
        print("Install them using: pip install imageio imageio-ffmpeg")
    except FileNotFoundError:
        print("\nERROR: Failed to save animation.")
        print("Could not find ffmpeg.")
        print("Please ensure ffmpeg is installed and accessible in your system's PATH.")
    except Exception as e:
        print(f"\nAn error occurred during animation saving: {e}")
    
#%%
def save_3D_plot_evol(target):
    #pdb.set_trace()  
    
    #file1 = file+'/run1/base.txt'
    file1 = Output_File(target,'output',['base.txt'])
    epsilon = target.epsilon
    radius = target.d/2
    n_lat = target.x_res
    n_lon = target.y_res
   

    if os.path.exists(file1):
       
        w = np.loadtxt(file1,delimiter=",",dtype='float')
        
        theta = w[:,0]
        phi = w[:,1]
        base = epsilon*w[:,2]

        r = w[:,4]
        dr = w[:,5]

        metric = r**2 + dr**2
         
        rad = np.sqrt((2*r**2*np.sqrt(metric)*(base) + metric*((base)**2 + r**2  ))/metric)
        print("maximum radial distance is ", max(rad) )
        z = np.cos(theta)*r + 1/np.sqrt(metric)*(base*dr*np.sin(theta)+np.cos(theta)*r*base)
        theta_new = np.arccos(z/rad)
        
        x = radius*rad*np.sin(theta_new)*np.cos(phi)
        y = radius*rad*np.sin(theta_new)*np.sin(phi)
        z = radius*rad*np.cos(theta_new)
        

    initial_points = np.vstack((x,   y,  z)).T
    
    initial_coords_xyz = np.stack((x.reshape((n_lat,n_lon)), y.reshape((n_lat,n_lon)), z.reshape((n_lat,n_lon))), axis=-1)
    
    lons = phi.reshape((n_lat,n_lon))

    grid = pv.StructuredGrid()
    #grid.texture_map_to_sphere(inplace=True)  #it is just for a sphere
    grid.dimensions =(n_lat,n_lon,1)
    grid.points = initial_points
    initial_height = np.full_like(lons,0)
    initial_height = initial_height.ravel(order='C')
    grid.point_data['fluid_height'] = initial_height
    grid.point_data.active_scalars_name = 'fluid_height' # Explicitly set active scalar
    
    plotter =  setup_old_plotter(radius,grid,initial_coords_xyz)
    
    main_dir = Output_File(target,'output',['data'])
    

    # Get and sort files by number within the subdirectory
    dirFiles = sorted( [f for f in os.listdir(main_dir) if f.lower().endswith('.csv')], key=extract_number)
    
        
    print("Starting animation generation (saving to MP4)...")
        # Open a movie file (requires ffmpeg)
    try:
        movie = Output_File(target,'output',['animation.mp4'])
        plotter.open_movie(movie) # Adjust framerate as needed

        frame_dir = Output_File(target,'output',['frames'])
        if not os.path.exists(frame_dir):
            os.makedirs(frame_dir)
        i=0
        for file in dirFiles:
            if i%1==0:
                myfile = os.path.join(main_dir,file)
                w=np.loadtxt(myfile,delimiter=",",dtype=float)
                
                theta = w[:,0]
                phi = w[:,1]
                base = epsilon*w[:,2]
                height =  epsilon*w[:,3]
                r = w[:,4]
                dr = w[:,5]
    
                metric = r**2 + dr**2
                 
                rad = np.sqrt((2*r**2*np.sqrt(metric)*(base) + metric*((base)**2 + r**2  ))/metric)
                print(file )
                z = np.cos(theta)*r + 1/np.sqrt(metric)*(base*dr*np.sin(theta)+np.cos(theta)*r*base)
                theta_new = np.arccos(z/rad)
                
                x = radius*rad*np.sin(theta_new)*np.cos(phi)
                y = radius*rad*np.sin(theta_new)*np.sin(phi)
                z = radius*rad*np.cos(theta_new)
                
                current_points = np.vstack((x,   y,  z)).T
                
                v = ((r+height+base)*radius-radius)
                vmax = np.nanmax(v)
                vmin = np.nanmin(v)
                clim = [vmin, vmax]
                
                max_idx = np.argmax(height)
                current_peak_pos = current_points[max_idx]
    
                # Point the camera at the new peak location
                plotter.camera.focal_point =(0,0,0)
                
                # Update scalars on the mesh
                actor = plotter.actors['surface']
                mesh_in_plotter = actor.mapper.dataset
                mesh_in_plotter.points = current_points
                plotter.camera.elevation = -25 #-i*0.1
                plotter.camera.azimuth = 135 #- min(100,i*1.0)
                plotter.camera.distance = radius * 2.0
               # plotter.camera.zoom(1.25)
                plotter.render()
                
                mesh_in_plotter.point_data['fluid_height'] = v
                mesh_in_plotter.point_data.active_scalars_name = 'fluid_height'
                actor.mapper.scalar_range = clim
           
            
            
                print('saved frames:',i)
                filename = os.path.join(frame_dir, f"frame_{i:04d}.png")
                plotter.screenshot(filename, transparent_background=True)
            i=i+1
            # Write the current view as a frame to the movie file
        #    plotter.write_frame()

# ... rest of the loop ...
    
        # Finalize and close the movie file and plotter
        plotter.close()
        print("Animation saved successfully to Bennu.mp4")
    
    except ImportError:
        print("\nERROR: Failed to save animation.")
        print("Saving MP4 requires 'imageio' and 'imageio-ffmpeg'.")
        print("Install them using: pip install imageio imageio-ffmpeg")
    except FileNotFoundError:
        print("\nERROR: Failed to save animation.")
        print("Could not find ffmpeg.")
        print("Please ensure ffmpeg is installed and accessible in your system's PATH.")
    except Exception as e:
        print(f"\nAn error occurred during animation saving: {e}")
    
