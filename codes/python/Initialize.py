#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Mar 29 16:59:36 2025

@author: g
"""
from IO import Output_File, Exparameter, load_asteroid_data
import os
import subprocess
from collisions import G
import math
import shutil
import numpy as np
import time
from scipy.interpolate import RegularGridInterpolator
import matplotlib.pyplot as plt
import pdb
from gravity import Gravitycalc_2D, Gravitycalc
from Fit import Fit
from scipy.ndimage import gaussian_filter1d
#%%
"""
This is the first function that is called. It initializes multiple simulations. It is called by GUI at the 
beginning of the simulation. It first create a list named run. Each run corresponds to the number of simulation 
to be run. For each run it creates a parameters list. Then it creates the output directories. It exports the 
parameters.txt file for each run and also returns the parameters_list required for GUI.
"""

def Initialize_simulations(parameters,parameters_list=[]):  
    
    run = list(range(1, int(parameters['Number of simulations'])+1))  # List of input values
    for i in run:
        new_parameters={}
        for key in parameters:
            new_parameters[key]=parameters[key]
        new_parameters['run']=i
        mydir=Output_File (new_parameters,"output")
        if os.path.exists(mydir):
            subprocess.run(["rm", "-r", mydir])
        else:
            print("No directory exists")
        os.makedirs(mydir, exist_ok=True) 
        new_parameters['Data folder']=os.path.join(mydir,"data")
        os.makedirs(new_parameters['Data folder'], exist_ok=True) 
        Exparameter(new_parameters)
        parameters_list.append(new_parameters)
    
    Exparameter(parameters)


#%%
"""
    This initializes the simulation before each landslide. Following steps are done:
        1. Non-dimensionalizing the inertia, omega
        2. setting the no. of landslides count to zero
        3. Intializing the time to zero
        4. Initializing the dia variable to the input diameter
        5. Initializing the epsilon to zero for explicit case
        6. Export these changes to parameters file
        7. Write the initial data to the output.yorp file
        8. Initialize the grid and base
        9. Calculate the gravity
        10.Copy the executable file to the local location
        
"""

def Initialize(parameters,target):
    

    
    parameters['Maximum acceleration'] = 0
    parameters['slides']    = 0
    parameters['time']      = 0
    parameters["k_d"] = target.k_d/(G * (4/3) * math.pi * target.dens)**0.5
    parameters['omega']     = target.omega[2]/(G * 4 / 3 * math.pi * target.dens) ** 0.5   
    parameters['epsilon']   = float(parameters['epsilon'])   
    #parameters['Seismic_shaking_time'] = target.t_lan/(target.d/2/target.grav)**0.5
    parameters['Mass shed'] = 0
    parameters['Initial mass'] = 0
    mydir = Output_File(parameters,"output")
    
    
    file = Output_File(parameters, "output", ["output.yorp"])
    with open(file, "w") as output_file:  # Overwrites existing files
        output_file.write(
            f"{0 :12.6e} {target.omega[2]:12.8e} {target.obliq:12.8e}\n")
        
     
    grid_uniform(parameters,target)
    
    if target.landslide:
        
        dimension = parameters['Dimension'] 
              
        


        executable_file = Output_File(parameters,"codes",['build',f'landslides_{dimension}',parameters["executable"]])
        
        try:
            shutil.copy(executable_file, mydir)
            print("Files copied successfully.")
        except FileNotFoundError:
            print("Source/executable file not found.")
        except PermissionError:
            print("Permission denied.")
        except Exception as e:  
            print("An error occurred:", e)
     
    parameters['jinertia1'] = target.jinertia[0] / (target.d / 2)**5 / target.dens
    parameters['jinertia']  = target.jinertia[2] / (target.d / 2)**5 / target.dens
    parameters['Current diameter'] = target.d
    Exparameter(parameters)
    
#%%


class Crater_Impactor:
    
    def __init__(self,explicit=True):
        if explicit:
            self.phi   = math.acos(1 - 2 * np.random.random()) / 2
            self.vel   = 5.5e+3#scipy.stats.maxwell.ppf(random.random(), scale=3.232)*1000
            self.d     = 5 
            self.theta = 2 * math.pi * np.random.random()
            self.Phi = np.pi
            self.Theta   = np.pi/3
        
            self.dens  = 1500

        self.M         = (math.pi / 6) * self.dens * self.d**3
        self.explicit  = explicit
        

def grid_uniform(parameters,target):
    
    def dunes(x,y):
        x_peak = np.deg2rad(45)
        y_peak = np.deg2rad(180)
        return np.exp(-((x - x_peak)**2 + (y - y_peak)**2) / 0.05)
    
    def crater(): #Not required. Will have to find a way to include
        nonlocal mydir,parameters,target
        impactor = Crater_Impactor()
        
     #   base = Crater(parameters,target,impactor)
         
        i=0
        with open(mydir+"/base.txt", "w") as file:
            for x in x_values:
                for y in y_values:
        
                    file.write(f"{x:.12f},{y:.12f},{base[i,2]:.12f},1\n")  # Format to 6 decimal places for precision
                    i=i+1
                    
        print("Grid data with crater height profile saved to base.txt")
        
        
    mydir = Output_File(parameters,"output")

    if parameters['Dimension'] == '2D':
        

    # Define the grid ranges
 
        x_values, y_values, r, dr, ddr, base = grid(parameters,target) 
                
        base = base/float(parameters['epsilon'])

    # Open file to write the grid data
        i=0
        j=0
        with open(mydir+"/base.txt", "w") as file:
            for x in x_values:
                for y in y_values:
                    
                    height = dunes(x,y) #change

                    file.write(f"{x:.12f},{y:.12f},{base[i]:.12f},{height:.12f},{r[j]:.12f},{dr[j]:.12f},{ddr[j]:.12f}\n") 
                    i+=1
                j+=1
        print("Grid data with Bennu basal profile and a dune saved to base.txt")
        
        #crater()
 

    if parameters['Dimension'] == '1D':
        #pdb.set_trace() 
        res = int(parameters["Resolution"])
        dia = target.d


        base = np.zeros((res,5))
        offset=float(parameters["offset"])
        dx=(math.pi-2*offset)/res 
        parameters['dx'] =dx
        #c = 500   #uncomment for single run with given shape
        #a = 1000
        #D = (a**2*c)**(1/3) 
        #parameters['Current diameter'] = D
        for i in range(res):
            
            base[i,0] = offset + dx * (i+ 0.5)   
           # base [i,2] =1 
        axisymmetric_asteroid(parameters,target,base)
        parameters["Initial mass"] =  2*np.pi/3*np.trapezoid(((base[2:res-2,1])**3*
                                                   np.sin(base[2:res-2,0])),base[2:res-2,0])*(dia/2)**3
        #base = Fit(parameters,base_old=base)
        #axisymmetric(base,a/2,c/2,D/2)
        
        np.savetxt(mydir+"/base.txt",base,delimiter=",")
        
        target.rgrav,target.tgrav =  Gravitycalc(parameters,target) 
        
def grid_2D(parameters):
    
    nx, ny = int(parameters['X Resolution']), int(parameters['Y Resolution'])
    
# Define the grid ranges
    offset =float(parameters["offset"])
    dx=(np.pi-2*offset)/(nx)
    parameters['dx'] = dx
    
    dy=(2*np.pi)/(ny)
    parameters['dy'] = dy
    
    x_values = np.linspace(offset+dx/2, np.pi-offset-dx/2, nx)
    y_values = np.linspace(0+dy/2, 2 * np.pi-dy/2, ny)
    
    return x_values, y_values
    
        
 
def grid(parameters, target):
    #pdb.set_trace()
    if parameters['Dimension'] == '2D':
      #  pdb.set_trace()
        nx, ny = int(parameters['X Resolution']), int(parameters['Y Resolution'])

        asteroid_file = Output_File(parameters,"input",[f"{parameters['Asteroid']}.npz"])
        radius, lats, lons, devs = load_asteroid_data(asteroid_file)
        
        
        
        grid1_lats = lats[:, 0]
        grid1_lons = lons[0, :]
        devs_axi =  np.mean(devs,axis=1)

        
        x_values, y_values = grid_2D(parameters)

        devs_new = regrid_height_profiles(radius,devs, grid1_lats, grid1_lons, x_values, y_values)
        devs_new_axi =  np.mean(devs_new,axis=1)
        devs_new_1 = devs_new -devs_new_axi.reshape(-1,1)
        r_vals = 1 + devs_new_axi/radius
    
        # Calculate first derivative: dr/dtheta
        r_smooth = gaussian_filter1d(r_vals, sigma=2)
        dr = np.gradient(r_smooth, x_values, edge_order=2)
        
        dr_smooth = gaussian_filter1d(dr, sigma=2)
        # Calculate second derivative: d^2r/dtheta^2
        d2r = np.gradient(dr_smooth, x_values, edge_order=2)
       # print(np.max(d2r),np.min(d2r))
        #to accomodate for the fact that the normal to axisymmtric surface and sphere is not aligned
        cos_theta = r_vals/np.sqrt(r_vals**2+dr**2)
        sin_theta = dr/np.sqrt(r_vals**2+dr**2)
        d_surface = np.gradient(devs_new, x_values, axis=0)/(devs_new + radius)
        cos_theta_matrix =  np.tile(cos_theta.reshape(-1,1),(1,ny))
        sin_theta_matrix =  np.tile(sin_theta.reshape(-1,1),(1,ny))
        devs_new_2 = devs_new_1/(cos_theta_matrix + d_surface*sin_theta_matrix)
        
        
        target.d = 2 * radius*1000
        target.rgrav,target.tgrav = Gravitycalc_2D(parameters,grid1_lats,grid1_lons,devs_axi,radius,x_values,y_values,devs_new_axi,r_vals,dr)
        devs_new_non = devs_new_2.ravel()/radius
        
        return x_values, y_values, r_vals, dr, d2r, devs_new_non
    
    
    
def regrid_height_profiles(radius,grid1_data, grid1_lats, grid1_lons, grid2_lats, grid2_lons):
    """
    Regrids height profile data from Grid 1 to Grid 2 using interpolation.
    (This function remains unchanged - it handles the last dimension correctly)

    Args:
        grid1_data (np.ndarray): Data on Grid 1 (num_lat1, num_lon1, profile_length).
                                 *Must be 3D, even if profile_length is 1.*
        grid1_lats (np.ndarray): Latitudes of Grid 1 (sorted, 1D).
        grid1_lons (np.ndarray): Longitudes of Grid 1 (sorted, 1D).
        grid2_spec (dict): Dictionary defining Grid 2.

    Returns:
        tuple: (grid2_data_interp, grid2_lats, grid2_lons)
            grid2_data_interp (np.ndarray): Interpolated data on Grid 2
                                           (num_lat2, num_lon2, profile_length).
            grid2_lats (np.ndarray): Latitudes of Grid 2.
            grid2_lons (np.ndarray): Longitudes of Grid 2.
    """
    print("Starting regridding process...")
    start_time = time.time()

    # --- Check if Grid 2 is a subset of Grid 1 ---
    if (grid2_lats.min() < grid1_lats.min() or grid2_lats.max() > grid1_lats.max() or
            grid2_lons.min() < grid1_lons.min() or grid2_lons.max() > grid1_lons.max()):
        print("Warning: Grid 2 extends beyond the boundaries of Grid 1. "
              "Interpolation outside Grid 1 will result in NaN values "
              "(or the specified fill_value).")

    # --- Prepare target points for interpolation ---
    grid2_lat_mesh, grid2_lon_mesh = np.meshgrid(grid2_lats, grid2_lons,indexing = 'ij')
    points_to_interpolate = np.vstack((grid2_lat_mesh.ravel(), grid2_lon_mesh.ravel())).T
    print(f"Total points to interpolate for Grid 2: {points_to_interpolate.shape[0]}")

    # --- Create the interpolator ---
    # grid1_data must be 3D here, even if last dim is 1

    try:
        interpolator = RegularGridInterpolator(
            points=(grid1_lats, grid1_lons), # Pass the 1D coordinate arrays
            values=grid1_data,               # Pass the 3D data array (..., 1)
            method='linear',
            bounds_error=False,
            fill_value=np.nan
        )
    except ValueError as e:
         print(f"Error creating interpolator. Check grid dimensions and data consistency: {e}")
         print(f"Is grid1_lats strictly increasing? {np.all(np.diff(grid1_lats) > 0)}")
         print(f"Is grid1_lons strictly increasing? {np.all(np.diff(grid1_lons) > 0)}")
         return None, None, None

    # --- Perform interpolation ---
    print("Performing interpolation...")
    grid2_data_flat = interpolator(points_to_interpolate)

    # --- Reshape the result ---
    # Result will be (n_points, 1), reshape to (n_lat2, n_lon2, 1)
    n_lat2=grid2_lats.shape[0]
    n_lon2=grid2_lons.shape[0]
    grid2_data_interp = grid2_data_flat.reshape(n_lat2, n_lon2)
    #plot_deviation_map( grid2_lat_mesh, grid2_lon_mesh, grid2_data_interp )
    #visualization_3D(radius, grid2_lat_mesh, grid2_lon_mesh, grid2_data_interp)
    end_time = time.time()
    print(f"Regridding completed in {end_time - start_time:.2f} seconds.")

    return grid2_data_interp
    
def plot_deviation_map(lats_rad, lons_rad, deviations, title="Asteroid Surface Deviation"):
    """Plots the deviation data as a 2D map using pcolormesh.
    Note: lats_rad_ignored and lons_rad_ignored are no longer directly used
          for plotting but kept for API consistency if needed elsewhere.
    """

    if deviations is None or np.all(np.isnan(deviations)):
        print("Warning: Deviation data is empty or all NaN. Skipping plot.")
        return


    n_lat, n_lon = deviations.shape
    print(f"Plotting deviations with shape (n_lat={n_lat}, n_lon={n_lon})") # Debug print



    # --- Create the plot ---
    fig, ax = plt.subplots(figsize=(10, 5))


    im = ax.pcolormesh(lons_rad, lats_rad, deviations,
                       shading='auto', # 'flat' needs dimensions (M+1, N+1) for X,Y and (M,N) for C
                       cmap='coolwarm',
                       vmin=np.nanmin(deviations), vmax=np.nanmax(deviations))

    ax.set_xlabel("Longitude (degrees)")
    ax.set_ylabel("Latitude (degrees)")
    ax.set_title(title)

    ax.set_xlim(lons_rad[0,0], lons_rad[0,n_lon-1])
    ax.set_ylim(lats_rad[0,0], lats_rad[n_lat-1,0])
    ax.set_xticks(np.linspace(lons_rad[0,0], lons_rad[0,n_lon-1], 7))
    ax.set_yticks(np.linspace(lats_rad[0,0], lats_rad[n_lat-1,0], 7))

    cbar = fig.colorbar(im, ax=ax)
    cbar.set_label("Deviation (Surface Distance - Sphere Radius) km")

    plt.tight_layout()
    plt.show()   
    
    
def axisymmetric_asteroid(parameters,target,base):
    
    asteroid_name = parameters['Asteroid']
    
    if asteroid_name == 'Spherical':
        
        base[:,1] = 1
       # base[:,2] = 1 # It should only be kept for debug run
        
    else:
        
        asteroid_file = Output_File(parameters,"input",[asteroid_name+'.npz'])
        radius, lats, lons, devs = load_asteroid_data(asteroid_file)
        
        grid1_lats = lats[:, 0]
        grid1_lons = lons[0, :]
    
    
        x_values = base[:,0]
        y_values = grid1_lons
        devs_new = regrid_height_profiles(radius, devs, grid1_lats, grid1_lons, x_values, y_values)
        devs_new_axi =  np.mean(devs_new,axis=1)
        
        r_vals = devs_new_axi
    
        
        target.d = 2 * radius*1000
        # Calculate first derivative: dr/dtheta
        dr_dtheta = np.gradient(r_vals, x_values, edge_order=2)
        
        # Calculate second derivative: d^2r/dtheta^2
        d2r_dtheta2 = np.gradient(dr_dtheta, x_values, edge_order=2)
        base[:,1] = 1+r_vals/radius
        base[:,3] = dr_dtheta/radius
        base[:,4] = d2r_dtheta2/radius  
    return 
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    