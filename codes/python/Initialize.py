#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Mar 29 16:59:36 2025

@author: g
"""
from IO import Output_File, Exparameter
import os
import subprocess
from collisions import G
import math
import shutil
import numpy as np
from Crater import Crater
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
    

    parameters['jinertia1'] = target.jinertia[0] / (target.d / 2)**5 / target.dens
    parameters['jinertia']  = target.jinertia[2] / (target.d / 2)**5 / target.dens
    parameters['slides']    = 0
    parameters['time']      = 0
    parameters['omega']     = target.omega[2]/(G * 4 / 3 * math.pi * target.dens) ** 0.5
    parameters['current diameter'] = target.d 
    parameters['epsilon']   = float(parameters['epsilon'])   
    parameters['Seismic_shaking_time'] = target.t_lan/(target.d/2/target.grav)**0.5
    mydir = Output_File(parameters,"output")
    
    Exparameter(parameters)
    file = Output_File(parameters, "output", ["output.yorp"])
    with open(file, "w") as output_file:  # Overwrites existing files
        output_file.write(
            f"{0 :12.6e} {target.omega[2]:12.8e} {target.obliq:12.8e}\n")
        
     
    grid_uniform(parameters,target)
    
    if target.landslide:
        
        dimension = parameters['Dimension'] 
              
       # target.rgrav,target.tgrav =  Gravitycalc(parameters) 

        executable_file=Output_File(parameters,"codes",['build',f'landslides_{dimension}',parameters["executable"]])
        
        try:
            shutil.copy(executable_file, mydir)
            print("Files copied successfully.")
        except FileNotFoundError:
            print("Source/executable file not found.")
        except PermissionError:
            print("Permission denied.")
        except Exception as e:  
            print("An error occurred:", e)


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
        

def grid_uniform(parameters,target): #change remove dunes
    
    def dunes(x,y):
        x_peak = np.deg2rad(35)
        y_peak = np.deg2rad(180)
        return np.exp(-((x - x_peak)**2 + (y - y_peak)**2) / 0.05)
    
    def crater():
        nonlocal mydir,parameters,target
        impactor = Crater_Impactor()
        
        base = Crater(parameters,target,impactor)
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
 
        x_values, y_values = grid(parameters) 
                


    # Open file to write the grid data
        with open(mydir+"/base.txt", "w") as file:
            for x in x_values:
                for y in y_values:
                    height = 1#dunes(x,y)
                    file.write(f"{x:.12f},{y:.12f},0,{height:.12f}\n")  # Format to 6 decimal places for precision
    
        print("Grid data with sand dune height profile saved to base.txt")
        
        crater()
        
        with open(mydir+"/grav.txt", "w") as file: 
            for x in x_values:
                for y in y_values:
                    file.write("-1 0 0\n")  # Format to 6 decimal places for precision

        print("Gravity data saved to grav.txt")

    if parameters['Dimension'] == '1D':
        
        res=int(parameters["Resolution"])
        base=np.zeros((res,3))
        base[:,2]=1
        offset=float(parameters["offset"])
        dx=(math.pi-2*offset)/res
        for i in range(res):
            base[i,0] = offset + dx * (i+ 0.5)
        np.savetxt(mydir+"/base.txt",base,delimiter=",")

        grav=np.zeros((res,4))
        grav[:,0] = base[:,0]
        grav[:,1] = -1
        np.savetxt(mydir+"/grav.txt",base,delimiter=",")
        
def grid(parameters):
    
    if parameters['Dimension'] == '2D':
        
        nx, ny = int(parameters['X Resolution']), int(parameters['Y Resolution'])

    # Define the grid ranges
        offset =float(parameters["offset"])
        dx=(np.pi-2*offset)/(nx)
        dy=(2*np.pi)/(ny)
        x_values = np.linspace(offset+dx/2, np.pi-offset-dx/2, nx)
        y_values = np.linspace(0+dy/2, 2 * np.pi-dy/2, ny)
        
        return x_values, y_values