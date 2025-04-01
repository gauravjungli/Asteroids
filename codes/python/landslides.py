#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 18 11:13:39 2023

@author: kumargaurav
"""

import math
import subprocess

import numpy as np
import os
import multiprocessing 
from collisions import G
import time
from Diffusion_spherical import  compute_height
from IO import Parameter, Exparameter, Output_File
from scipy.ndimage import gaussian_filter1d
from gravity import Gravitycalc
from Crater import Crater

    
#%%


"""
    It is used for the best fit circle. Currently we are only using the average of min and max.
"""

def Fit(parameters):
    
    
    Gamma=float(parameters["Gamma"])
    file=Output_File(parameters,"output",["base.txt"])
    try:
        w=np.loadtxt(file,dtype=float,delimiter=",")
    except:
        print("cannot import shape")
        return
    
    # w_new = gaussian_filter1d(w[:,1], sigma=20)

    # vol =  np.trapz(np.sin(w[:,0])*np.power((1+Gamma*w[:,1]),3),w[:,0])
    
    # vol_new =  np.trapz(np.sin(w[:,0])*np.power((1+Gamma*w_new),3),w[:,0])
    
    # w[:,1] = (vol/vol_new)**(1/3)*w_new
        
    r = 1+Gamma*(min(w[:,1])+max(w[:,1]))/2
    w[:,1] = ((1+Gamma*w[:,1])-r)/(Gamma*r)
    if parameters['Shallowness'] == 'Variable':
        Gamma_new = Gamma*np.average(np.abs(w[:,1]))
        if Gamma_new>0:
            w[:,1] = Gamma/Gamma_new*w[:,1]
        else:
            print("The Gamma value is zero and hence setting the basal topography to zero")
            w[:,1] = 0
        parameters["Gamma"] = Gamma_new
    parameters["current diameter"] = float(parameters["current diameter"])*r
    parameters["jinertia"] = float(parameters["jinertia"])/r**5
    parameters["jinertia1"] = float(parameters["jinertia1"])/r**5
    Exparameter(parameters)
    np.savetxt(file,w,delimiter=",") 
    print ("The best fit value of r is ", r)
  



  
 #%%  
"""
    This calculates the failure height of the landslide. First check whether landslide model is included in the 
    simulation or not. If yes then check whether it is initiated by an impact or it is a rotational failure. 
    For impacts, update epsilon for energy dependent simulations and select a failure profile from Gaussian and
    uniform. Update these in the base array and save it in a text file that will be later used by the C++ module.
""" 

def Height(parameters,target,impactor=None):
    
    if  not target.landslide:
        return
    
    mydir=Output_File (parameters,"output",["base.txt"])
    Gamma=float(parameters["Gamma"])
    epsilon=float(parameters["epsilon"])
    base=np.loadtxt(mydir,delimiter=",",dtype=float)
    min_epsilon = float(parameters['Minimum epsilon'])
    max_epsilon = float(parameters['Maximum epsilon'])
    
     
    if parameters['Dimension'] == '1D':
        
        res=int(parameters["Resolution"])
        height =np.ones(res)
        
        if impactor:
            
            H = compute_height(target=target,impactor=impactor)
            epsilon = np.clip(H,min_epsilon,max_epsilon)
            
            print(f'the destablization height is {epsilon} and the impactor diameter is {impactor.d}\
                  and impactor velocity is {impactor.vel}' )
        
        else:
       
            print("No impactor found to inititate landslide. Probably its a rotational failure")
            epsilon=10*min_epsilon
            
        base[:,1]=(base[:,1]*Gamma-epsilon*height)/np.abs(Gamma)
        base[:,2]=height[:]        
    
    else:
        base = Crater(parameters,target,impactor) #Fourth column contains the flag to check whether the point lies inside the crater
#change
        H=compute_height(target=target,impactor=impactor,grid=base)
        epsilon = (np.min(H) + np.max(H))/2
        epsilon = np.clip(epsilon,min_epsilon,max_epsilon)
        base[:,2] = (base[:,2]*Gamma-H[:])/np.abs(Gamma)
        base[:,3] = H[:]/epsilon  
    parameters['epsilon'] = epsilon   
    print(f"The epsilon: {parameters['epsilon']} and Gamma: {parameters['Gamma']}")

    

    np.savetxt(mydir,base,delimiter=",")



#%%
"""
This is the main function which calls the Landslide module:
    1. Update YORP if sochastic YORP is simulated. 
    2. Nondimensionalize all the relevant values such as omega, and inertia.
    3. Update slide count, impact time and dump directory of the landslide module.
    4. Export updated parameters
    5. Check whether the executable exist, the directories exist and create relevant directories
    6. Run the landslide executables
    7. Import the update values in the parameters
    8. Create the best fit circle and update the parameters dict
    9. Update omega,diameter and inertia of the target asteroid.
    10. Update gravity if necessary.
    
"""

def Landslides(target,parameters,impacttime,myomega):
    
    
    if  not target.landslide:
        return
    if float(parameters['epsilon'])<float(parameters['Minimum epsilon']):
        print(f"Value of epsilon {float(parameters['epsilon'])} too small to simulate landslide")
        return

    parameters["omega"] = target.omega[2] /(G * (4/3) * math.pi * target.dens)**0.5
    parameters['Seismic_shaking_time'] = target.t_lan/(target.d/2/target.grav)**0.5
    parameters["slides"] = int(parameters["slides"])+1
    parameters["current diameter"] = target.d
    parameters["jinertia"]=target.jinertia[2]/(target.d/2)**5/target.dens
    parameters["jinertia1"]=target.jinertia[0]/(target.d/2)**5/target.dens
    parameters["time"]=impacttime
    
    if bool(parameters['verbose']):
        parameters['verbose_dir']=Output_File(parameters,"output",['data',f"landslides_{parameters['slides']}"])
       
    try:
        Exparameter(parameters)
    except IOError:
        print("Failed in exporting the data")
    
   # exitcode = subprocess.run(["/home/g/Asteroids/build/asteroid"], cwd=".", stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL).returncode


    # 1. Ensure Executable Path & Permissions
    executable_path=Output_File(parameters,"output",[parameters["executable"]])
    if not os.path.isfile(executable_path) or not os.access(executable_path, os.X_OK):
        raise FileNotFoundError(f"Executable not found or not executable: {executable_path}")
    
    # 2. Adjust Working Directory if Necessary
    working_directory = Output_File(parameters,"output")  # Update if the executable is elsewhere
    if not os.path.isdir(working_directory):
        raise NotADirectoryError(f"Working directory not found: {working_directory}")
    
    try:
        # 3. Capture Output for Debugging (initially)
        os.mkdir(parameters['verbose_dir'])
        result = subprocess.run([executable_path], cwd=working_directory, capture_output=True, text=True) 
    
        if result.returncode != 0:
            print("Aborting due to error in running the executable landslide")
            raise subprocess.CalledProcessError(result.returncode, executable_path, output=result.stderr)
            
        else:
            print("Ran successfully slide", parameters["slides"])
    
        # If successful, you can later redirect output to DEVNULL 
        # result = subprocess.run([executable_path], cwd=working_directory, 
        #                        stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        
    except subprocess.CalledProcessError as e:
        print(f"Error running command: {e.cmd}")
        print(f"Return code: {e.returncode}")
        print(f"Error output:\n{e.stderr}")  
    

    Parameter(parameters,"output")
    
    if parameters['Dimension'] == '1D': #expand later for 2D
        Fit(parameters)
        Parameter(parameters,"output")
    
        target.omega[2] = float(parameters["omega"]) * (G * (4/3) * math.pi * target.dens )**0.5
        target.d = float( parameters["current diameter"])
    
        r = target.d / 2
        target.jinertia[2] = float(parameters["jinertia"]) * r**5 * target.dens
        target.jinertia[0] = target.jinertia[1] = float(parameters["jinertia1"]) * r**5 * target.dens
    
        if int(parameters["slides"])%10==0 or float(parameters["omega"])>0.9:
            print("Updating gravity")
            target.rgrav, target.tgrav = Gravitycalc(parameters)
    
    
    myomega.append([impacttime,target.omega[2]])
    print("Omega after the Landslides", target.omega[2])





    
