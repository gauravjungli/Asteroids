#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jul  6 01:33:23 2024

@author: g

This is made for debugging the codes. These are for the 1D cases only.
debug_main: This is used  to create and save files required for debugging the  main.py codes
debug_plot: This us used for creating the parameter files required for making plots from already available simulations
debug_cpp: For debugging the cpp files. It initializes the simulation as well as saves all the files like gravity and basal topography.
debug_post_process: For doing post processing to save the mean and std shapes
"""
from IO import Exparameter, Output_File_old, read_xlsx_to_input_field_dict, Parameter
from Initialize import Initialize_simulations, Initialize
import sys
from Target import Target
import os
import math
import numpy as np

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
import pdb

def Parameters(data_dict,parameters):
    for name in data_dict:
       # parameters[name]="Yes"
        inputs=data_dict[name]
        for Input in inputs:
            parameters[Input.Name]=Input.Value


class Impactor:
    
    def __init__(self,explicit=True):
        if explicit:
            self.phi   = math.acos(1 - 2 * np.random.random()) / 2
            self.vel   = 5.5e+3#scipy.stats.maxwell.ppf(random.random(), scale=3.232)*1000
            self.d     = 5 
            self.theta = 2 * math.pi * np.random.random()
            self.Phi = 2 * math.pi * np.random.random()
            self.Theta   = math.acos(1 - 2 * np.random.random())
        
            self.dens  = 1500

        self.M         = (math.pi / 6) * self.dens * self.d**3
        self.explicit  = explicit
        
        
# To debug the main script        
def debug_main():   
   # pdb.set_trace()
    parameters={"run":0}
    inputfile = Output_File_old(parameters, "input" ,["parameters.xlsx"])
    data_dict=read_xlsx_to_input_field_dict(inputfile)
    Parameters(data_dict,parameters)
    Initialize_simulations(parameters)
    return parameters

#call this command from the console
#  %debugfile /home/g/Asteroids/codes/python/main.py --wdir --args "'Output folder' 'Debug' 'run' '1'"

#for plotting the data
def debug_plot():
 #   pdb.set_trace()
    parameters ={}
    parameters['run'] = 1
    parameters['Output folder'] = "RF25"
    Parameter(parameters,'output')
    return parameters
    
# for doing post processing to save the mean and std shapes
def debug_post_process():  
    parameters ={}
    parameters["run"] = 0
    inputfile = Output_File_old(parameters, "input" ,["parameters.xlsx"])
    data_dict=read_xlsx_to_input_field_dict(inputfile)
    Parameters(data_dict,parameters)
    return parameters
   
 
#impactor =Impactor()
#for debugging cpp codes
def debug_cpp():
    
    #pdb.set_trace()
    parameters={"run":0}
    inputfile = Output_File_old(parameters, "input" ,["parameters.xlsx"])
    data_dict=read_xlsx_to_input_field_dict(inputfile)
    Parameters(data_dict,parameters)
    

   
   # parameters['Output folder'] = f'Spherical_{parameters["Rotation period"]}_{parameters["Static Friction angle"]}'
    Initialize_simulations(parameters)
    
    parameters['run'] = 1
    Parameter(parameters,"output")
    target = Target(parameters)
    Initialize(parameters, target)
    
    parameters['slides'] = 1
   # parameters['omega'] =0
    
    parameters['verbose_dir'] = Output_File_old(parameters,"output",['data',f"landslides_{parameters['slides']}"])
    

    Exparameter(parameters)  
    os.mkdir(parameters['verbose_dir'])
    return parameters
    
# data = Crater(parameters, target,impactor)
# mydir=Output_File (parameters,"output",["base.txt"])
# base=np.loadtxt(mydir,delimiter=",",dtype=float)
# base[:,2] = data[:,2]
# np.savetxt(mydir,base,delimiter=",")
# sys.argv = ['main.py', 'Output folder', 'dunes' , 'run','1']



if __name__ == "__main__":
    parameters = {}
    parameters = debug_plot()
    target = Target(parameters)
