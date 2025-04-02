#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jul  6 01:33:23 2024

@author: g
"""
from plots import post_process, show_omega
from IO import Exparameter, Output_File, read_xlsx_to_input_field_dict, Parameter
from Initialize import Initialize_simulations, Initialize
import sys
from Target import Target
import os
import math
import numpy as np
from landslides import Height
from Crater import Crater
#from main import main
def Parameters(data_dict,parameters):
    for name in data_dict:
        parameters[name]="Yes"
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
        
    
parameters={"run":0}
inputfile = Output_File(parameters, "input" ,["parameters.xlsx"])
data_dict=read_xlsx_to_input_field_dict(inputfile)
Parameters(data_dict,parameters)
Initialize_simulations(parameters)




# parameters['run'] = 1
# impactor =Impactor()
# Parameter(parameters,'output')
# target = Target(parameters)
# Initialize(parameters, target)
# if bool(parameters['verbose']):
#     parameters['slides'] = 1
#     parameters['Seismic_shaking_time'] = 5
#     parameters['verbose_dir']=Output_File(parameters,"output",['data',f"landslides_{parameters['slides']}"])
# Exparameter(parameters)  
# os.mkdir(parameters['verbose_dir'])
# data = Crater(parameters, target,impactor)
# mydir=Output_File (parameters,"output",["base.txt"])
# base=np.loadtxt(mydir,delimiter=",",dtype=float)
# base[:,2] = data[:,2]
# np.savetxt(mydir,base,delimiter=",")
# sys.argv = ['main.py', 'Output folder', 'dunes' , 'run','1']

#call this command from the console
#  %debugfile /home/g/Asteroids/codes/python/main.py --wdir --args "'Output folder' 'debug' 'run' '1'"


