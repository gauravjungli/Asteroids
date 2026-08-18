#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jun 16 12:53:34 2024

@author: g

This is the main script that should be run to launch the GUI. 

# %%
It first creates a parameter dict which is initialized by reading the parameters.xslx file.

Then it starts a GUI
"""
import sys
import os
sys.path.append(os.path.abspath(os.path.dirname(__file__)))
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__),'script_1D')))
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__),'script_2D')))
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__),'Failure_profile')))
from IO import read_xlsx_to_input_field_dict, Output_File_old
from GUI_asteroid import GUI


#%%

if __name__ == "__main__":
    
    #create parameters dict
    parameters={"run":0}
    #name of the input file for creating the dict
    inputfile = Output_File_old(parameters, "input" ,["parameters.xlsx"])
    # read the input file 
    data_dict=read_xlsx_to_input_field_dict(inputfile)
    # create a GUI object using the dictionary
    gui=GUI(data_dict)
    # calls the mainloop of the GUI. GUI starts at this point.
    gui.root.mainloop()
    print("End of GUI")
    
