#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jun 16 12:53:34 2024

@author: g
"""
from gaurav import read_xlsx_to_input_field_dict, Output_File, Exparameter
from GUI_asteroid import GUI
import os, time, subprocess

#%%


  

if __name__ == "__main__":
    
    
    parameters={"run":0}
    inputfile = Output_File(parameters, "input" ,["parameters.xlsx"])
    data_dict=read_xlsx_to_input_field_dict(inputfile)
    gui=GUI(data_dict)
    gui.root.mainloop()
    start_time=time.time()
     
    end_time=time.time()
    elapsed_time=end_time-start_time

    
