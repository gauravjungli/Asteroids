#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jun 16 12:53:34 2024

@author: g
"""
from IO import read_xlsx_to_input_field_dict, Output_File
from GUI_asteroid import GUI


#%%


if __name__ == "__main__":
    
    
    parameters={"run":0}
    inputfile = Output_File(parameters, "input" ,["parameters.xlsx"])
    data_dict=read_xlsx_to_input_field_dict(inputfile)
    gui=GUI(data_dict)
    gui.root.mainloop()
    print("End of GUI")
    

    
