#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Mar 29 16:43:12 2025

@author: g
"""
from openpyxl import Workbook, load_workbook
from openpyxl.styles import Alignment, Font, PatternFill, Border, Side
import re
import os
import copy
import numpy as np

#%%
""" The data structure is used for reading inputs from the xslx file """

class InputField:
    def __init__(self, Name, Value=None, Type="str", Options=None, Help=None):
        self.Name = Name
        self.Value = Value
        self.Type = Type
        self.Options = Options  # Ensure options is always a list
        self.Help= Help

    def __str__(self):
        return f"InputField(name='{self.name}', value='{self.value}', type='{self.type}', options={self.options}, help_text='{self.help_text}')"



#%%   
 
""" Writes InputField data to an XLSX workbook, each sheet named after a dictionary key, 
    with custom column widths.
"""
def write_input_fields_to_xlsx(data_dict, filename):
    

    try:
       wb = load_workbook(filename)
    except FileNotFoundError:
       wb = Workbook()  # Create new workbook if file doesn't exist

    
    column_widths = {
        "Options": 30,
        "Help": 60,
        "Type": 15
    }

    for sheet_name, input_fields in data_dict.items():

        if sheet_name in wb.sheetnames:  
            ws = wb[sheet_name]
            ws.delete_rows(2, ws.max_row) 
        else:
            ws = wb.create_sheet(sheet_name)

        # Header Styling
        header_font = Font(bold=True)
        header_fill = PatternFill(start_color="00FF00", end_color="00FF00", fill_type="solid")
        thin_border = Border(left=Side(style='thin'), right=Side(style='thin'), top=Side(style='thin'), bottom=Side(style='thin'))

        # Write header and set column widths
        header_attributes = ["Name", "Type", "Value", "Options","Help"] 
        for col_num, attr in enumerate(header_attributes, 1):  
            cell = ws.cell(row=1, column=col_num, value=attr.capitalize())
            cell.font = header_font
            cell.fill = header_fill
            cell.border = thin_border
            cell.alignment = Alignment(horizontal='center', vertical='center')

            # Set column width based on header attribute
            column_letter = ws.cell(row=1, column=col_num).column_letter
            ws.column_dimensions[column_letter].width = column_widths.get(attr.capitalize(), 20)  # Default 15 if not in dict

        # Write data
        for row_num, field in enumerate(input_fields, 2): 
            for col_num, attr in enumerate(header_attributes, 1):
                value = getattr(field, attr, "")
                cell_value = ", ".join(value) if isinstance(value, list) else value  # Handle lists
                cell = ws.cell(row=row_num, column=col_num, value=cell_value)
                cell.alignment = Alignment(wrapText=True, horizontal='center', vertical='center')
                cell.border = thin_border
                lines_needed = 5
                ws.row_dimensions[row_num].height = lines_needed * 12.75  # 12.75 is a rough approximation 
        

    # Remove the default sheet created by openpyxl
    #del wb["Sheet"]
    
    # Save workbook
    wb.save(filename)
    
#%% 
"""
    Reads an XLSX file with multiple sheets and returns a dictionary
    where keys are sheet names and values are lists of InputField objects.
"""

def read_xlsx_to_input_field_dict(filename):

    wb = load_workbook(filename, data_only=True)
    #wb.data_only=True
    data_dict = {}

    # Iterate over each sheet in the workbook
    for sheet_name in wb.sheetnames:
        ws = wb[sheet_name]  # Get the worksheet object
        input_fields = []  # List to store InputFields for this sheet
        header = [cell.value for cell in ws[1]]
        # Iterate over data rows, starting from the second row
        for row in ws.iter_rows(min_row=2):
            values = [cell.value for cell in row]
            input_field = InputField(
                Name=values[header.index("Name")],    Type=values[header.index("Type")], 
                        Value=values[header.index("Value")],
                        Options=values[header.index("Options")].split(",") if values[header.index("Options")] else [],
                        Help=values[header.index("Help")]
            )
            # Strip leading/trailing spaces from options
            input_field.Options = [option.strip() for option in input_field.Options]

            input_fields.append(input_field)

        # Add the list of InputFields to the dictionary using the sheet name as the key
        data_dict[sheet_name] = input_fields

    return data_dict


#%%   
"""
    It is the function which creates the parameters dict (which is the backbone of the data structure ).
"""

def Parameter(parameters,filetype):

    inputfile = Output_File(parameters, filetype ,["parameters"])
    with open(inputfile, "r") as file:
        for line in file:
            line = line.strip()
            if '--' in line:
                continue
            if line:
                values = re.split(r'\t+',line)
                key = values[0].strip()
                value = values[1].strip()
                parameters[key] = value
    return parameters
              
#%%  
"""
Writes the parameter dict to the parameters file. 
"""

def Exparameter(parameters, filetype="output"):
       
    mydir=Output_File (parameters,filetype,["parameters"])
    with open(mydir,"w") as f:
        for key in parameters.keys():
            f.writelines(["-"*100,"\n"])
            f.writelines(f'{key.ljust(40)}\t{parameters[key]}\n') 
        f.writelines(["-"*100])

#%%
""" Returns the folder location for writng and reading"""

def Output_File (parameters,filetype="",filenames=[]):
    
    directory = os.path.dirname(os.path.dirname(os.getcwd()))
    filenames1=copy.deepcopy(filenames)
    if filetype=="output":
        output_folder=parameters['Output folder']
        if int(parameters['run'])>0:
            filenames1.insert(0,"run"+str(parameters["run"]))
        filenames1.insert(0,output_folder)
        
    return os.path.join(directory,filetype,*filenames1)

#%% 
""" Cumulative distribution of the astroid population."""

def Cumdistr(parameters):
    
    pathcum = Output_File(parameters,"input",["cumpopulation", f"{parameters['Cummulative distribution']}.txt"])
    cumdistr = np.loadtxt(pathcum)
    cumdistr[:,0]=cumdistr[:,0]*1000
    return cumdistr

#%% 
""" Exports value of omega to the file. """

def ExportOmega(myomega,parameters):
    filename=f'Omega_{"C" if parameters["Collision"]=="Yes" else ""}{"L" if parameters["Landslide"]=="Yes" else ""}{"Y" if parameters["YORP"]=="Yes" else ""}.txt'
    mydir=Output_File(parameters,"output",[filename])
    
    resultExport = ""
    with open(mydir, "w") as file:
        resultExport = file.write("\n".join(["\t".join(map(str, omega)) for omega in myomega]))
    if resultExport == -1:
        print("Failed in exporting the data")
