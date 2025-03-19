#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 18 11:13:39 2023

@author: kumargaurav
"""

import math
import subprocess
from scipy.interpolate import make_interp_spline, CubicSpline
import numpy as np
from circle_fit import taubinSVD
import re
import os
from scipy.special import ellipk, ellipe,elliprf,elliprj, jv, jvp,lpn, hyp2f1
import multiprocessing 
from collisions import G,getdiaf,qstarf,probi,astnum, wobblecalcf
import time
import shutil
import copy
from openpyxl import Workbook, load_workbook
from openpyxl.styles import Alignment, Font, PatternFill, Border, Side
from scipy.optimize import root_scalar
from Diffusion_spherical import frequency, find_roots, find_roots_parallel, parallel_root_computation, compute_energy, compute_height
import numba
from scipy.ndimage import gaussian_filter1d


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
""" Represents a target object with physical properties and dynamics. It is the main data struture that holds 
    all the information about the target asteroid. The attributes of these data structure are as follows:
        d: Diameter
        atype: type of the asteroid (C-type or S-type)
        delta: Friction angle
        landslide, YORP, collision: Flags to include these processes in the simulation
        dens: Density 
        M: mass of the asteoid
        jinertia: The inertia tensor in principal coordinates
        omega: the angular velocity of the asteroid
        K: Thermal inertia
        Kvg: Holsapple parameter
        grav: Gravity field due to a sphere
        obliq: Is the obliquity paremeter used in the YORP calculation
        dstarave: Is the diameter for the catastrophic disruption
        coeff_f,coeff_g: Are the coefficients used in the YORP calculation for stochasticity
        sma: Semi major axis
        f_spline, g_spline: Used for YORP calcution. These are the spline fits for the YORP data
        k_s: is the diffusivity constant
        efficiency: is the seismic efficiency
        f: is the frequency
        Q: is the quality factor of the seismic waves
        N: Number of grid points at which the seismic energy is calculated
        roots: is the root of the bessel equations 
        theta: is the array of the latitudes at which the seismic energy is calculated
        energy: is the array of seismic energy when the impact energy is 1 Joule.
"""

class Target:
    


    def __init__(self, parameters):

        self.d=float(parameters["Diameter"])
        self.atype=parameters["atype"]
        self.delta=float(parameters['Friction angle'])
        self.landslide = True if parameters["Landslide"].lower()=='yes' else False
        self.YORP = True if parameters["YORP"].lower()=='yes' else False
        self.collision = True if parameters["Collision"].lower()=='yes' else False
        self.dens = float(parameters["Density"])
        if self.atype == "S-Type":
            self.mu = 0.55
            self.Y0, self.d0strength   =  1.44e7, 0.1
            self.nsize = 3  # strength decreases with size as 1/nsize
            self.qconst1, self.qconst2 =  1e3, 1e6
            self.k1, self.k2 = 0.06, 1
        else:
            # otherwise - C-Type
            self.mu = 0.41
            self.Y0, self.d0strength   = 1e5,  0.1
            self.nsize = 3  
            self.qconst1, self.qconst2 = 2e3, 4e5 
            self.k1, self.k2 = 0.15, 1
            
        self.M= (math.pi / 6) * self.dens * self.d**3
        self.jinertia= [2/5 * self.M * (self.d/2)**2]*3
        
        self.omega=[0,0,2*np.pi/(float(parameters['Rotation period'])*3600)]
        self.kvg=0.3
        self.K=float(parameters["K"])
        self.grav=G*self.M/(self.d/2)**2
        self.obliq=float(parameters["Obliquity"])
        velave = float(parameters["Impactor velocity"])
        self.dstarave=qstarf(self, math.pi / 4, velave)[2]
        # Set the seed for NumPy's random number generator
        np.random.seed(int(time.time()/float(parameters['run'])))
        self.coeff_f,self.coeff_g=shape_gen(self.K)
        self.sma= float(parameters['Semi major axis'])
        self.f_spline, self.g_spline = read_f_g_spline(parameters)
        self.wave_speed = float(parameters["P wave speed"])*(self.grav/1.6*1500/self.dens)**(1/4) 
        self.k_s = 1/3*100*self.wave_speed
        self.efficiency = float(parameters["Seismic efficiency"])
        self.f = float(parameters["Frequency"])
        self.Q = float(parameters["Q"])
        self.N = 50
        self.roots = parallel_root_computation(300,50,self.d)
        self.theta = np.linspace(np.pi/10, np.pi,self.N)
        self.energy, self.t_lan = compute_energy(self)
        self.cohesion_cons = float(parameters["Cohesion constant"])
        self.cohesion_linear = float(parameters["Cohesion linear"])
        
        self.rgrav=None 
        self.tgrav = None
        
    
    """ Not currently in use. Using new parallel version from difusion_spherical.py"""
    def Roots(self):
        n=300
        m=50
        start =time.time()
        roots = np.zeros((n+1,m))

        #Find roots
        for i in range(0,n+1):
            
            roots[i,:] = find_roots_parallel(i,num_roots=m)/(self.d/2)
        end =time.time()
        print(f"Time taken in finding roots:{end-start}")
        return roots
      
""" 
    Is the class of impactors. Its attributes are as follows:
        explict: IS the flag to check whether an impact is big enough to cause landslides
        impacttime: is the time of impact of each impact
        M: is the masss of each impactor
        dia: is the function that generates the random number which represents the number of asteroid greater
            than a specific dia. This dia is returned by the function getdiaf
        d: is the diameter of the impactor
        theta, Theta, Phi: are the angles of the impact.
 
"""

class Impactor:
    
    def __init__(self,tmaxby,velave,low,high,cumdistr,explicit=True):
        if explicit:
            self.phi   = math.acos(1 - 2 * np.random.random()) / 2
            self.vel   = velave#scipy.stats.maxwell.ppf(random.random(), scale=3.232)*1000
            self.d     = self.dia(low,high,cumdistr)   
            self.theta = 2 * math.pi * np.random.random()
            self.Theta = 2 * math.pi * np.random.random()
            self.Phi   = math.acos(1 - 2 * np.random.random())
        else:
            self.d     = math.exp((math.log(low) + math.log(high)) / 2)
            self.phi   = math.pi / 4
            self.vel   = velave
            self.theta = math.pi
            self.Theta = 0
            self.Phi   = math.pi / 2
        if np.random.randint(1, 4) == 1:
            self.dens  = 2500
        else:
            self.dens  = 1500
        self.impacttime= np.random.uniform(0, tmaxby)
        self.M         = (math.pi / 6) * self.dens * self.d**3
        self.explicit  = explicit
        
    def dia(self,low,high,cumdistr):
        d=np.random.randint(low=low, high=high)
        return getdiaf(d,cumdistr)

    
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
    
    #uncomment for the special fit
    #res=int(parameters["res"])
    #x=np.sin(w[:,0])*(1+Gamma*w[:,1])
    #y=np.cos(w[:,0])*(1+Gamma*w[:,1])
     
    #point = []
    #for i in range(res):
     #   point.append([x[i],y[i]])
     #   point.append([-x[i],y[i]])
    #xc, yc, r, sigma = taubinSVD(point)
    w[:,1] = gaussian_filter1d(w[:,1], sigma=20)
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
    return r


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
    parameters['current diameter']       = target.d 
    parameters['Gamma']     =  0
    parameters['epsilon']   =  0   
    mydir=Output_File (parameters,"output")

    
    Exparameter(parameters)
    file = Output_File(parameters, "output", ["output.yorp"])
    with open(file, "w") as output_file:  # Overwrites existing files
        output_file.write(
            f"{0 :12.6e} {target.omega[2]:12.8e} {target.obliq:12.8e}\n")
   
    res=int(parameters["Resolution"])
    base=np.zeros((res,3))
    base[:,2]=1
    offset=float(parameters["offset"])
    dx=(math.pi-2*offset)/res
    for i in range(res):
        base[i,0] = offset + dx * (i+ 0.5)
    np.savetxt(mydir+"/base.txt",base,delimiter=",")
    
    if target.landslide:
        target.rgrav,target.tgrav =  Gravitycalc(parameters) 

        executable_file=Output_File(parameters,"build",[parameters["executable"]])
        
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
    res=int(parameters["Resolution"])
    offset=float(parameters["offset"])
    Gamma=float(parameters["Gamma"])
    epsilon=float(parameters["epsilon"])
    base=np.loadtxt(mydir,delimiter=",",dtype=float)
    min_epsilon = float(parameters['Minimum epsilon'])
    max_epsilon = float(parameters['Maximum epsilon'])
    height =np.zeros(res)
    def gaussian(x, mean, std_dev):
        return np.exp(-((x - mean)**2) / (2 * std_dev**2)) / (std_dev * np.sqrt(2 * np.pi))
    

    if impactor:
        
            
        if parameters["Failure depth"]=='Energy dependent':
            H=compute_height(target=target,impactor=impactor)
        if parameters["Shallowness"] == 'Variable':
            epsilon=np.clip(H,min_epsilon,max_epsilon)


        print(f'the destablization height is {epsilon} and the impactor diameter is {impactor.d} and impactor velocity is {impactor.vel}' )
        if parameters['Failure profile'].lower()=="gaussian":
            
            #for Gaussian profiles
            mean = impactor.Phi        # Mean of the distribution in degrees
            std_dev =   math.pi/10     # Standard deviation of the distribution in degrees
        
            gaussian_values = gaussian(base[:,0], mean, std_dev)
        
            dx=(math.pi-2*offset)/res
            area_under_curve = sum((1+Gamma*base[:,1])**2*gaussian_values*np.sin(base[:,0])*dx)
            scaling_factor = sum((1+Gamma*base[:,1])**2*np.sin(base[:,0])*dx)/ area_under_curve
            height = gaussian_values * scaling_factor
            
        elif parameters['Failure profile'].lower()=="uniform":#for uniform profiles
            height= np.ones(res)
        
    else:
       
        print("No impactor found to inititate landslide. Probably its a rotational failure")
        epsilon=2*min_epsilon
        height=np.ones(res)
    if parameters["Shallowness"] == 'Variable':
        parameters['epsilon'] = epsilon
        parameters['Gamma']  = np.abs(Gamma +epsilon)
        print(f"The epsilon: {parameters['epsilon']} and Gamma: {parameters['Gamma']}")

        base[:,1]=(base[:,1]*Gamma-epsilon*height)/np.abs(Gamma+epsilon)
        base[:,2]=height[:]
    else:
        base[:,1]=(base[:,1]*Gamma-epsilon*height)/np.abs(Gamma)
        base[:,2]=height[:]*np.min(H/epsilon,max_epsilon/epsilon)
    #base=np.column_stack((base,height))
    np.savetxt(mydir,base,delimiter=",")

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
"""
Calculates gravity for the axisymmetric body.
"""

def Gravitycalc(parameters):
    
    start =time.time()
    Res=int(parameters["Resolution"])
    epsilon=0.001 #This is a different epsilon
    Gamma=float(parameters["Gamma"])
    density=float(parameters["Density"])
    rad=float(parameters["current diameter"])/2
    
    file=Output_File(parameters,"output",["base.txt"])

    try:
        w=np.loadtxt(file,dtype=float,delimiter=",")
    except:
            print("No file available for the fit")
            return
        
    R=rad*np.sin(w[:,0])*(1+Gamma*(w[:,1]))
    Z=rad*np.cos(w[:,0])*(1+Gamma*(w[:,1]))
    
    fR=make_interp_spline(w[:,0],R)
    fZ=make_interp_spline(w[:,0],Z)

    res=10000

    theta=np.linspace(0,math.pi,res)
    r = fR(theta)
    z = fZ(theta)
    
    num_processes = multiprocessing.cpu_count()
    pool = multiprocessing.Pool(processes=num_processes)
    
    arguments=[(R[i]+epsilon/2*rad*np.sin(w[i,0]),Z[i]+ epsilon/2*rad*np.cos(w[i,0]),r,z) for i in range(int(Res/2))]
    
    grav = pool.starmap(Gravity, arguments)
    grav=np.array(grav)
    grav1=np.array([(grav[i,0],-grav[i,1]) for i in range(round(Res/2)-1,-1,-1)])
    grav=np.vstack((grav,grav1))
    
    R_grav = -G*density*(grav[:,0]*np.sin(w[:,0])+grav[:,1]*np.cos(w[:,0]))
    T_grav = -G*density*(grav[:,0]*np.cos(w[:,0])-grav[:,1]*np.sin(w[:,0]))
    r_grav = R_grav/(4/3*np.pi*density*rad*G)
    t_grav = T_grav/(4/3*np.pi*density*rad*G)
   # plt.plot(w[:,0],r_grav)
   # plt.plot(w[:,0],t_grav)
    print("gravity updated")
    grav=np.hstack((r_grav.reshape(-1,1),t_grav.reshape(-1,1)))
    file=Output_File(parameters,"output",["grav.txt"])
    np.savetxt(file,grav)
    
    pool.close()
    pool.join()       
    end =time.time()
    print(f"Time taken in calculating gravity:{end-start}")
    return R_grav, T_grav

def Gravity(R, Z,r,z): 
    
    r_grav=0
    z_grav=0
    for i in range(1,len(r)-1):
        a = r[i] # radius of disc being integrated
        zeta = Z - z[i] # vertical disctance of disc from point of evaluation
        delta = np.sqrt((a + R)**2 + (zeta)**2) # parameter for elliptic integrals
        k = 2 * np.sqrt(a * R) / delta  
        m = 2 * np.sqrt(a * R) / (a + R)
        if (R < a):
            eps = 1
        elif (R > a):
            eps = 0
        else:
            eps = 0.5
        if (k>=1 or m>=1):
            print("k = ",k," m = ",m,Z,z[i],a,R)
            m=min(m,1-1e-6)
        ks = ellipk(k**2)
        es = ellipe(k**2)
        pi=elliprf(0,1-k**2,1)+1/3*m**2*elliprj(0,1-k**2,1,1-m**2)

        r_grav=r_grav + np.abs((z[i]-z[i-1]))*(2 * delta * ((1 - k**2 / 2) * ks - es) / R)
        z_grav=z_grav + np.abs((z[i]-z[i-1]))*(2 * np.pi * np.sign(zeta) * eps + 2 * zeta * ((R - a)/(R + a) * pi - ks) / delta)
    return r_grav, z_grav    
    
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
    
    if parameters["stoc_yorp"].lower()=='yes':
        target.coeff_f, target.coeff_g = shape_gen(target.K)
    if  not target.landslide:
        return
    if float(parameters['epsilon'])<float(parameters['Minimum epsilon']):
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
            print("Aborting due to error in running the executable gaurav")
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
        
    min_epsilon=float(parameters['Minimum epsilon'])
    

    Parameter(parameters,"output")

    fit=Fit(parameters)
    Parameter(parameters,"output")
    
    target.omega[2] = float(parameters["omega"]) * (G * (4/3) * math.pi * target.dens )**0.5
    target.d = float( parameters["current diameter"])
    
    r = target.d / 2
    target.jinertia[2] = float(parameters["jinertia"]) * r**5 * target.dens
    target.jinertia[0] = target.jinertia[1] = float(parameters["jinertia1"]) * r**5 * target.dens
    
    if int(parameters["slides"])%10==0 or float(parameters["omega"])>0.85:
        print("Updating gravity")
        target.rgrav, target.tgrav = Gravitycalc(parameters)
    
    
    myomega.append([impacttime,target.omega[2]])
    print("Omega after the Landslides", target.omega[2])

#%% 

"""old version of YORP which uses orbit9, a fortran code, Please use the newer version of the YORP implemented in YORP.py """

def Yorp(target,parameters,impacttime,oldtime,myomega):
    
    wobblecalcf(target,impacttime,oldtime)  #which omega to use because it is being changed by the yorp
    while impacttime > oldtime + 10:
        subprocess.run(["make"], cwd="../OrbFit/tests/gaurav",stdout=subprocess.DEVNULL,
    stderr=subprocess.STDOUT)
        try:
            yark = []
            with open("../OrbFit/tests/gaurav/yarkovsky.in", "r") as file:
                yark = [list( line.strip().split()) for line in file]
        except FileNotFoundError:
            print("Failed in importing the data from yarkovsky.in")
            return
        yark[0][5] = target.obliq
        yark[0][6] = 2 * math.pi / (target.omega[2] * 3600)
        try:
            with open("../OrbFit/tests/gaurav/yarkovsky.in", "w") as file:
                for row in yark:
                    file.write("\t".join(map(str, row)) + "\n")
        except IOError:
            print("Failed in exporting the yarkovsky.in")
            return
        exitcode = subprocess.run(["./orbit9.x"], cwd="../OrbFit/tests/gaurav").returncode
        if exitcode != 0:
            print("Failed in running orbit9")
            return
        try:
            omegOrb = []
            with open("../OrbFit/tests/gaurav/clo0.yorp", "r") as file:
                omegOrb = [list( line.strip().split()) for line in file]
        except FileNotFoundError:
            print("Failed in importing the data from orbit9")
            return
        omegaLimit = (G*4/3* math.pi*target.dens)**0.5
        index = min(round((impacttime - oldtime) / 50), 2000)
        target.omega[2] = 2 * math.pi / (float(omegOrb[index][1]) * 3600)
        target.obliq = float(omegOrb[index][2])
        oldtime = oldtime + 10**5
        myomega.append([min(impacttime,oldtime), target.omega[2]])
        if target.omega[2] > 0.9*omegaLimit:
            print("Too fast spinning causing landslides")
            #parameters["uni_h"]=min(max((target.omega[2]-0.9*omegaLimit)/(omegaLimit)*(0.2/float(parameters["epsilon"])),1),10)
            #print(parameters["uni_h"])
            Height(parameters,target)
            Landslides(target,parameters,min(impacttime,oldtime),myomega)
        
        
    print("Omega after the yorp effect:", target.omega[2])


#%% 
""" Cumulative distribution of the astroid population."""

def Cumdistr(parameters):
    
    pathcum = Output_File(parameters,"input",["cumpopulation", f"{parameters['Cummulative distribution']}.txt"])
    cumdistr = np.loadtxt(pathcum)
    cumdistr[:,0]=cumdistr[:,0]*1000
    return cumdistr

#%% verified
""" This creates the collision history. Using Poisson's distribution for determining the number of impactors and deducing th
    the collision time """

def Istuff(parameters,target,tmaxby,cumdistr):
    prob = probi * tmaxby * (target.d/2) ** 2
    explicit_cutoff = float(parameters['Explicit cutoff'])
    numgtd = round(astnum(explicit_cutoff*target.dstarave,cumdistr)[0])
    dialittle = cumdistr[-1][0]
    velave = float(parameters['Impactor velocity'])
    energy_cons =  (math.pi* target.efficiency*velave**2*target.dens/12)*target.energy[-1]
    energy_min = target.cohesion_cons**2/(2*target.wave_speed**2*target.dens*np.tan(target.delta*math.pi/180)**2)
    dexplicit =  (energy_min/energy_cons)**(1/3)
    dexplicit = max(1.1*dialittle, dexplicit)
    binexplicit = astnum(dexplicit,cumdistr)[1]
    nexplicit = round(astnum(dexplicit,cumdistr)[0])
    
    nexpimpactors = round(prob * (nexplicit-numgtd))
    nexpimpactors = poisson_from_exponential(nexpimpactors)
    print(f'Number of expected impactor is {nexpimpactors}')
    density = 1500
    dimplicit = ((G**2*target.dens**3*target.d**5)/(9*target.efficiency*density*velave**2*target.f**2))**(1/3)*(
                np.exp(2*math.pi*target.f*target.d**2/(target.k_s*math.pi**2*target.Q)))
    binimplicit = astnum(dimplicit,cumdistr)[1]

    istuff=[]
    for j in range(nexpimpactors):
        istuff.append(Impactor(tmaxby=tmaxby,velave=velave,low=numgtd,high=nexplicit,cumdistr=cumdistr,explicit=True))


    for j in range(binexplicit+1, min(binimplicit,len(cumdistr)-1)):
        istuff.append(Impactor(tmaxby=tmaxby,velave=velave,low=cumdistr[j,0],high=cumdistr[j+1,0],cumdistr=cumdistr,explicit=False))

    istuff.sort(key=lambda x: x.impacttime)
    print(f"The minimum diameter for creating landslides is {dexplicit}")
    return istuff


def poisson_from_exponential(lambda_rate, max_time=1):
    """
    Simulate a Poisson-distributed random variable using exponential inter-arrival times.
    
    Parameters:
    - lambda_rate: The rate (λ) of the Poisson process.
    - max_time: The maximum time interval fGammaor the simulation (default is 1).
    
    Returns:
    - Number of events (Poisson-distributed) in the interval [0, max_time].
    """
    num_events = 0
    cumulative_time = 0

    # Generate exponential random variables until the cumulative time exceeds max_time
    while cumulative_time <= max_time:
        # Generate the time between the next event (Exponential random variable)
        inter_arrival_time = np.random.exponential(1 / lambda_rate)
        
        # Update cumulative time
        cumulative_time += inter_arrival_time
        
        # If the event happens within the time interval, increase event count
        if cumulative_time <= max_time:
            num_events += 1

    return num_events


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

def read_f_g_spline(parameters):
    """
    Reads data from files, creates spline interpolations for f and g functions.

    This function assumes the following file structure:

    - input/yorp_f.txt: Contains gamma (in degrees) and corresponding f values.
    - input/yorp_g.txt: Contains gamma (in degrees) and corresponding g values.

    Returns:
        tuple: Two CubicSpline objects representing the spline interpolations for f and g.
    """

    # --- Read data for the f function ---
    f = []
    mydir=Output_File(parameters,"input",["yorp_f.txt"])
    with open(mydir, 'r') as file:
        for line in file:
            x, y = map(float, line.split(","))
            f.append([x,y])
          
    f=np.array(f)

    # Create cubic spline interpolation for f
    f_spline = CubicSpline(f[:,0], f[:,1])

    g= []
    mydir=Output_File(parameters,"input",["yorp_g.txt"])
    with open(mydir, 'r') as file:
        for line in file:
            x, y = map(float, line.split(","))
            g.append([x,y])
    
    g=np.array(g)

    # Create cubic spline interpolation for g
    g_spline= CubicSpline(g[:,0], g[:,1])

    return f_spline, g_spline


#%%

def shape_gen(K):
    """
    generate random coefficients to determine the functions f,g. To
    this purpose, we use the statistics presented in Capek & Vokrouhlicky 2004

    Args:
        K (float): Input value used for determining probabilities.

    Returns:
        tuple: A tuple containing the generated values of coeff1 and coeff2.
    """

    K_t = 0.005
    max_g = 1.8 / 1.1
    min_g = 0.4 / 1.1
    std_g = abs(1.1 - 1.8 / 1.1) / 3.0
    max_f = 3.0 / 2.0
    min_f = -3.0 / 2.0
    std_f = 0.5**2.0

    if K <= K_t:
        # 80% probability to reach 0/180 (g, coeff2)
        # 40% probability to accelerate (f, coeff1)
        coeff2 = np.random.normal(1.0, std_g)
        coeff2 = np.clip(coeff2, min_g, max_g)  # Ensure coeff2 is within bounds

        if np.random.rand() > 0.8:
            coeff2 = -coeff2  # Switch to reaching 90 degrees

        # ! If the asymptotic state is 90 (i.e. coeff2 < 0), then we always decelerate
        # ! the rotation rate! See Capek & Vokrouhlicky 2004, Fig 7.
        if coeff2 < 0.0:
            # Always decelerate if asymptotic state is 90
            sgn = -1.0
        else:
            # 60/40 probability to decelerate/accelerate for 0/180
            sgn = 1.0 if np.random.rand() <= 0.6 else -1.0 

        coeff1 = sgn * (np.random.normal(1.0, std_f))
        coeff1 = np.clip(coeff1, min_f, max_f)

    else:  # K > K_t
        # 100% probability to reach 0/180 (g, coeff2)
        # 50% probability to accelerate (f, coeff1)

        sgn = 1.0 if np.random.rand() > 0.5 else -1.0  # Choose sign for coeff1
        coeff1 = sgn * (np.random.normal(1.0, std_f))
        coeff1 = np.clip(coeff1, min_f, max_f)
       
        coeff2 = np.random.normal(1.0, std_g)
        coeff2 = np.clip(coeff2, min_g, max_g)
    print(f"The value of the coefficients are {coeff1} and {coeff2}")
    return coeff1, coeff2

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