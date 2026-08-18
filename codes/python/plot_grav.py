#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 22 16:22:43 2025

@author: g
"""
import time
from IO import Output_File
import re
import pandas as pd
import matplotlib.pyplot as plt
import glob
import os
import numpy as np
import numpy as np
from scipy.interpolate import make_interp_spline, CubicSpline
import multiprocessing
from collisions import G
from scipy.special import ellipk, ellipe,elliprf,elliprj
from gravity import Gravity
from debug import debug_plot
from matplotlib.ticker import ScalarFormatter
#%%

def Gravitycalc(target,file,rad):
    
    start =time.time()
    Res=target.res
    epsilon=1e-3 #This is a different epsilon
    density=target.dens
    
    try:
        w=np.loadtxt(file,dtype=float,delimiter=",")
    except:
        print("No file available for the fit")
        return   

    R=rad*np.sin(w[:,0])*w[:,1]
    Z=rad*np.cos(w[:,0])*w[:,1]
    fR=make_interp_spline(w[2:Res-2,0],R[2:Res-2]) 
    fZ=make_interp_spline(w[2:Res-2,0],Z[2:Res-2])

    res=10000

    theta=np.linspace(0,np.pi,res)
    r = fR(theta)
    z = fZ(theta)
    
    num_processes = 20# multiprocessing.cpu_count()#change for a single run
    pool = multiprocessing.Pool(processes=num_processes)
    
    arguments=[(R[i]+epsilon/2*rad*np.sin(w[i,0]),Z[i]+ epsilon/2*rad*np.cos(w[i,0]),r,z) for i in range(int(Res/2))]#changed it from Res/2
    
    grav = pool.starmap(Gravity, arguments)
    grav=np.array(grav)
    grav1=np.array([(grav[i,0],-grav[i,1]) for i in range(round(Res/2)-1,-1,-1)])#change_P uncomment if using equator symmetry
    grav=np.vstack((grav,grav1))
    
    R_grav = -G*density*(grav[:,0]*np.sin(w[:,0])+grav[:,1]*np.cos(w[:,0]))
    T_grav = -G*density*(grav[:,0]*np.cos(w[:,0])-grav[:,1]*np.sin(w[:,0]))
    
    #Only for the axisymmetric case
    theta = w[2:Res-2,0]
    base = w[2:Res-2,1]
    dbase = w[2:Res-2,3]

    metric = base**2 + dbase**2
    for i in range(2,Res-2):
        R_grav_temp = R_grav[i]
        R_grav[i] = (R_grav[i]*base[i-2]-dbase[i-2]*T_grav[i])/np.sqrt(metric[i-2])
        T_grav[i] = (T_grav[i]*base[i-2]+dbase[i-2]*R_grav_temp)/np.sqrt(metric[i-2])

 
    #plt.plot(w[:,0],t_grav)
    print("gravity updated")

    pool.close()
    pool.join()       
    end =time.time()
    print(f"Time taken in calculating gravity:{end-start}")
    return R_grav, T_grav

#%%
# 1. Find all Excel files in the current directory
def plot1():
    filename = "/home/g/Asteroids/output/New_1/run2/data"
    excel_files = glob.glob(filename+"/landslides_*")
    excel_files.sort(key=lambda f: [int(num) for num in re.findall(r'\d+', f)])
    
    # 2. Setup the plot style for a journal
    #plt.figure(figsize=(8, 6))
    plt.rcParams.update({'font.size' : 18})
    dia = pd.read_csv(filename+"/dia.txt",sep='\s+', engine='python',header=None)
    count=0
    N=33
    y_data=np.zeros((N,1))
    for i, file in enumerate(excel_files[0:N]):
        # Read the excel file
        # header=None if there are no headers, otherwise remove the argument
        df = pd.read_csv(file+'/field_0.csv',header=None)
        
        # Column 1 (index 0) is X, Column 8 (index 7) is Y
    
        y_data[i] = df.iloc[500, 7]*(dia.iloc[i,1]/dia.iloc[0,1])
    
        # Plot data
    plt.plot(dia.iloc[0:N-1,0], y_data[1:N],marker='o', markersize=4, linewidth=2)
    
    plt.xlabel("X-Axis Label", fontsize=18)
    plt.ylabel("Y-Axis Label (log)", fontsize=18)
    plt.grid(True)
    plt.tight_layout()
    plt.xticks(fontsize=18)
    plt.yticks(fontsize=18)

#%%
def plot2():
    # 1. Find all Excel files in the current directory
    filename = "/home/g/Asteroids/output/New_1/run1/data"
    excel_files = glob.glob(filename+"/landslides_*")
    excel_files.sort(key=lambda f: [int(num) for num in re.findall(r'\d+', f)])
    
    # 2. Setup the plot style for a journal
    #plt.figure(figsize=(8, 6))
    plt.rcParams.update({'font.size' : 18})
    dia = pd.read_csv(filename+"/dia.txt",sep='\s+', engine='python',header=None)
    count=0
    N=18
    y_data=np.zeros((N,1))
    #for i, file in enumerate(excel_files[0:N]):
        # Read the excel file
        # header=None if there are no headers, otherwise remove the argument
    df = pd.read_csv(filename+'/field_400.csv',header=None)
    
    #plt.clf()
    y_data = df.iloc[:, 7]
    x_data = df.iloc[:,0]
    
    plt.plot(x_data, y_data, linewidth=2)
    plt.pause(0.5)
    # Column 1 (index 0) is X, Column 8 (index 7) is Y
    
    
    
    # Plot data
    
    #%%
    
if __name__ == "__main__":
    file ='/home/g/Asteroids/output/New_1/run1/data/dia.txt'
    w =np.loadtxt(file,dtype=float)
    dia =w[:,1]
    steps = [1,100,200,300,400]
    
    parameters =debug_plot()
    fig, ax = plt.subplots()
    for step in steps:
        file =f'/home/g/Asteroids/output/New_1/run1/data/field_{step}.csv'
        w =np.loadtxt(file,dtype=float,delimiter=",")
        R_grav, T_grav = Gravitycalc(parameters,file,dia[step-1])
        ax.plot(w[:,0],(T_grav),linewidth=2,label = f'{step}') 
    xticks = [0,np.pi/12, np.pi/6, np.pi/4, np.pi/3, np.pi/2, 2*np.pi/3, 3*np.pi/4, 5*np.pi/6, np.pi]
    labels = ['0',r'$\dfrac{\pi}{12}$', r'$\dfrac{\pi}{6}$', r'$\dfrac{\pi}{4}$', r'$\dfrac{\pi}{3}$', r'$\dfrac{\pi}{2}$', r'$\dfrac{2\pi}{3}$', r'$\dfrac{3\pi}{4}$', r'$\dfrac{5\pi}{6}$', r'$\pi$']
    plt.xticks(xticks,labels)
    plt.xlim([0,np.pi/2])
    formatter = ScalarFormatter(useMathText=True)
    formatter.set_scientific(True)
    formatter.set_powerlimits((-2, 2))
    ax.yaxis.set_major_formatter(formatter)
    plt.grid(True)
    plt.tight_layout()