#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 17 12:24:53 2024

@author: g
"""

import numpy as np
import math
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick
from pathlib import Path
import glob
import os
import re
from gaurav import Parameter
import seaborn as sns


def show_shape(parameters):
    plt.close()
    fig,ax = plt.subplots()
    
    for count in range(0,100):
        
        file1=parameters['Data folder']
        
        file=os.path.join(file1,f'field_{(count+1)}.csv')
        if os.path.exists(file):
            
            w=np.loadtxt(file,delimiter=",",dtype=float)
            print(file)
           # if count!=slides-1 and count!=0:
            #    continue
            #plt.clf()
            x=np.sin(w[:,0])*(1+(float(parameters["epsilon"])*w[:,1]+float(parameters["Gamma"])*w[:,2]))
            y=np.cos(w[:,0])*(1+(float(parameters["epsilon"])*w[:,1]+float(parameters["Gamma"])*w[:,2]))
            ax.plot(x,y,'-r',linewidth=4)
            x=-np.sin(w[:,0])*(1+(float(parameters["epsilon"])*w[:,1]+float(parameters["Gamma"])*w[:,2]))
            ax.plot(x,y,'-r',linewidth=4)
            ax.set_aspect('equal')
            title=f'lanslide number={count+1}'
            ax.set_title(title)
            fig.canvas.draw()
            fig.canvas.flush_events()
            plt.show(block=False)
            fig.canvas.draw_idle()
            plt.pause(0.1)
            #plt.savefig(file1+"/img_"+str(count+1)+".svg",dpi=300,bbox_inches="tight")
        else:
            break
    #plt.pause(10)
    #plt.close(fig)
    

def show_omega(parameters):
    plt.close()
    fig,ax = plt.subplots()
    file2='/home/g/Asteroids/output'
    file1=f'Omega_{"C" if parameters["Collision"]=="Yes" else ""}{"L" if parameters["Landslide"]=="Yes" else ""}{"Y" if parameters["YORP"]=="Yes" else ""}.txt'
    file=os.path.join(file2,parameters['Output folder'],f'run{parameters["run"]}',file1)
    if os.path.exists(file):
      
        w=np.loadtxt(file,dtype=float)
        x=w[:,0]
        y=2*np.pi/w[:,1]/3600
        ax.plot(x,y,'-r',linewidth=2)
        title="Evolution of time period"
        ax.set_ylabel("Time period (hrs)")
        ax.set_xlabel("Simulation time")
        ax.set_title(title)
        fig.canvas.draw()
        fig.canvas.flush_events()
        plt.show(block=False)
        fig.canvas.draw_idle()
        
    plt.pause(0.1)
    #plt.close(fig)
     
     