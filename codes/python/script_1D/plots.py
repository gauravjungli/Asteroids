#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 17 12:24:53 2024

@author: g
"""
"""
This script is used to show results while running simulations.
"""
import numpy as np

import matplotlib.pyplot as plt
import matplotlib.ticker as mtick
from matplotlib.figure import Figure
from pathlib import Path
import os
import pdb
from IO import Output_File_old, ExportOmega
import seaborn as sns
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import tkinter as tk
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import matplotlib.gridspec as gridspec
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.animation import FuncAnimation
from matplotlib.colors import LinearSegmentedColormap, Normalize
from matplotlib.cm import ScalarMappable
from scipy.ndimage import gaussian_filter1d
from Fit import Fit, Fit_radius

def show_shape(parameters):
    #pdb.set_trace()
    plt.close()
    fig,ax = plt.subplots()
    file1=parameters['Data folder']
    file2=os.path.join(file1,'dia.txt')
    epsilon=np.loadtxt(file2,dtype=float)[:,2]

    for count in range(0,10000,1):
        
        
        file=os.path.join(file1,f'field_{(count+1)}.csv')
        
        if os.path.exists(file):
           
            w=np.loadtxt(file,delimiter=",",dtype=float)

            ax.cla()#change
            x = np.sin(w[:,0])*w[:,1]
            y = np.cos(w[:,0])*w[:,1]
            ax.plot(x,y,'-r',linewidth=2)
            x =-np.sin(w[:,0])*w[:,1]
            ax.plot(x,y,'-r',linewidth=2)
            ax.set_aspect('equal')
            title=f'landslide number={count+1}'
            ax.set_title(title)
            fig.canvas.draw()
            #fig.canvas.flush_events()
            plt.show(block=False)
            #fig.canvas.draw_idle()
            plt.pause(0.1)
         #   plt.savefig()

        else:
            break
    #plt.pause(10)
    #plt.close(fig)
    

def show_omega(parameters):
    plt.close()
    fig,ax = plt.subplots()
    file2 = "output.yorp"
    file1 = f'Omega_{"C" if parameters["Collision"]=="Yes" else ""}{"L" if parameters["Landslide"]=="Yes" else ""}{"Y" if parameters["YORP"]=="Yes" else ""}.txt'
    file = Output_File_old(parameters,'output',[file1])
    
    if not os.path.exists(file):
        file = Output_File_old(parameters,'output',[file2])
   
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
        
    else:
        
        print("No omega file exists")
    
    plt.pause(0.1)
    #plt.close(fig)


     
def post_process(parameters):
    #pdb.set_trace()
    length=101
    if parameters["YORP"]=="Yes":
        file1="output.yorp"
    else:   
        file1=f'Omega_{"C" if parameters["Collision"]=="Yes" else ""}{"L" if parameters["Landslide"]=="Yes" else ""}{"Y" if parameters["YORP"]=="Yes" else ""}.txt'
    
    N=int(parameters['Number of simulations'])
    T=float(parameters['Simulation period'])*1e+6
    res=int(parameters['Resolution'])
 
    x=np.zeros((length,N+1))
    
    x[:,0]=np.linspace(0, T,num=length, endpoint=True)
    
    for i in range (1,N+1):
        file = Output_File_old(parameters=parameters,filetype='output',filenames=[f'run{i}',file1])     
        try:
            w=np.loadtxt(file,dtype=float)
        except FileNotFoundError:
            print(f"File {file} does not exist. Post processing can't be done")
            return
        x[:,i] = 2*np.pi/np.interp(x[:,0], w[:,0], w[:,1])/3600   
        
    myomega=[[x[i,0],np.mean(x[i,1:N+1]),np.std(x[i,1:N+1])] for i in range(length) ]
    ExportOmega(parameters = parameters, myomega=myomega)
    
    shape=np.ones((length+1,res+1,N))
        
    for i in range (0,N):
        
        file2 = Output_File_old(parameters=parameters,filetype='output',filenames=[f'run{i+1}','base.txt']) #change
        try:
            grid = np.loadtxt(file2,dtype=float,delimiter=",")
        except FileNotFoundError:
            print(f"File {file2} does not exist. Post processing can't be done")
            return
        
        shape[0,1:res+1,i]=grid[:,0]
        shape[1:length+1,0,i]=x[0:length,0]
        
        if parameters["Landslide"]=="Yes":
            
            file2 = Output_File_old(parameters=parameters,filetype='output',filenames=[f'run{i+1}','data','dia.txt'])
            try:
                slides=np.loadtxt(file2,dtype=float,ndmin = 2)
                epsilon = slides[:,2]
                dia =  slides[:,1]
            except FileNotFoundError:
                print(f"File {file2} does not exist. Post processing can't be done")
                return
            
            base=np.zeros((len(slides),res))
            
            for j in range(0,len(slides)):
                
                file2 = Output_File_old(parameters=parameters,filetype='output',
                                    filenames=[f'run{i+1}','data',f'field_{j+1}.csv'])
                try:
                    w = np.loadtxt(file2,delimiter=",",dtype=float)
                except FileNotFoundError:
                    print(f"File {file2} does not exist. Post processing can't be done")
                    return
                
                base[j,:] = Fit_radius(w,epsilon[j])*dia[j]/2

            for j in range(1,length+1):
                index=np.searchsorted(slides[:,0], x[j-1,0],side='right')
                if index==0:
                    shape[j,1:res+1,i] = float(parameters['Diameter'])/2
                else:
                    shape[j,1:res+1,i] = base[index-1,:]   
                  
    # First row of mean shape contains theta and first column contains time      
    mean_shape=np.power(np.mean(np.power(shape,3),axis=2),1/3)
    std_shape=np.std(shape,axis=2)
    std_shape[0,1:res+1]=mean_shape[0,1:res+1]
    std_shape[1:length+1,0]=mean_shape[1:length+1,0]
    
    file = Output_File_old(parameters=parameters,filetype='output',filenames=['mean_shape.csv']) 
    np.savetxt(file, mean_shape, delimiter=',', fmt='%f')
    file = Output_File_old(parameters=parameters,filetype='output',filenames=['std_shape.csv']) 
    np.savetxt(file, std_shape, delimiter=',', fmt='%f')
    
""" For old style 2D plots"""

def show_plots(parameters, root):
    file = Output_File_old(parameters=parameters, filetype='output', filenames=['mean_shape.txt']) 
    w = np.loadtxt(file, dtype=float, delimiter=",")
    file1=f'Omega_{"C" if parameters["Collision"]=="Yes" else ""}{"L" if parameters["Landslide"]=="Yes" else ""}{"Y" if parameters["YORP"]=="Yes" else ""}.txt'
    file = Output_File_old(parameters=parameters, filetype='output', filenames=[file1]) 
    omega = np.loadtxt(file, dtype=float)
    
    length, width = np.shape(w)

    fig = Figure(figsize=(16, 9), dpi=100)
    ax1 = fig.add_subplot(121)
    ax2 = fig.add_subplot(122)

    # Explicitly specify the height and width of each subplot
   # ax1.set_position([0.10, 0.2, 0.35, 0.7])  # [left, bottom, width, height]
   # ax2.set_position([0.52, 0.12, 0.42, 0.84])  # [left, bottom, width, height]


    canvas = FigureCanvasTkAgg(fig, master=root)
    canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=True)

    def update_plot(i):
        # First subplot
        
        ax2.cla()
        x2 = np.sin(w[0, 1:width]) * w[i, 1:width]
        y2 = np.cos(w[0, 1:width]) * w[i, 1:width]
        ax2.plot(x2, y2, '-r', linewidth=2)
        x2 = -np.sin(w[0, 1:width]) * w[i, 1:width]
        ax2.plot(x2, y2, '-r', linewidth=2)
        ax2.set_aspect('equal')
        ax2.set_title(f'Time = {w[i,0]}')
        
        ax1.cla()
        
        x1 =  omega[0:i-1, 0]
        y1 =  omega[0:i-1, 1]
        ax1.plot(x1, y1, '-b', linewidth=2)
        ax1.set_title('Rotation period')
        
        pos1 = ax1.get_position()  # Get position of the first subplot
        pos2 = ax2.get_position()  # Get position of the second subplot

# Set the position of the second subplot to match the first subplot
        ax1.set_position([pos1.x0, pos2.y0, pos1.width, pos2.height])



        canvas.draw()

    def animate():
        for i in range(1, length):
            update_plot(i)
            root.update_idletasks()
            root.after(100)

    root.after(0, animate)


# Function to animate the plot


def show_3D_plots(gui):
    
    parameters = gui.parameters
    root = gui.root
    plt.rcParams.update({'font.size' : 11})
    file = Output_File_old(parameters=parameters, filetype='output', filenames=['mean_shape.csv']) 
    w = np.loadtxt(file, dtype=float, delimiter=",")
    file1=f'Omega_{"C" if parameters["Collision"]=="Yes" else ""}{"L" if parameters["Landslide"]=="Yes" else ""}{"Y" if parameters["YORP"]=="Yes" else ""}.txt'
    file = Output_File_old(parameters=parameters, filetype='output', filenames=[file1]) 
    omega = np.loadtxt(file, dtype=float)
    
    length, width = np.shape(w)
    length -= 1       # First row of mean shape contains theta and first column contains time   
    width  -= 1
    gui.fig = Figure(figsize=(16, 9), dpi=300) 
    fig = gui.fig
    fig.patch.set_facecolor('black')
    ax2 = fig.add_axes([0.27, -0.3, 0.7, 1.6], projection='3d') 
    ax1 = fig.add_axes([0.10, 0.25, 0.3,0.6])
    cbar_ax = fig.add_axes([0.87, 0.1, 0.03, 0.75])


    ax2.set_facecolor('black')
    ax1.set_facecolor('black')
    canvas = FigureCanvasTkAgg(fig, master=root)
    canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=True)
    ax1.xaxis.label.set_color('white')
    ax1.yaxis.label.set_color('white')
    ax1.title.set_color('white')
    ax1.tick_params(axis='x', colors='white')
    ax1.tick_params(axis='y', colors='white')
    ax1.yaxis.get_offset_text().set_color('white')
    colors = [(0, 0, 1), (1, 0, 0)]  # Dark blue to green to yellow
    n_bins = 100  # Number of bins in the colormap
    cmap_name = 'my_custom_cmap'
    cm = LinearSegmentedColormap.from_list(cmap_name, colors, N=n_bins)
    norm = Normalize(vmin=w[1:length,1:width].min(), vmax=w[1:length,1:width].max())
    cbar = None
    sm = ScalarMappable(cmap=cm, norm=norm)
    sm.set_array([])
    cbar = fig.colorbar(sm, shrink=0.5, aspect=6, cax=cbar_ax)
    cbar.ax.yaxis.set_tick_params(color='white')
    cbar.outline.set_edgecolor('white')
    cbar.set_label('Radial distance', color='white')
    cbar.ax.yaxis.get_offset_text().set_color('white')
    plt.setp(cbar.ax.yaxis.get_ticklabels(), color='white')
    surf = None
    plot_step = 5
    anim_step = 5


    def update_plot(i):
        
        nonlocal cbar,surf
            
        ax2.clear()

        R, Theta = np.meshgrid(w[i+1, 1:width+1: plot_step], w[0, 1:width+1: plot_step])
        phi = np.linspace(0, 2 * np.pi, int(width/ plot_step))
        Phi, R = np.meshgrid(phi, w[i+1, 1:width+1: plot_step])
        X = R * np.sin(Theta) * np.cos(Phi)
        Y = R * np.sin(Theta) * np.sin(Phi)
        Z = R * np.cos(Theta)
        color_values=norm(R)
        facecolors = cm(color_values)
        surf = ax2.plot_surface(X, Y, Z,edgecolor='k',linewidth=0.1, antialiased=True,facecolors=facecolors,shade=False)
        #ax2.view_init(elev=0)
        ax2.text2D(0.5, 0.2, f'Rotation period = {omega[i,1]:.2f} hrs', transform=ax2.transAxes, fontsize=12, color='white', ha='center', va='center')
        ax2.set_aspect('equal')
        fig.suptitle(f'Time = {w[i+1,0]/1e+6} Myrs',color='white',fontsize=12)
        ax2.set_axis_off()
        ax2.grid(False)
        sm.set_array(R)
 
        ax1.clear()
        
        x1 =  omega[0:i, 0]/1e+6
        y1 =  omega[0:i, 1]
        ax1.plot(x1, y1, '-b', linewidth=2)

        ax1.spines['top'].set_visible(False)
        ax1.spines['right'].set_visible(False)
        ax1.spines['left'].set_color('white')
        ax1.spines['bottom'].set_color('white')
        ax1.set_xlabel('Time in Myrs',color='white')
        ax1.set_ylabel('Rotation period in hrs',color='white')
        ax1.set_xlim([0,max(omega[:,0])/1e+6])
        ax1.set_ylim([min(omega[:,1]),max(omega[:,1])])

        canvas.draw_idle()

    azim=0
    
    def rotate(frames):
        nonlocal azim
        azim+=int(1/omega[int(frames/anim_step)%length,1]*10)
        azim %= 360
        if ax2:
            ax2.view_init(elev=0,azim=azim)
        if frames%int(anim_step)==0 :
            gui.plot_index+=1
            gui.plot_index%=length
            update_plot(gui.plot_index)
            root.update_idletasks()
        return ax2,
    
    def animate():
        if not gui.is_paused:
            gui.plot_index+=1
            gui.plot_index%=length
            update_plot(gui.plot_index)
            root.update_idletasks()
            
            
        
        elif gui.next_frame:
            gui.plot_index%=length
            update_plot(gui.plot_index)
            root.update_idletasks()
            gui.next_frame = False
        
        elif gui.is_anim:
            gui.progress['maximum'] = int((length-1)*anim_step)
            gui.is_anim = False
            gui.anim = FuncAnimation(fig, rotate, frames=int((length-1)*anim_step), interval=50,blit=False,repeat=False)
            
            #root.after(100,animate) 
        root.after(100,animate)    
        
    root.after(0, animate)
    #gui.anim = FuncAnimation(fig, rotate, frames=360, interval=20,blit=False)
 



    