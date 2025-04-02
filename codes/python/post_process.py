#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import math
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick
from matplotlib.figure import Figure
from matplotlib.backends.backend_agg import FigureCanvasAgg
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import matplotlib.gridspec as gridspec
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.colors import LinearSegmentedColormap, Normalize
from matplotlib.cm import ScalarMappable
from pathlib import Path
import glob
import os
import re
from IO import Parameter
import seaborn as sns
from scipy.ndimage import gaussian_filter1d
from IO import Output_File

#from mpl_toolkits.mplot3d import Axes3D

def show_shape(parameters):
    plt.close()
    fig,ax = plt.subplots()
    epsilon = float(parameters["epsilon"])
    Gamma = float(parameters["Gamma"])
    for count in range(0,100):
        
        file1=parameters['Data folder']
        
        file=os.path.join(file1,f'field_{(count+1)}.csv')
        if os.path.exists(file):
            
            w=np.loadtxt(file,delimiter=",",dtype=float)
           # if count!=slides-1 and count!=0:
            #    continue
            #plt.clf()
            x=np.sin(w[:,0])*(1+epsilon*w[:,2]+Gamma*w[:,1])
            y=np.cos(w[:,0])*(1+epsilon*w[:,2]+Gamma*w[:,1])
            ax.plot(x,y,'-r',linewidth=4)
            x=-np.sin(w[:,0])*(1+epsilon*w[:,2]+Gamma*w[:,1])
            ax.plot(x,y,'-r',linewidth=2)
            ax.set_aspect('equal')
            title=f'landslide number={count+1}'
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
     
     



def backup():
    plt.rcParams.update({'font.size' : 14})
    colors=sns.color_palette("rocket",7)
    
    parameters={}
    Parameter(parameters)
    omega=float(parameters["omega_in"])
    delta=float(parameters["Friction angle"])
    slides=int(parameters["slides"])
    epsilon=float(parameters["epsilon"])
    Gamma=float(parameters["Gamma"])
    res=int(parameters["res"])
    offset=float(parameters["offset"])
    dx=(math.pi-2*offset)/res
    #omega=0.65
    #delta=30

#file1="output/omega_15_0.65_0.002"
# omega=np.loadtxt(file1+"/omega.txt",delimiter=" ")
#%%
def show_plot(parameters):
    fig = plt.figure(figsize=(6,6))
    
    for count in range(0,int(parameters["slides"])):
        file1=parameters['Data folder']
        file=glob.glob(file1+"/field_"+str(count+1)+".csv",recursive=True)
        w=np.loadtxt(file[0],delimiter=",",dtype=float)
        print(file)
       # if count!=slides-1 and count!=0:
        #    continue

        x=w[:,0]#np.sin(w[:,0])*(1+(float(parameters["epsilon"])*w[:,2]+float(parameters["Gamma"])*w[:,1]))
        y=w[:,2]#np.cos(w[:,0])*(1+(float(parameters["epsilon"])*w[:,2]+float(parameters["Gamma"])*w[:,1]))
        plt.clf()
        #plt.axis('equal')
        plt.plot(x,y,'-r',linewidth=1)
       # x=-np.sin(w[:,0])*(1+(float(parameters["Gamma"])*w[:,1]+float(parameters["Gamma"])*w[:,1]))
       # plt.plot(x,y,'-r',linewidth=4)
        plt.title("lanslide number="+str(count+1))
        plt.pause(0.5)
        #plt.savefig(file1+"/img_"+str(count+1)+".svg",dpi=300,bbox_inches="tight")
    
#plt.close()
 #%%   
def show_individual_run(parameters):
    # Function to extract the numerical value from a name (e.g., "dir_10" -> 10)
    def extract_number(text):
        match = re.search(r'\d+', text)  # Find first occurrence of a number
        return int(match.group()) if match else float('inf')  # Default to large number if no match

    main_dir=parameters['Data folder']
    file1=os.path.join(main_dir,'dia.txt')

    epsilon=np.loadtxt(file1,dtype=float)[:,2]
    Gamma = np.loadtxt(file1,dtype=float)[:,3]
 
    dx=(math.pi-2*float(parameters["offset"]))/float(parameters["Resolution"])
    fig = plt.figure(figsize=(6,6))  
    subdirs = sorted([d for d in os.listdir(main_dir) if os.path.isdir(os.path.join(main_dir, d))], key=extract_number)
    count=0
    
    for subdir in subdirs:
        
        subdir_path = os.path.join(main_dir, subdir)
        
        # Get and sort files by number within the subdirectory
        dirFiles = sorted( [f for f in os.listdir(subdir_path) if f.lower() != "log.txt"], key=extract_number)
   
        os.chdir(subdir_path)
        ang_mom=[]
        lin_mom=[]
        
        if count <0:
            count +=1
            continue
        for file in dirFiles:
            

            w=np.loadtxt(file,delimiter=",",dtype=float)
          #  print(file)
            w_filtered = gaussian_filter1d(w[:,2], sigma=20)
            plt.clf()
            x=(w[:,0])
            y=(w[:,3])#+Gamma[count]/epsilon[count]*w[:,1])
            # y=(w[:,3])
            ang_mom.append([count,sum(w[:,4])*dx])
            lin_mom.append([count,sum(w[:,3])*dx])
            
            #plt.plot(x,w_filtered)
            plt.plot(x,y)
            #plt.plot(x, w_filtered)
            plt.title("Time="+str(count))
            plt.pause(0.1) 
           # if count > -1:
            #    break
            #plt.ylim(0.975,1.025)
            
        count +=1
        print(count, subdir)
        if count>4:
            return
        # print(sum(w[:,1]))  
    # ang_mom=np.array(ang_mom)
    # plt.plot(ang_mom[:,0],ang_mom[:,1])
    # lin_mom=np.array(lin_mom)
    # plt.plot(lin_mom[:,0],lin_mom[:,1])
    os.chdir("../..")
#%%
def plot_grav():
    fig = plt.figure(figsize=(14,6)) 
    x=w[:,0]
    grav=np.loadtxt(file1+"/grav.txt",delimiter=" ")
    #plt.clf()
    plt.plot(x,grav[:,1],linewidth=2,markersize=8)
    
    plt.grid()
    plt.xlabel(r'$\theta$')
    plt.ylabel('Non-dimensionalized Normal Gravity')
    plt.xlim([0,3.14])
    plt.ylim([-0.25,0.25])
    plt.minorticks_on()
    plt.tick_params(direction='in',right=True, top=True, left=True, bottom=True)
    plt.tick_params(labelsize=14)
    plt.tick_params(labelbottom=True, labeltop=False, labelright=False, labelleft=True)
    plt.tick_params(direction='in',which='minor', length=5, bottom=True, top=True, left=True, right=True)
    plt.tick_params(direction='in',which='major', length=10, bottom=True, top=True, left=True, right=True)
    
    plt.savefig('output/grav_20_static.svg', dpi=300,bbox_inches="tight")


#%%
#fig = plt.figure(figsize=(14,6)) 
# omega=np.loadtxt(file1+"/omega_L.txt",delimiter="\t")


# #plt.plot(omega[:,0]/1e+6,(2*math.pi/omega[:,3]/3600),linewidth=2,marker='^',mfc='w',markersize=8,color=colors[0],linestyle='solid',label='L')
# plt.plot(omega[:,0]/1e+6,(2*math.pi/omega[:,3]/3600),linewidth=2,markersize=8,color=colors[2],linestyle='dotted',label='L')

# plt.xlabel('Time (Myr)')
# plt.ylabel('Time Period (hr)')
# plt.xlim([0,1])
# plt.ylim([2.5,4.5])
# plt.minorticks_on()
# plt.tick_params(direction='in',right=True, top=True, left=True, bottom=True)
# plt.tick_params(labelsize=14)
# plt.tick_params(labelbottom=True, labeltop=False, labelright=False, labelleft=True)
# plt.tick_params(direction='in',which='minor', length=5, bottom=True, top=True, left=True, right=True)
# plt.tick_params(direction='in',which='major', length=10, bottom=True, top=True, left=True, right=True)
# plt.legend(loc='best')
# plt.legend(fontsize=14) 
# plt.savefig('output/Omega.svg', dpi=300,bbox_inches="tight")

#%%
#script for plot 1
def omega_comparison_plot():
    plt.rcParams.update({'font.size' : 14})
    colors=sns.color_palette("rocket",7)
    file1="/home/g/Asteroids/output/omega"
    fig = plt.figure(figsize=(6,6)) 
    
    
    omega=np.loadtxt(file1+"/Omega_L.txt",delimiter="\t")
    plt.semilogy(omega[0:-1:3,0]/1e+6,omega[0:-1:3,1],linewidth=2,marker='o',mfc='w',markersize=8,color=colors[2],linestyle='None',label='L')
    
    
    omega=np.loadtxt(file1+"/Omega_CY.txt",delimiter="\t")
    plt.semilogy(omega[0:-1:3,0]/1e+6,omega[0:-1:3,1],linewidth=2,marker='s',mfc='w',markersize=8,color=colors[0],linestyle='None',label='CY')
    
    omega=np.loadtxt(file1+"/Omega_CL.txt",delimiter="\t")
    plt.semilogy(omega[1:-1:3,0]/1e+6,omega[1:-1:3,1],linewidth=2,marker='^',mfc='w',markersize=8,color=colors[4],linestyle='None',label='CL')
    
    omega=np.loadtxt(file1+"/Omega_C.txt",delimiter="\t")
    plt.semilogy(omega[:,0]/1e+6,omega[:,1],linewidth=2,markersize=8,color=colors[0],linestyle='dashdot',label='C')
    
    
    omega=np.loadtxt(file1+"/Omega_LY.txt",delimiter="\t")
    plt.semilogy(omega[:,0]/1e+6,omega[:,1],linewidth=2,markersize=8,color=colors[5],linestyle='solid',label='LY')
    
    omega=np.loadtxt(file1+"/Omega_CLY.txt",delimiter="\t")
    plt.semilogy(omega[:,0]/1e+6,omega[:,1],linewidth=2,markersize=8,color=colors[1],linestyle='dotted',label='CLY')
    
    omega=np.loadtxt(file1+"/Omega_Y.txt",delimiter="\t")
    plt.semilogy(omega[:,0]/1e+6,omega[:,1],linewidth=2,markersize=8,color=colors[6],linestyle='dashed',label='Y')
    
    yticks = np.arange(2.5,6.5,0.5)
    #plt.plot(omega[:,0]/1e+6,(2*math.pi/omega[:,3]/3600),linewidth=2,marker='^',mfc='w',markersize=8,color=colors[0],linestyle='solid',label='L')
    plt.gca().yaxis.set_major_formatter(mtick.ScalarFormatter())
    plt.yticks(yticks)
    plt.xlabel('Time (Myr)')
    plt.ylabel('Time Period (hr)')
    plt.xlim([0,0.5])
    plt.ylim([3.0,6.0])
    plt.minorticks_on()
    plt.tick_params(direction='in',right=True, top=True, left=True, bottom=True)
    plt.tick_params(labelsize=14)
    plt.tick_params(labelbottom=True, labeltop=False, labelright=False, labelleft=True)
    plt.tick_params(direction='in',which='minor', length=5, bottom=True, top=True, left=True, right=True)
    plt.tick_params(direction='in',which='major', length=10, bottom=True, top=True, left=True, right=True)
    plt.legend(loc='best')
    plt.legend(fontsize=14) 
    plt.savefig(file1+'/Omega.svg', dpi=300,bbox_inches="tight")



#%%

#script for plot 3
def plot4():
    file1="/home/g/Documents/Asteroids/output/saved_data/files_15.000000_0.650000_C"
    
    fig = plt.figure(figsize=(6,6)) 
    omega=np.loadtxt(file1+"/omega_T.txt",delimiter="\t")
    plt.semilogy(omega[:,0]/1e+6,((2*math.pi/omega[:,1]/3600)),linewidth=2,markersize=8,color=colors[0],linestyle='dashdot',label='Y')
    
    
    
    
    yticks = [3,4,5,10,15,20,25]
    #plt.plot(omega[:,0]/1e+6,(2*math.pi/omega[:,3]/3600),linewidth=2,marker='^',mfc='w',markersize=8,color=colors[0],linestyle='solid',label='L')
    plt.yticks(yticks)
    plt.gca().yaxis.set_major_formatter(mtick.ScalarFormatter())
    plt.gca().yaxis.set_minor_formatter(mtick.NullFormatter())
    plt.xlabel('Time (Myr)')
    plt.ylabel('Time Period (hr)')
    plt.xlim([0,5])
    plt.ylim([2.8,25.5])
    plt.minorticks_on()
    plt.tick_params(direction='in',right=True, top=True, left=True, bottom=True)
    plt.tick_params(labelsize=14)
    plt.tick_params(labelbottom=True, labeltop=False, labelright=False, labelleft=True)
    plt.tick_params(direction='in',which='minor', length=5, bottom=True, top=True, left=True, right=True)
    plt.tick_params(direction='in',which='major', length=10, bottom=True, top=True, left=True, right=True)
    plt.legend(loc='best')
    plt.legend(fontsize=14) 
    plt.savefig('output/Omega_3.svg', dpi=300,bbox_inches="tight")
    
    
def save_3D_omega(parameters):

    file1=f'Omega_{"C" if parameters["Collision"]=="Yes" else ""}{"L" if parameters["Landslide"]=="Yes" else ""}{"Y" if parameters["YORP"]=="Yes" else ""}.txt'
    file = Output_File(parameters=parameters, filetype='output', filenames=[file1]) 
    omega = np.loadtxt(file, dtype=float)
    
    length, width = np.shape(omega)

    fig = plt.figure(figsize=(12, 5), dpi=150) 
    fig.patch.set_facecolor('black')
    ax1 = fig.add_axes([0.10, 0.12, 0.85,0.85])

    ax1.set_facecolor('black')

    ax1.xaxis.label.set_color('white')
    ax1.yaxis.label.set_color('white')
    ax1.title.set_color('white')
    ax1.tick_params(axis='x', colors='white')
    ax1.tick_params(axis='y', colors='white')
    ax1.tick_params(direction='in',which='minor', length=5, bottom=True, top=True, left=True, right=True)
    ax1.tick_params(direction='in',which='major', length=10, bottom=True, top=True, left=True, right=True)
    ax1.yaxis.get_offset_text().set_color('white')
    

    

    def update_plot(i):
        
        ax1.clear()
        
        x1 =  omega[0:i, 0]/1e+6
        y1 =  omega[0:i, 1]
        z1 = omega[0:i, 1]- omega[0:i,2]
        z2 = omega[0:i, 1]+ omega[0:i,2]
        ax1.plot(x1, y1, '-b', linewidth=2)
        ax1.plot(x1,z1,color='gray',linewidth=2)
        ax1.plot(x1,z2,color='gray',linewidth=2)
        ax1.spines['top'].set_visible(True)
        ax1.spines['right'].set_visible(True)
        ax1.spines['left'].set_color('white')
        ax1.spines['bottom'].set_color('white')
        ax1.spines['right'].set_color('white')
        ax1.spines['top'].set_color('white')
        ax1.set_xlabel('Simulation time (Myrs)',color='white')
        ax1.set_ylabel('Rotational time period  (hrs)',color='white')
        ax1.set_xlim([0,1.0])#max(omega[:,0])/1e+6])
        ax1.set_ylim([3.8,5.3])
        img = Output_File(parameters=parameters, filetype='output', filenames=['plots','omega.svg']) 
        fig.savefig(img, dpi=150, bbox_inches='tight')  # Save as PNG
   #    canvas.draw_idle()
        #canvas = FigureCanvasAgg(fig)
        #canvas.draw()
        #plt.imshow(canvas.buffer_rgba()) 
        plt.show()
    for i in range(length):
        if i%(length-1) ==0:
            #print(i)
           # print(omega[i,0])
            update_plot(i)
    
def save_3D_shapes(parameters):
    
    file = Output_File(parameters=parameters, filetype='output', filenames=['mean_shape.txt']) 
    w = np.loadtxt(file, dtype=float, delimiter=",")
    file1=f'Omega_{"C" if parameters["Collision"]=="Yes" else ""}{"L" if parameters["Landslide"]=="Yes" else ""}{"Y" if parameters["YORP"]=="Yes" else ""}.txt'
    file = Output_File(parameters=parameters, filetype='output', filenames=[file1]) 
    omega = np.loadtxt(file, dtype=float)
    
    length, width = np.shape(w)
    length -= 1       # First row of mean shape contains theta and first column contains time   
    width  -= 1
    fig1 = plt.figure(figsize=(7, 5), dpi=150) 
    fig1.patch.set_facecolor('black')
    ax2 = fig1.add_axes([0, 0.2, 0.7, 0.7], projection='3d') 
    cbar_ax = fig1.add_axes([0.80, 0.1, 0.03, 0.75])


    ax2.set_facecolor('black')

   
    colors = [(0, 0, 1), (1, 0, 0)]  # Dark blue to green to yellow
    n_bins = 50  # Number of bins in the colormap
    cmap_name = 'my_custom_cmap'
    cm = LinearSegmentedColormap.from_list(cmap_name, colors, N=n_bins)
    norm = Normalize(vmin=w[1:length,1:width].min(), vmax=w[1:length,1:width].max())
    cbar = None
    sm = ScalarMappable(cmap=cm, norm=norm)
    sm.set_array([])
    cbar = fig1.colorbar(sm, shrink=0.5, aspect=6, cax=cbar_ax)
    cbar.ax.yaxis.set_tick_params(color='white')
    cbar.outline.set_edgecolor('white')
    cbar.set_label('Radial distance', color='white')
    cbar.ax.yaxis.get_offset_text().set_color('white')
    plt.setp(cbar.ax.yaxis.get_ticklabels(), color='white')
    surf = None
    plot_step = 5
    i=0

    

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

        ax2.text2D(0.5, 0.2, f'Rotation period = {omega[i,1]:.2f} hrs', transform=ax2.transAxes, fontsize=12, color='white', ha='center', va='center')
        ax2.set_aspect('equal')
        fig1.suptitle(f'Time = {w[i+1,0]/1e+6} Myrs',color='white')
        ax2.set_axis_off()
        ax2.grid(False)
        sm.set_array(R)
 
       
        img = Output_File(parameters=parameters, filetype='output', filenames=['plots',f'shape{i}.svg']) 
        fig1.savefig(img, dpi=150, bbox_inches='tight')  # Save as PNG
   #    canvas.draw_idle()
        #canvas = FigureCanvasAgg(fig)
        #canvas.draw()
        #plt.imshow(canvas.buffer_rgba()) 
        plt.show()
    for i in range(length):
        if i%(1) ==0:
            update_plot(i) 

def show_shape_deviation(parameters):
    plt.close()
    fig,ax = plt.subplots()
    file = Output_File(parameters=parameters, filetype='output', filenames=['mean_shape.txt']) 
    w = np.loadtxt(file, dtype=float, delimiter=",")
    file = Output_File(parameters=parameters, filetype='output', filenames=['std_shape.txt']) 
    std = np.loadtxt(file, dtype=float, delimiter=",")
    file1=f'Omega_{"C" if parameters["Collision"]=="Yes" else ""}{"L" if parameters["Landslide"]=="Yes" else ""}{"Y" if parameters["YORP"]=="Yes" else ""}.txt'
    file = Output_File(parameters=parameters, filetype='output', filenames=[file1]) 
    omega = np.loadtxt(file, dtype=float)
    length, width = np.shape(w)
    length -= 1       # First row of mean shape contains theta and first column contains time   
    width  -= 1
    for i in range(int(length/2)):
        R1 = w[i+1, 1:width+1:]#-std[i+1, 1:width+1:]
        R2 = w[i+1, 1:width+1:]#+std[i+1, 1:width+1:]
        Theta =  w[0, 1:width+1]



        ax.cla()
        x1=np.sin(Theta)*R1
        x2=np.sin(Theta)*R2
        y1=np.cos(Theta)*R1
        y2=np.cos(Theta)*R2
        ax.plot(x1,y1,'-',color='red',linewidth=2)
        #ax.plot(x2,y2,'-',linewidth=1,color='gray')
        x1=-np.sin(Theta)*R1
        x2=-np.sin(Theta)*R2
        ax.plot(x1,y1,'-r',linewidth=2,color='red')
        #ax.plot(x2,y2,'-r',linewidth=1, color='gray')
        ax.set_aspect('equal')
        title=f'landslide number={i+1}'
        ax.set_title(title)
        fig.canvas.draw()
        #fig.canvas.flush_events()
        plt.show(block=False)
        #fig.canvas.draw_idle()
        plt.pause(0.1)
        plt.savefig("/home/g/Asteroids/plots_stoc"+"/img_"+str(i)+".svg",dpi=300,bbox_inches="tight")
    #plt.pause(10)
    #plt.close(fig)
    
    
def save_2D_omega(parameters):

    file1=f'Omega_{"C" if parameters["Collision"]=="Yes" else ""}{"L" if parameters["Landslide"]=="Yes" else ""}{"Y" if parameters["YORP"]=="Yes" else ""}.txt'
    file = Output_File(parameters=parameters, filetype='output', filenames=[file1]) 
    omega = np.loadtxt(file, dtype=float)
    
    length, width = np.shape(omega)

    fig = plt.figure(figsize=(12, 5), dpi=150) 
    fig.patch.set_facecolor('white')
    ax1 = fig.add_axes([0.10, 0.12, 0.85,0.85])

    ax1.set_facecolor('white')

    ax1.xaxis.label.set_color('black')
    ax1.yaxis.label.set_color('black')
    ax1.title.set_color('black')
    ax1.tick_params(axis='x', colors='black')
    ax1.tick_params(axis='y', colors='black')
    ax1.tick_params(direction='in',which='minor', length=5, bottom=True, top=True, left=True, right=True)
    ax1.tick_params(direction='in',which='major', length=10, bottom=True, top=True, left=True, right=True)
    ax1.yaxis.get_offset_text().set_color('black')
    

    

    def update_plot(i):
        
        ax1.clear()
        
        x1 =  omega[0:i, 0]/1e+6
        y1 =  omega[0:i, 1]
        z1 = omega[0:i, 1]- omega[0:i,2]
        z2 = omega[0:i, 1]+ omega[0:i,2]
        ax1.plot(x1, y1, '-b', linewidth=2)
        ax1.plot(x1,z1,color='gray',linewidth=2)
        ax1.plot(x1,z2,color='gray',linewidth=2)
        ax1.spines['top'].set_visible(True)
        ax1.spines['right'].set_visible(True)
        ax1.spines['left'].set_color('black')
        ax1.spines['bottom'].set_color('black')
        ax1.spines['right'].set_color('black')
        ax1.spines['top'].set_color('black')
        ax1.set_xlabel('Simulation time (Myrs)',color='black')
        ax1.set_ylabel('Rotational time period  (hrs)',color='black')
        ax1.set_xlim([0,1.0])#max(omega[:,0])/1e+6])
        ax1.set_ylim([3.0,6.0])
        img = Output_File(parameters=parameters, filetype='output', filenames=['plots','omega.svg']) 
        fig.savefig(img, dpi=150, bbox_inches='tight')  # Save as PNG
   #    canvas.draw_idle()
        #canvas = FigureCanvasAgg(fig)
        #canvas.draw()
        #plt.imshow(canvas.buffer_rgba()) 
        plt.show()
    for i in range(length):
        if i%(length-1) ==0:
            #print(i)
           # print(omega[i,0])
            update_plot(i)