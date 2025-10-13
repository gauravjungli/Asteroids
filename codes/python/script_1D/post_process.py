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
import os
import seaborn as sns
from IO import Output_File, extract_number
import pdb
from matplotlib.ticker import ScalarFormatter
from Fit import Fit
#from mpl_toolkits.mplot3d import Axes3D

"""
Main post processing script.
"""
def create_fill(ax,x,y,x1,y1):
    
    
    ax.fill_betweenx(y,x,0,color='lightcoral')
    ax.fill_betweenx(y1,x1,0,color='lightblue')
    x=-x
    x1=-x1
    ax.fill_betweenx(y,x,0,color='lightcoral')
    ax.fill_betweenx(y1,x1,0,color='lightblue')
    
    arrow_x_location = 0      # The x-coordinate where the arrow will be
    arrow_y_start = -1.90      # The y-coordinate where the arrow starts
    arrow_length = 3.55     # The length of the arrow

# 3. Draw the vertical arrow using ax.arrow()
    ax.arrow(
        arrow_x_location,      # x-start
        arrow_y_start,         # y-start
        0,                     # dx (change in x) -> 0 for vertical
        arrow_length,          # dy (change in y) -> length of the arrow
        head_width=0.05,        # width of the arrowhead
        head_length=0.05,       # length of the arrowhead
        fc='black',              # face color (fill) of the arrow
        ec='black',              # edge color of the arrow
        linewidth=1,           # width of the arrow line
        zorder=10              # zorder > 0 to ensure it's on top of other elements
        )



def show_shape (parameters):
    #plt.close()
    #pdb.set_trace()
    fig,ax = plt.subplots()
    file1 = parameters['Data folder']
    Res = int(parameters['Resolution'])
    file2 = os.path.join(file1,'dia.txt')
    epsilon = np.loadtxt(file2,dtype=float,ndmin=2)[:,2]
    dia = np.loadtxt(file2,dtype=float,ndmin=2)[:,1]
    

    for count in range(1,1000):        
        #pdb.set_trace()
        file = os.path.join(file1,f'field_{(count)}.csv')
        print(file)
        if os.path.exists(file):
           
            w = np.loadtxt(file,delimiter=",",dtype='float')
            ax.cla()
            theta = w[2:Res-2,0]
            base = w[2:Res-2,1]
            height = epsilon[count-1]*w[2:Res-2,2]
            dbase = w[2:Res-2,3]
        
            metric = base**2 + dbase**2
             
            rad = np.sqrt((2*base**2*np.sqrt(metric)*height + metric*(height**2 + base**2  ))/metric)
            print("maximum radial distance is ", max(rad) )
            z = np.cos(theta)*base + 1/np.sqrt(metric)*(height*dbase*np.sin(theta)+np.cos(theta)*base*height)
            theta_new = np.arccos(z/rad)
            
            x = dia[count-1]*rad*np.sin(theta_new)/2
            y = dia[count-1]*rad*np.cos(theta_new)/2
           # x1 = dia[count-1]*base*np.sin(theta)/2
           # y1 = dia[count-1]*base*np.cos(theta)/2
            ax.plot(x,y,color='black')
            ax.plot(-x,y,color='black')

            #create_fill(ax,x,y,x1,y1)
            # Set x and y limits
            ax.set_xlim(-300, 300)
            ax.set_ylim(-300, 300)
            ax.set_aspect('equal')
            title=f'landslide number={count}'
            #ax.set_title(title)
            #fig.canvas.draw()
            fig.canvas.flush_events()
            plt.show(block=False)
            #fig.canvas.draw_idle()
            plt.pause(0.1)
            plt.savefig(file1+f"/img/img_{count}.svg", dpi=100)


        else:
            break

#%%    



def show_omega(parameters):
    #pdb.set_trace()
    plt.close()
    plt.rcParams.update({'font.size' : 14})
    fig,ax1 = plt.subplots()
    file2='/home/g/Asteroids/output'
    file1=f'Omega_{"C" if parameters["Collision"]=="Yes" else ""}{"L" if parameters["Landslide"]=="Yes" else ""}{"Y" if parameters["YORP"]=="Yes" else ""}.txt'
    file=os.path.join(file2,parameters['Output folder'],f'run{parameters["run"]}',file1)
    if os.path.exists(file):
      
        w=np.loadtxt(file,dtype=float)
        x=w[:,0]
        y=2*np.pi/w[:,1]/3600
        ax1.plot(x/1000,y,'-r',linewidth=2)
        #title="Evolution of time period"
        ax1.set_xlim(0,200)
        ax1.set_ylabel("Time period (hr)")
        ax1.set_xlabel(" Time (ky)")
        
        secax = ax1.secondary_yaxis('right', functions=(lambda y: 2*np.pi/(y*3600), lambda y: 2*np.pi/(y*3600)))
        secax.set_ylabel(r"$\omega$ (1/s)" )
        secax.tick_params(axis='y')
        plt.grid(True)
        #ax.set_title(title)
        fig.canvas.draw()
        fig.canvas.flush_events()
        plt.show(block=False)
        fig.canvas.draw_idle()
        plt.tight_layout()
        # Apply scientific notation formatting
        formatter = ScalarFormatter(useMathText=True)
        formatter.set_scientific(True)
        formatter.set_powerlimits((-2, 2))
        secax.yaxis.set_major_formatter(formatter)
        
        plt.pause(0.1)
    #plt.close(fig)
     
     



 #%%   
def show_individual_run(parameters):
    # Function to extract the numerical value from a name (e.g., "dir_10" -> 10)
    #pdb.set_trace()
    main_dir=parameters['Data folder']
    file1=os.path.join(main_dir,'dia.txt')
    res = int(parameters['Resolution'])
    dia = np.loadtxt(file1,dtype=float,ndmin=2)[:,1]
    epsilon = np.loadtxt(file1,dtype=float,ndmin=2)[:,2]
    Gamma = np.loadtxt(file1,dtype=float,ndmin=2)[:,3]
 
    dx=(math.pi-2*float(parameters["offset"]))/float(parameters["Resolution"])
    fig = plt.figure(figsize=(6,6))  
    subdirs = sorted([d for d in os.listdir(main_dir) if os.path.isdir(os.path.join(main_dir, d))], key=extract_number)
    count=0
    
    for subdir in subdirs:
        
        subdir_path = os.path.join(main_dir, subdir)
        
        # Get and sort files by number within the subdirectory
        dirFiles = sorted( [f for f in os.listdir(subdir_path) if f.lower() != "log.txt"], key=extract_number)
        print(subdir)
        os.chdir(subdir_path)
        ang_mom=[]
        lin_mom=[]
        
        if count<0 or count>12:
            count+=1
            continue
        
        for file in dirFiles:
            

            w=np.loadtxt(file,delimiter=",",dtype=float)
    
            theta  = w[:,0]
            base   = w[:,1]
            height = w[:,2]
            dbase  = w[:,3]
            ddbase = w[:,4]
            u = w[:,5]
            v = w[:,6]
            metric = np.sqrt(base**2 + dbase**2)
            sum1 = 0
            for j in range(2,res-2):
                sum1+= metric[j]*base[j]*height[j]*np.sin(theta[j])*(theta[3]-theta[2])

            print(height[2])

            plt.clf()

            plt.plot(theta[2:res-2],height[2:res-2])
           # plt.plot(theta[2:res-2],w[2:res-2,7])
 
           # plt.plot(x[2:res-2],z1[2:res-2])
            #plt.plot(x, w_filtered)
            plt.title("Time="+str(count))
            plt.pause(0.1) 

            #plt.ylim(0.975,1.025)
           # if count>-1:
           #     break
        count +=1
       # if count>10:
        #    return
        print(count, subdir)
        
        # print(sum(w[:,1]))  
    # ang_mom=np.array(ang_mom)
    # plt.plot(ang_mom[:,0],ang_mom[:,1])
    # lin_mom=np.array(lin_mom)
    # plt.plot(lin_mom[:,0],lin_mom[:,1])
    os.chdir("/home/g/Asteroids/codes/python")
    
#%%

def show_fit(parameters):
   # pdb.set_trace()
   # plt.close()


    fig,ax = plt.subplots()
    file1=parameters['Data folder']
    file2=os.path.join(file1,'dia.txt')
    epsilon=np.loadtxt(file2,dtype=float,ndmin=2)[:,2]
    
    for count in range(1000):
        
        
        file=os.path.join(file1,f'field_{(count+1)}.csv')
        print(file)
        if os.path.exists(file):
           
            w=np.loadtxt(file,delimiter=",",dtype=float)

            ax.cla()#change
            x = np.sin(w[:,0])*w[:,1]
            y = np.cos(w[:,0])*w[:,1]
            x =-np.sin(w[:,0])*w[:,1]
            w_new =  Fit(parameters,epsilon[count],w) 
            continue
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

#%%
def plot_grav(parameters):
   # pdb.set_trace()
    Res = int(parameters['Resolution'])
    fig = plt.figure(figsize=(14,6))
    file =Output_File(parameters,filetype='output',filenames=['base.txt'])
    w=np.loadtxt(file,delimiter=",",dtype=float)
    file1 =Output_File(parameters,filetype='output',filenames=['grav.txt'])
    x=w[:,0]
    grav=np.loadtxt(file1,delimiter=" ")
    #plt.clf()
    plt.plot(x[2:Res-2],grav[2:Res-2,0],linewidth=2,markersize=8)
    plt.plot(x[2:Res-2],grav[2:Res-2,1],linewidth=2,markersize=8)
    plt.grid()
    plt.xlabel(r'$\theta$')
    plt.ylabel('Non-dimensionalized Normal Gravity')
    plt.xlim([0,3.14])
    #plt.ylim([-0.25,0.25])
    plt.minorticks_on()
    plt.tick_params(direction='in',right=True, top=True, left=True, bottom=True)
    plt.tick_params(labelsize=14)
    plt.tick_params(labelbottom=True, labeltop=False, labelright=False, labelleft=True)
    plt.tick_params(direction='in',which='minor', length=5, bottom=True, top=True, left=True, right=True)
    plt.tick_params(direction='in',which='major', length=10, bottom=True, top=True, left=True, right=True)
    
    #plt.savefig('output/grav_20_static.svg', dpi=300,bbox_inches="tight")



#%%
#script for plot 1
def omega_comparison_plot():
    plt.rcParams.update({'font.size' : 14})
    colors=sns.color_palette("rocket",7)
    file1="/home/g/Asteroids/output/PRSA/omega"
    fig = plt.figure(figsize=(6,6)) 
    
    
    #omega=np.loadtxt(file1+"/Omega_L.txt",delimiter="\t")
    #plt.semilogy(omega[0:-1:3,0]/1e+6,omega[0:-1:3,1],linewidth=2,marker='o',mfc='w',markersize=8,color=colors[2],linestyle='None',label='L')
    
    
    omega=np.loadtxt(file1+"/Omega_CY.txt",delimiter="\t")
    plt.semilogy(omega[0:-1:3,0]/1e+6,omega[0:-1:3,1],linewidth=2,marker='o',mfc='w',markersize=8,color=colors[0],linestyle='None',label='CY')
    
    omega=np.loadtxt(file1+"/Omega_CL.txt",delimiter="\t")
    #plt.semilogy(omega[1:-1:3,0]/1e+6,omega[1:-1:3,1],linewidth=2,marker='^',mfc='w',markersize=8,color=colors[4],linestyle='None',label='CL')
    plt.semilogy(omega[:,0]/1e+6,omega[:,1],linewidth=2,markersize=8,color=colors[1],linestyle='solid',label='CL')
    
    omega=np.loadtxt(file1+"/Omega_C.txt",delimiter="\t")
    plt.semilogy(omega[:,0]/1e+6,omega[:,1],linewidth=2,markersize=8,color=colors[2],linestyle='dashdot',label='C')
    
    
    #omega=np.loadtxt(file1+"/Omega_LY.txt",delimiter="\t")
    #plt.semilogy(omega[:,0]/1e+6,omega[:,1],linewidth=2,markersize=8,color=colors[5],linestyle='solid',label='LY')
    
    omega=np.loadtxt(file1+"/Omega_CLY.txt",delimiter="\t")
    plt.semilogy(omega[:,0]/1e+6,omega[:,1],linewidth=2,markersize=8,color=colors[3],linestyle='dotted',label='CLY')
    
    omega=np.loadtxt(file1+"/Omega_Y.txt",delimiter="\t")
    plt.semilogy(omega[:,0]/1e+6,omega[:,1],linewidth=2,markersize=8,color=colors[4],linestyle='dashed',label='Y')
    
    yticks = np.arange(3.0,5,0.25)
    #plt.plot(omega[:,0]/1e+6,(2*math.pi/omega[:,3]/3600),linewidth=2,marker='^',mfc='w',markersize=8,color=colors[0],linestyle='solid',label='L')
    plt.gca().yaxis.set_major_formatter(mtick.ScalarFormatter())
    plt.yticks(yticks)
    plt.xlabel('Time (Myr)')
    plt.ylabel('Time Period (hr)')
    plt.xlim([0,0.25])
    plt.ylim([3.0,4.5])
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
    
    file = Output_File(parameters=parameters, filetype='output', filenames=['mean_shape.csv']) 
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