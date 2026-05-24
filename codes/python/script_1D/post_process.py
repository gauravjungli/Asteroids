#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import math
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick
from matplotlib.colors import LinearSegmentedColormap, Normalize
from matplotlib.cm import ScalarMappable
import os
import seaborn as sns
from IO import Output_File, extract_number
import pdb
from matplotlib.ticker import ScalarFormatter
from Fit import Fit, Fit_radius
from collisions import G
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
            #Fit(parameters,base_old=w)
            ax.cla()
            #continue
            
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
            ax.set_xlim(-350, 350)
            ax.set_ylim(-350, 350)
            ax.set_aspect('equal')
            ax.xaxis.set_visible(False)
            ax.yaxis.set_visible(False)
            title=f'landslide number={count}'
            #ax.set_title(title)
            #fig.canvas.draw()
            fig.canvas.flush_events()
            plt.show(block=False)
            #plt.grid()
            #fig.canvas.draw_idle()
            plt.pause(0.5)
           # plt.savefig(file1+f"/img/img_{count}.svg", dpi=300)


        else:
            break
#%%

def show_height (parameters):
    plt.close()
    pdb.set_trace()
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
            Fit(parameters,base_old=w)
            ax.cla()
            continue
            
            theta = w[2:Res-2,0]
            base = w[2:Res-2,1]
            height = w[2:Res-2,2]
            dbase = w[2:Res-2,3]
            psi = w[2:Res-2,7]
            x = theta
            y = dbase

            ax.plot(x,y,color='black')
            


            title=f'landslide number={count}'
            ax.set_title(title)
            #fig.canvas.draw()
            fig.canvas.flush_events()
            plt.show(block=False)
            #plt.grid()
            #fig.canvas.draw_idle()
            plt.pause(0.5)


        else:
            break

#%%    

def show_omega(parameters):
    pdb.set_trace()
    #plt.close()
    plt.rcParams.update({'font.size' : 17})
   # fig,ax1 = plt.subplots()
    file2='/home/g/Asteroids/output'
    file1=f'Omega_{"C" if parameters["Collision"]=="Yes" else ""}{"L" if parameters["Landslide"]=="Yes" else ""}{"Y" if parameters["YORP"]=="Yes" else ""}.txt'
    file=os.path.join(file2,parameters['Output folder'],f'run{parameters["run"]}',file1)
    if os.path.exists(file):
      
        w=np.loadtxt(file,dtype=float)
        x=w[:,0]
        y=2*np.pi/w[:,1]/3600
        plt.plot(x/1e+6,y,'-r',linewidth=2)
        return
        #title="Evolution of time period"
        ax1.set_xlim(0,2.5)
        ax1.set_ylabel("Time period (hr)")
        ax1.set_xlabel(" Time (Myr)")
        
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
    plt.rcParams.update({'font.size' : 14})
    colors=sns.color_palette("rocket",7)
    main_dir=parameters['Data folder']
    file1=os.path.join(main_dir,'dia.txt')
    res = int(parameters['Resolution'])
    dia = np.loadtxt(file1,dtype=float,ndmin=2)[:,1]
    epsilon = 0.0075#np.loadtxt(file1,dtype=float,ndmin=2)[:,2]
    Gamma = np.loadtxt(file1,dtype=float,ndmin=2)[:,3]
 
    dx=(math.pi-2*float(parameters["offset"]))/float(parameters["Resolution"])
    fig = plt.figure(figsize=(6,6))  
    subdirs = sorted([d for d in os.listdir(main_dir) if os.path.isdir(os.path.join(main_dir, d))], key=extract_number)
    count=0
    
    for subdir in subdirs:
        
        subdir_path = os.path.join(main_dir, subdir)
        dens =1250
        # Get and sort files by number within the subdirectory
        omega = 2*np.pi/(float(parameters["Rotation period"])*3600)/(G * (4/3) * np.pi * dens)**0.5
        print(omega)
        omega =  0.295239
        dirFiles = sorted( [f for f in os.listdir(subdir_path) if f.lower() != "log.txt"], key=extract_number)
        print(omega)
        print(subdir)
        os.chdir(subdir_path)
        ang_mom=[]
        lin_mom=[]
        
        if count<0 or count>1000:
            count+=1
            continue
        max_sum2=0
        min_sum2 = 1e+12
        count1 = 0
        N = 1000
        TH = np.zeros(N)
        TE = np.zeros(N)
        KE = np.zeros(N)
        PE = np.zeros(N)
        AM = np.zeros(N)
        AMS =np.zeros(N)
        time = np.zeros(N)
        for file in dirFiles:
     
      
            w=np.loadtxt(file,delimiter=",",dtype=float)
    
            theta  = w[2:res-2,0]
            base   = w[2:res-2,1]
            height = w[2:res-2,2]
            dbase  = w[2:res-2,3]
            ddbase = w[2:res-2,4]
            u = w[2:res-2,5]
            v = w[2:res-2,6]
            psi = w[2:res-2,7]
            metric = np.sqrt(base**2 + dbase**2)
            J =  metric*base*np.sin(theta)*(1+epsilon*height)
            R = base*np.sin(theta)*(1+epsilon*height/2)

                
            TE[count1] = np.trapezoid(J*height*(u**2+(v+omega*R)**2 + epsilon*height),theta)/2
            PE[count1] = np.trapezoid(J*height*(epsilon*height),theta)/2
            KE[count1] = np.trapezoid(J*height*(u**2+(v+omega*R)**2),theta)/2
            AM[count1] = np.trapezoid(J*R*(v+omega*R)*height,theta)
            TH[count1] =  np.trapezoid(J*height,theta)
            time[count1] = 0.000628318530717958618*100*(count1)

            plt.clf()
           # plt.title("Time="+str(time[count1]))
            count1+=1
           # if count1<N:
            #    continue
            plt.plot(theta,psi,linewidth=2,linestyle ='-',color='r')
            #plt.plot(theta[2:res-2],height[2:res-2]*u[2:res-2])
           # plt.plot(theta[2:res-2],w[2:res-2,7])
 
           # plt.plot(x[2:res-2],z1[2:res-2])
            #plt.plot(x, w_filtered)
            
            plt.grid()
            plt.tick_params(labelsize=20)
            plt.xlabel('Theta', fontsize=20)
            plt.ylabel('Height', fontsize=20)
            plt.xlim(0,3.14)
            #plt.savefig(f"/home/g/Asteroids/output/check_5/run1/img/img_height{count1}.svg", dpi=300)
            plt.pause(0.1) 
            
            
        count +=1
        
        plot_flag = False
        if plot_flag:
            plt.plot(time,TH/max(TH),linewidth=2,color='r',markersize=8,linestyle='dashed',label=' Mass')
            plt.plot(time,TE/max(TE),linewidth=2,markersize=8,color='b',linestyle='solid',label='Total energy')
            plt.plot(time,KE/max(TE),linewidth=2,markersize=8,color=colors[2],linestyle='dashdot',label='Kinetic energy')
            plt.plot(time,PE/max(TE),linewidth=2,markersize=8,color=colors[1],linestyle='dotted',label='Potential energy')
            plt.plot(time[::20], AM[::20]/max(AM),linewidth=2,marker='o',mfc='w',markersize=4,color=colors[3],linestyle='None',label='Angular momentum')
            plt.xlabel('Non-dimensional time', fontsize=20)
            plt.ylabel('Non-dimensional quantities', fontsize=20)
            plt.legend(fontsize=20)
            plt.tick_params(labelsize=20)
            plt.xlim(0,200)
            plt.grid()
       # if count>10:
        #    return
        print(count, subdir)
       # print("The ratio of energy is: ", 1 - min_sum2/max_sum2)
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
    plt.rcParams.update({'font.size' : 18})
    colors=sns.color_palette("rocket",7)
    file1="/home/g/Asteroids/output/thesis/comparison_4"
    fig,ax = plt.subplots(figsize=(16,9),dpi=150) 
    
    
    #omega=np.loadtxt(file1+"/Omega_L.txt",delimiter="\t")
    #plt.semilogy(omega[0:-1:3,0]/1e+6,omega[0:-1:3,1],linewidth=2,marker='o',mfc='w',markersize=8,color=colors[2],linestyle='None',label='L')
    
    
    omega=np.loadtxt(file1+"/Omega_CY.txt",delimiter="\t")
    plt.plot(omega[0:-1:3,0]/1e+6,omega[0:-1:3,1],linewidth=3,marker='o',mfc='w',markersize=8,color=colors[0],linestyle='None',label='CY')
    
    omega=np.loadtxt(file1+"/Omega_CL.txt",delimiter="\t")
    #plt.semilogy(omega[1:-1:3,0]/1e+6,omega[1:-1:3,1],linewidth=2,marker='^',mfc='w',markersize=8,color=colors[4],linestyle='None',label='CL')
    plt.plot(omega[:,0]/1e+6,omega[:,1],linewidth=3,markersize=8,color=colors[1],linestyle='solid',label='CL')
    
    omega=np.loadtxt(file1+"/Omega_C.txt",delimiter="\t")
    plt.plot(omega[:,0]/1e+6,omega[:,1],linewidth=3,markersize=8,color=colors[2],linestyle='dashdot',label='C')
    
    
    #omega=np.loadtxt(file1+"/Omega_LY.txt",delimiter="\t")
    #plt.semilogy(omega[:,0]/1e+6,omega[:,1],linewidth=2,markersize=8,color=colors[5],linestyle='solid',label='LY')
    
    omega=np.loadtxt(file1+"/Omega_CLY.txt",delimiter="\t")
    plt.plot(omega[:,0]/1e+6,omega[:,1],linewidth=3,markersize=8,color=colors[3],linestyle='dotted',label='CLY')
    
    omega=np.loadtxt(file1+"/Omega_Y.txt",delimiter="\t")
    plt.plot(omega[:,0]/1e+6,omega[:,1],linewidth=3,markersize=8,color=colors[4],linestyle='dashed',label='Y')
    
    # 3. Define the threshold for the horizontal line
    threshold = 6
    # 5. Color the regions
    # axhspan(ymin, ymax, ...) fills the area between two y-values across the whole x-axis
    # Region below the line
    plt.axhspan(ymin=0, ymax=threshold, facecolor='lightblue', alpha=0.5)

    # Region above the line
    plt.axhspan(ymin=threshold, ymax=10, facecolor='lightgreen', alpha=0.5 )

# 4. Draw the horizontal line parallel to the x-axis
    #plt.axhline(y=threshold, color='black', linestyle='-')
    
    # plt.text(1.5, 6.5, "Landslides not significant", 
    #       color='black',
    #      ha='center', va='center')

    # Text for the Lower Region (centered horizontally at 5, vertically at 2.5)
    # plt.text(0.75, 5, "Landslides are important", 
    #       color='black',
    #      ha='center', va='center')
    
    yticks = np.arange(2,10,0.5)
    plt.yticks(yticks)
    
    #plt.plot(omega[:,0]/1e+6,(2*math.pi/omega[:,3]/3600),linewidth=2,marker='^',mfc='w',markersize=8,color=colors[0],linestyle='solid',label='L')
    plt.gca().yaxis.set_major_formatter(mtick.ScalarFormatter())
    plt.xlabel('Time (Myr)')
    plt.ylabel('Time Period (hr)')
    plt.xlim([0,2])
    plt.ylim([2.9,5.2])
    secax = ax.secondary_yaxis('right', functions=(lambda y: 2*np.pi/(y*3600), lambda y: 2*np.pi/(y*3600)))
    secax.set_ylabel(r"$\omega$ (1/s)" )
    secax.tick_params(axis='y')
    #secax.set_yticks(ax.get_yticks())
    plt.minorticks_on()
    secax.minorticks_on()
    # Apply scientific notation formatting
    formatter = ScalarFormatter(useMathText=True)
    formatter.set_scientific(True)
    formatter.set_powerlimits((-2, 2))
    secax.yaxis.set_major_formatter(formatter)
    secax.tick_params(direction='in',which='minor', length=5, width=1, color='black', bottom=False, top=False, left=False, right=True)
    secax.tick_params(direction='in',which='major', length=10, width=1, color='black', bottom=False, top=False, left=False, right=True)
    plt.tick_params(labelbottom=True, labeltop=False, labelright=False, labelleft=True)
    plt.tick_params(direction='in',which='minor', length=5, bottom=True, top=True, left=True, right=False)
    plt.tick_params(direction='in',which='major', length=10, bottom=True, top=True, left=True, right=False)
    plt.legend(loc='best')
    plt.savefig(file1+'/Omega.svg', dpi=300,bbox_inches="tight")


#%%

    
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


    colors = [(0, 0, 1),(0,0.5,0), (1, 0, 0)]  # Dark blue to green to yellow
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
        surf = ax2.plot_surface(X, Y, Z,linewidth=0.1, antialiased=True,facecolors=facecolors,shade=False)

        ax2.text2D(0.5, 0.2, f'Rotation period = {omega[i,1]:.2f} hrs', transform=ax2.transAxes, fontsize=12, color='white', ha='center', va='center')
        ax2.set_aspect('equal')
        fig1.suptitle(f'Time = {w[i+1,0]/1e+6} Myrs',color='white')
        ax2.set_axis_off()
        ax2.grid(False)
        ax2.view_init(elev=20)
        sm.set_array(R)
 
       
        img = Output_File(parameters=parameters, filetype='output', filenames=['plots',f'shape{i}.svg']) 
        fig1.savefig(img, dpi=300, bbox_inches='tight')  # Save as PNG
   #    canvas.draw_idle()
        #canvas = FigureCanvasAgg(fig)
        #canvas.draw()
        #plt.imshow(canvas.buffer_rgba()) 
        plt.show()
    for i in range(length):
        if i%(20) ==0:
            update_plot(i) 

def show_shape_deviation(parameters):
    plt.close()
    fig,ax = plt.subplots()
    file = Output_File(parameters=parameters, filetype='output', filenames=['mean_shape.csv']) 
    w = np.loadtxt(file, dtype=float, delimiter=",")
    file = Output_File(parameters=parameters, filetype='output', filenames=['std_shape.csv']) 
    std = np.loadtxt(file, dtype=float, delimiter=",")
    file1=f'Omega_{"C" if parameters["Collision"]=="Yes" else ""}{"L" if parameters["Landslide"]=="Yes" else ""}{"Y" if parameters["YORP"]=="Yes" else ""}.txt'
    file = Output_File(parameters=parameters, filetype='output', filenames=[file1]) 
    omega = np.loadtxt(file, dtype=float)
    length, width = np.shape(w)
    length -= 1       # First row of mean shape contains theta and first column contains time   
    width  -= 1
    for i in range(int(length)):
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
        ax.set_xlim([-300,300])#max(omega[:,0])/1e+6])
        ax.set_ylim([-300,300])
        ax.set_aspect('equal')
        title=f'landslide number={i+1}'
        ax.set_title(title)
        fig.canvas.draw()
        #fig.canvas.flush_events()
        plt.show(block=False)
        #fig.canvas.draw_idle()
        plt.pause(0.1)
        plt.savefig("/home/g/Asteroids/output/Ishan_2/plots"+"/img_"+str(i)+".svg",dpi=300,bbox_inches="tight")
    #plt.pause(10)
    #plt.close(fig)
    
    
def save_2D_omega(parameters):
    #pdb.set_trace()
    file1=f'Omega_{"C" if parameters["Collision"]=="Yes" else ""}{"L" if parameters["Landslide"]=="Yes" else ""}{"Y" if parameters["YORP"]=="Yes" else ""}.txt'
    file = Output_File(parameters=parameters, filetype='output', filenames=[file1]) 
    omega = np.loadtxt(file, dtype=float)
        
    length, width = np.shape(omega)

    fig, ax = plt.subplots()
    fig.patch.set_facecolor('white')


    ax.set_facecolor('white')
    plt.rcParams.update({'font.size' : 18})
    ax.xaxis.label.set_color('black')
    ax.yaxis.label.set_color('black')
    ax.title.set_color('black')
    ax.tick_params(axis='x', colors='black')
    ax.tick_params(axis='y', colors='black')
    ax.tick_params(direction='in',which='minor', length=5, bottom=True, top=True, left=True, right=True)
    ax.tick_params(direction='in',which='major', length=10, bottom=True, top=True, left=True, right=True)
    ax.yaxis.get_offset_text().set_color('black')   

     
    x1 =  omega[:, 0]/1e+6
    y1 =  omega[:, 1]
    z1 = omega[:, 1]- omega[:,2]
    z2 = omega[:, 1]+ omega[:,2]
    ax.plot(x1, y1, '-r', linewidth=2)
    ax.fill_between(x1,z1,z2,color='#D3D3D3',alpha=0.5)
    ax.spines['top'].set_visible(True)
    ax.spines['right'].set_visible(True)
    ax.spines['left'].set_color('black')
    ax.spines['bottom'].set_color('black')
    ax.spines['right'].set_color('black')
    ax.spines['top'].set_color('black')
    ax.set_xlabel('Time (Myr)',color='black',fontsize=18)
    ax.set_ylabel('Rotational time period  (hr)',color='black',fontsize=18)
    ax.set_xlim([0,2])
    ax.set_ylim([3.98,4.75])
    img = Output_File(parameters=parameters, filetype='output', filenames=['plots','omega.svg']) 
    fig.savefig(img, dpi=150, bbox_inches='tight')  # Save as PNG
    
    
    
    plt.show()
    
#%%

from matplotlib.patches import ConnectionPatch

def save_2D_plot(parameters):
   # pdb.set_trace()
    file1=f'Omega_{"C" if parameters["Collision"]=="Yes" else ""}{"L" if parameters["Landslide"]=="Yes" else ""}{"Y" if parameters["YORP"]=="Yes" else ""}.txt'
    file = Output_File(parameters=parameters, filetype='output', filenames=[file1]) 
    omega = np.loadtxt(file, dtype=float)
        
    length, width = np.shape(omega)

    fig = plt.figure(figsize=(16, 9),dpi=150)
    margin_l =0.07
    margin_b =margin_l*16/9
    box_l =1-0.03 - margin_l
    box_w = 1-0.03 - margin_b
    ax = fig.add_axes([margin_l, margin_b, box_l, box_w])
    
    fig.patch.set_facecolor('white')

    plt.rcParams.update({'font.size' : 18})
    ax.tick_params(direction='in',which='minor', length=5, bottom=True, top=True, left=True, right=True)
    ax.tick_params(direction='in',which='major', length=10, bottom=True, top=True, left=True, right=True)
 

     
    x1 =  omega[:, 0]/1e+6
    y1 =  omega[:, 1]
    z1 = omega[:, 1]- omega[:,2]
    z2 = omega[:, 1]+ omega[:,2]
    ax.plot(x1, y1, '-b', linewidth=2)
    ax.fill_between(x1,z1,z2,color='#B0E0E6',alpha=0.5,zorder=10)

    ax.set_xlabel('Time (Myr)',color='black',fontsize=18)
    ax.set_ylabel('Rotational time period  (hr)',color='black',fontsize=18)
    xlim = [0,2.0]
    ax.set_xlim(xlim)
    ylim = [3.98,5.4]
    ax.set_ylim(ylim)
    
    #plt.tight_layout()
        
    file = Output_File(parameters=parameters, filetype='output', filenames=['mean_shape.csv']) 
    w = np.loadtxt(file, dtype=float, delimiter=",")

    length, width = np.shape(w)

    # Explicitly specify the height and width of each subplot
   # ax1.set_position([0.10, 0.2, 0.35, 0.7])  # [left, bottom, width, height]
   # ax2.set_position([0.52, 0.12, 0.42, 0.84])  # [left, bottom, width, height]

    custom_gray = LinearSegmentedColormap.from_list("my_gray", [ "#F5F5F5", "#7F7F7F"])
    def update_plot(ax,i):
        # First subplot
        
        ax.cla()
        x = np.sin(w[0, 1:width]) * w[i, 1:width]
        y = np.cos(w[0, 1:width]) * w[i, 1:width]
        ax.plot(x, y, '-k', linewidth=2)
        #ax.fill_betweenx(y,x,-x,color='lightcoral',alpha =0.3)
        # 2. Create the "Stencil" (The path to fill)
        # We create a dummy fill and extract its path
        path = ax.fill_betweenx(y, x, -x, color='none').get_paths()[0]
        
        # 3. Create a Gradient Image
        # 'Greys' is a built-in colormap; we create a 2D array for the gradient
        gradient = np.linspace(0, 1, 256).reshape(-1, 1)
        
        # 4. Display the gradient and CLIP it to the path

        img = ax.imshow(gradient, aspect='auto', 
                        extent=[-x.max(), x.max(), -y.max(), y.max()],
                        cmap=custom_gray, origin='lower', alpha=1,zorder=1)
        
        img.set_clip_path(path, transform=ax.transData)
        
        x = -np.sin(w[0, 1:width]) * w[i, 1:width]
        ax.plot(x, y, '-k', linewidth=2)

        # Set x and y limits
        ax.set_xlim(-1.02*np.max(w[:, 1:width]), 1.02*np.max(w[:, 1:width]))
        ax.set_ylim(-1.02*np.max(w[:, 1:width]), 1.02*np.max(w[:, 1:width]))
        ax.set_aspect('equal')
        ax.set_axis_off()

    box_size = 0.4
    axes =[]
    ax1 = [0, 0 , box_size, box_size]
    axes.append(ax1)
    ax1 = [-0.08, 0.6 , box_size, box_size]
    axes.append(ax1)
    ax1 = [0.23, 0 , box_size, box_size]
    axes.append(ax1)
    ax1 = [0.12, 0.6 , box_size, box_size]
    axes.append(ax1)
    ax1 = [0.46, 0 , box_size, box_size]
    axes.append(ax1)
    ax1 = [0.46, 0.35 , box_size, box_size]
    axes.append(ax1)
    ax1 = [0.69, 0. , box_size, box_size]
    axes.append(ax1)
    ax1 = [0.69, 0.35 , box_size, box_size]
    axes.append(ax1)
    print(np.max(w[:, 1:width]))
    for i in range(1,9):

        ax1 = ax.inset_axes(axes[i-1])
        ax1.patch.set_alpha(0)
        ax1.set_zorder(5)
        j= 1+ math.ceil((length-2)*xlim[1]/2.0/8*i)
        update_plot(ax1,j)
        plt.minorticks_on()

# xyA is the point on the inset, xyB is the point on the main plot
        con = ConnectionPatch(xyA=(0.5, 0.5), xyB=(i/8, 0+(i-1)%2), 
                      coordsA=ax1.transAxes, # Relative to Inset (0,0 is bottom-left)
                      coordsB=ax.transAxes,  # Relative to Main Plot Data
                      arrowstyle="->", color="red", lw=1.5)

# 2. Add it to the FIGURE (not the axis) so it stays on top
        fig.add_artist(con)

# 3. Ensure it's in the foreground
        con.set_zorder(100)
        
    #img = Output_File(parameters=parameters, filetype='output', filenames=['plots','omega.svg']) 
    ax.grid(zorder = -1)
    #fig.savefig(img, dpi=150)  # Save as PNG
    plt.show()
    return axes


