#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 20 18:16:49 2025

@author: g
"""
from IO import Output_File
from circle_fit import taubinSVD
import numpy as np
import matplotlib.pyplot as plt
import pdb
from scipy.interpolate import  interp1d 
from scipy.ndimage import gaussian_filter1d
import os
#%%


"""
    It is used for the best fit circle. Currently we are only using the average of min and max.
"""
def Fit_radius(base_old,epsilon):
    
    Res = base_old.shape[0]
    
    theta = base_old[2:Res-2,0]
    base  = base_old[2:Res-2,1]    
    height = epsilon*base_old[2:Res-2,2]
    
    height = gaussian_filter1d(height, sigma=int(len(theta)/200))

    dbase = base_old[2:Res-2,3]

    metric = base**2 + dbase**2

    rad = np.sqrt(( 2*base**2*np.sqrt(metric)*height + metric*(height**2 + base**2 ) )/metric)
   # rad_mean = np.mean(rad)
    #print("maximum radial distance is ", max(rad) )
    z = np.cos(theta)*base + 1/np.sqrt(metric)*(height*dbase*np.sin(theta)+np.cos(theta)*base*height)
    theta_new = np.arccos(z/rad)
    #print("Theta values lie in ",min(theta_new), " ",max(theta_new))
    sorted_indices = np.argsort(theta_new)
    theta_sorted = theta_new[sorted_indices]
    rad_sorted = rad[sorted_indices]

    f_extrapolate = interp1d(theta_sorted, rad_sorted, kind='linear', fill_value="extrapolate")

    return f_extrapolate(base_old[:,0])



def Fit(target,epsilon=None,base_old=None):
    
    dia = target.d
    Res =target.res
    
    epsilon = target.epsilon 
        
    initial_mass = target.initial_mass
    shed_mass = target.shed_mass
    file = Output_File(target,"output",["base.txt"])
    max_curv =50
    max_slope = 10
    
    #Flag to check whether it is called after an impact or just to fit the shape
    impact_flag = False
    if base_old is not None:
        impact_flag = True
    if not impact_flag:
        try:
            base_old = np.loadtxt(file,dtype=float,delimiter=",")
        except:
            print("cannot import shape")
            return
        
    theta = base_old[2:Res-2,0]
    base  = base_old[2:Res-2,1]    
    height = epsilon*base_old[2:Res-2,2]
    jinertia =  np.trapezoid(2*np.pi/5*np.sin(base_old[1:Res-1,0])**3*base_old[1:Res-1,1]**5,base_old[1:Res-1,0])
    
    if not impact_flag:
        
        fig, ax = plt.subplots()
        ax.plot(theta,height) 
    
    height = gaussian_filter1d(height, sigma=int(Res/200))
    
    if not impact_flag:
        
        ax.plot(theta,height)
        save_folder =  Output_File(target,"output",["height_profile"])
        os.makedirs(save_folder, exist_ok=True) 
        fig.savefig(save_folder + f"/height_{target.slides}.png")
        plt.close(fig)
    
    plot = False
    if plot:
        plt.plot(theta,height/epsilon,'-r') 
        plt.plot(theta,base_old[2:Res-2,2])
        plt.show()
        plt.pause(0.5)    
    
    dbase = base_old[2:Res-2,3]

    metric = base**2 + dbase**2
    
    shape =np.zeros((Res,5))
    shape[:,0] = base_old[:,0]
    shape[2:Res-2,1] = Fit_radius(base_old,epsilon)[2:Res-2]
    
    #print("New maximum radial distance is ", max(shape[2:Res-2,1]) )

    new_mass = 2*np.pi/3*np.trapezoid(((shape[2:Res-2,1])**3*np.sin(theta)),theta)*(dia/2)**3
    
    
    if impact_flag:
        
        new_mass += 2*np.pi*pow(dia/2,3)*np.trapezoid(shape[2:Res-2,1]*np.sin(theta)*np.abs(height)*np.sqrt(metric),theta)
        
    dia_ratio = ((initial_mass-shed_mass)/new_mass)**(1/3)
    
    shape[2:Res-2,1] = (shape[2:Res-2,1])*dia_ratio
    print("The dia ratio is ", dia_ratio)
    shape[0,1]= shape[2,1]
    shape[1,1]= shape[0,1]
    shape[Res-1,1] = shape[1,1]
    shape[Res-2,1] = shape[0,1]
    
    # Derivatives
    spline_deriv1_smooth =gaussian_filter1d(shape[2:Res-2,1], sigma =int(Res/100) )
    
    shape[2:Res-2,3] =  np.gradient(spline_deriv1_smooth, shape[2:Res-2,0]) 
   
    
    if not impact_flag:
        fig, ax = plt.subplots()
        ax.plot(shape[2:Res-2,0],np.gradient(shape[2:Res-2,1], shape[2:Res-2,0]),'-r')
        save_folder =  Output_File(target,"output",["first_derivative"])
        os.makedirs(save_folder, exist_ok=True) 
        ax.plot(shape[2:Res-2,0],shape[2:Res-2,3])
        fig.savefig(save_folder + f"/derivative_{target.slides}.png")
        plt.close(fig)

   
    spline_deriv2_smooth = gaussian_filter1d(shape[2:Res-2,3], sigma=int(Res/100))

    shape[2:Res-2,4] = np.gradient(spline_deriv2_smooth,shape[2:Res-2,0])


    if not impact_flag:
        fig, ax = plt.subplots()
        ax.plot(shape[2:Res-2,0],np.gradient(np.gradient(shape[2:Res-2,1], shape[2:Res-2,0]),shape[2:Res-2,0]),'-r')
        save_folder =  Output_File(target,"output",["Second_derivative"])
        os.makedirs(save_folder, exist_ok=True) 
        ax.plot(shape[2:Res-2,0],shape[2:Res-2,4])
        fig.savefig(save_folder + f"/second_derivative_{target.slides}.png")
        plt.close(fig)
        
    
    if max(shape[:,4])>max_curv or min(shape[:,4])<-max_curv:
        print("Too large value of curvature",max(shape[2:Res-2,4]), "   ", min(shape[2:Res-2,4]) )
    shape[:,3] = np.clip(shape[:,3],-max_slope,max_slope)
    shape[:,4] = np.clip(shape[:,4],-max_curv,max_curv)
    
    
    x=np.sin(shape[:,0])*(shape[:,1])
    y=np.cos(shape[:,0])*(shape[:,1])

# fig = plt.figure(figsize=(6,6))
# plt.axis('equal')
# plt.plot(x,y)


    point = []
    for i in range(Res):
        point.append([x[i],y[i]])
        point.append([-x[i],y[i]])
    xc, yc, radius, sigma = taubinSVD(point)
    
    
    print("The maximum value of second derivative is ", max(shape[:,4]))
    print("The minimum value of second derivative is ", min(shape[:,4]))


    #Updating the diameter. It is only place where the diameter is updated
    target.d = float(radius)*dia

    jinertia =  np.trapezoid(2*np.pi/5*np.sin(shape[1:Res-1,0])**3*shape[1:Res-1,1]**5,shape[1:Res-1,0])
    jinertia1 = np.trapezoid(np.pi/5*(-np.sin(shape[1:Res-1,0])**3+2*np.sin(shape[1:Res-1,0]))*shape[1:Res-1,1]**5,shape[1:Res-1,0])
   # print(jinertia,'      ' ,jinertia1)
    target.jinertia[2] = jinertia*(target.d/2)**5 * target.dens
    target.jinertia[0] = target.jinertia[0]  = jinertia1*(target.d/2)**5 * target.dens

    if not impact_flag:
        np.savetxt(file,shape,delimiter=",") 
    else:
        return shape
   
  

#To be completed later
def Fit_2D(target):
    
    file = Output_File(target,"output",["base.txt"])

    #Flag to check whether it is called after an impact or just to fit the shape
   
    
    base = np.loadtxt(file,dtype=float,delimiter=",")
    
    base[:,2] =base[:,2] + base[:,3]
    base[:,3] =0
    np.savetxt(file,base,delimiter=",")
    
    #pdb.set_trace()
    # if parameters['Dimension'] == '2D':
    #   #  pdb.set_trace()
    #     nx, ny = int(parameters['X Resolution']), int(parameters['Y Resolution'])
    #     epsilon = float(parameters['epsilon'])
    #     radius =target.d/2
    #     reshaped_base = base.reshape(nx, ny, 7)

    # # 2. Unpack the last axis into three separate (N, M) arrays
    #     x_vals = reshaped_base[:, :, 0]
    #     y_vals = reshaped_base[:, :, 1]
    #     b_vals = reshaped_base[:, :, 2]
    #     h_vals = reshaped_base[:, :, 3]
    #     r_vals = reshaped_base[:, :, 4]
    #     dr_vals = reshaped_base[:, :, 5]
    #     ddr_vals = reshaped_base[:, :, 6]
            


    #     devs_new = (b_vals +h_vals)*epsilon
    #     devs_new_axi =  np.mean(devs_new,axis=1)
    #     devs_new_1 = devs_new -devs_new_axi.reshape(-1,1)
    #     r_vals = r_vals + -devs_new_axi.reshape(-1,1)
    
    #     # Calculate first derivative: dr/dtheta
    #     r_smooth = gaussian_filter1d(r_vals, sigma=2)
    #     dr = np.gradient(r_smooth, x_vals, edge_order=2)
        
    #     dr_smooth = gaussian_filter1d(dr, sigma=2)
    #     # Calculate second derivative: d^2r/dtheta^2
    #     d2r = np.gradient(dr_smooth, x_values, edge_order=2)
    #     print(np.max(d2r),np.min(d2r))
    #     #to accomodate for the fact that the normal to axisymmtric surface and sphere is not aligned
    #     cos_theta = r_vals/np.sqrt(r_vals**2+dr**2)
    #     sin_theta = dr/np.sqrt(r_vals**2+dr**2)
    #     d_surface = np.gradient(devs_new, axis=0)/(devs_new + radius)
    #     cos_theta_matrix =  np.tile(cos_theta.reshape(-1,1),(1,ny))
    #     sin_theta_matrix =  np.tile(sin_theta.reshape(-1,1),(1,ny))
    #     devs_new_2 = devs_new_1/(cos_theta_matrix + d_surface*sin_theta_matrix)
        
        
    #     target.d = 2 * radius*1000
    #     target.rgrav,target.tgrav = Gravitycalc_2D(parameters,grid1_lats,grid1_lons,devs_axi,radius,x_values,y_values,devs_new_axi,r_vals,dr)
    #     devs_new_non = devs_new_2.ravel()/radius
        
    #     return x_values, y_values, r_vals, dr, d2r, devs_new_non