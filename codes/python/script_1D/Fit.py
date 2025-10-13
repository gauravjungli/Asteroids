#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 20 18:16:49 2025

@author: g
"""

from circle_fit import taubinSVD
import numpy as np
import matplotlib.pyplot as plt

from IO import Parameter, Exparameter, Output_File
from scipy.optimize import least_squares
import os

import pdb
from scipy.interpolate import UnivariateSpline, CubicSpline  , interp1d 
from scipy.signal import savgol_filter
from scipy.ndimage import gaussian_filter1d

#%%


"""
    It is used for the best fit circle. Currently we are only using the average of min and max.
"""


def Fit(parameters,epsilon=None,base_old=None):
    

    dia = float(parameters["Current diameter"])
    Res = int(parameters["Resolution"])
    
    if epsilon is None:
        epsilon = float(parameters["epsilon"])
        
    min_epsilon = float(parameters['Minimum epsilon'])
    initial_mass = float(parameters['Initial mass'])
    shed_mass = float(parameters['Mass shed'])
    file = Output_File(parameters,"output",["base.txt"])
    max_curv =5 #change
    max_slope = 1
    
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
    
   # height = gaussian_filter1d(height, sigma=2)

    # plt.plot(theta,height/epsilon,'-r')
    # plt.plot(theta,base_old[2:Res-2,2])
    # plt.draw()
    # plt.pause(0.5)    

    dbase = base_old[2:Res-2,3]

    metric = base**2 + dbase**2

    rad = np.sqrt((2*base**2*np.sqrt(metric)*height + metric*(height**2 + base**2  ))/metric)
   # rad_mean = np.mean(rad)
    print("maximum radial distance is ", max(rad) )
    z = np.cos(theta)*base + 1/np.sqrt(metric)*(height*dbase*np.sin(theta)+np.cos(theta)*base*height)
    theta_new = np.arccos(z/rad)
    print("Theta values lie in ",min(theta_new), " ",max(theta_new))
    sorted_indices = np.argsort(theta_new)
    theta_sorted = theta_new[sorted_indices]
    rad_sorted = rad[sorted_indices]

   # rad_sorted_dev = (rad_sorted - rad_mean)/epsilon

    
    # s = 1e-6
    # while(True):
    
    #      spline_smooth = UnivariateSpline(theta_sorted,rad_sorted_dev, s=s, k=3) # k=3 for cubic spline
    #      spline_deriv2_smooth = spline_smooth.derivative(n=2)
    #      ddbase = spline_deriv2_smooth(theta)*epsilon #change
    #      spline_deriv1_smooth = spline_smooth.derivative(n=1)
    #      dbase = spline_deriv1_smooth(theta)*epsilon
    #      # plt.clf()
    #      # plt.plot(theta,rad)
    #      # plt.plot(theta,spline_smooth(theta)*epsilon + rad_mean )
    #      # plt.draw()
         
    #      print(f'Extreme values of curvature are  {max(ddbase)}  and {min(ddbase)}')
    #      if ((max(ddbase)<max_curv and min(ddbase)>-max_curv) or s>Res):
    #          break
    #      else:
    #          s= 1.1*s
    # print(f"Using s value as {s}")
    
    #spline_smooth = UnivariateSpline(theta_sorted,rad_sorted_dev, s=s, k=3) # k=3 for cubic spline
    f_extrapolate = interp1d(theta_sorted, rad_sorted, kind='linear', fill_value="extrapolate")
    shape =np.zeros((Res,5))
    shape[:,0] = base_old[:,0]
    shape[2:Res-2,1] = f_extrapolate(theta) #gaussian_filter1d(f_extrapolate(theta), sigma=int(2))
    #shape[2:Res-2,1] = savgol_filter(rad, window_length=int(Res/50), polyorder=3)#spline_smooth(theta)*epsilon + rad_mean 
    # plt.plot(theta, shape[2:Res-2,1],'-r')
    # plt.plot(theta,rad)
    # plt.draw()
    # plt.pause(0.5)

    
    print("New maximum radial distance is ", max(shape[2:Res-2,1]) )

    new_mass = 2*np.pi/3*np.trapezoid(((shape[2:Res-2,1])**3*np.sin(theta)),theta)*(dia/2)**3
    
    
    if impact_flag:
        
        new_mass += 2*np.pi*pow(dia/2,3)*np.trapezoid(shape[2:Res-2,1]*np.sin(theta)*np.abs(height)*np.sqrt(metric),theta)
        
    dia_ratio = ((initial_mass-shed_mass)/new_mass)**(1/3)
    
    shape[2:Res-2,1] = (shape[2:Res-2,1])*dia_ratio
    print("The dia ratio is ", dia_ratio)
    shape[0,1]= float(parameters["Diameter"])/dia
    shape[1,1]= shape[0,1]
    shape[Res-1,1] = shape[1,1]
    shape[Res-2,1] = shape[0,1]
    
    # Derivatives
    spline_deriv1_smooth = np.gradient(shape[2:Res-2,1], shape[2:Res-2,0])#spline_smooth.derivative(n=1)
    #shape[2:Res-2,3] =  savgol_filter(spline_deriv1_smooth, window_length=int(Res/50), polyorder=3) #*epsilon
    shape[2:Res-2,3] =  gaussian_filter1d(spline_deriv1_smooth, sigma=int(2))
    spline_deriv2_smooth = np.gradient(shape[2:Res-2,3],shape[2:Res-2,0])#spline_smooth.derivative(n=2)
    
    
    #shape[2:Res-2,4] = savgol_filter(spline_deriv2_smooth, window_length=int(Res/50), polyorder=3)#*epsilon 
    shape[2:Res-2,4] =gaussian_filter1d(spline_deriv2_smooth, sigma=int(2))
    
    #change comment it for a normal run 
    plt.plot(shape[2:Res-2,0],spline_deriv2_smooth,'-r')
    plt.plot(shape[2:Res-2,0],shape[2:Res-2,4])
    plt.draw()
    plt.pause(0.5)
    
    if max(shape[:,4])>max_curv or min(shape[:,4])<-max_curv:
        print("Too large value of curvature",max(shape[2:Res-2,4]),"   ",min(shape[2:Res-2,4]) )
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
    print("The maximum value of first derivative is ", max(shape[:,3]))
    print("The minimum value of first derivative is ", min(shape[:,3]))

    parameters["Current diameter"] = float(radius)*dia

    jinertia =  np.trapezoid(2*np.pi/5*np.sin(shape[1:Res-1,0])**3*shape[1:Res-1,1]**5,shape[1:Res-1,0])
    jinertia1 = np.trapezoid(np.pi/5*(-np.sin(shape[1:Res-1,0])**3+2*np.sin(shape[1:Res-1,0]))*shape[1:Res-1,1]**5,shape[1:Res-1,0])
    print(jinertia,'      ' ,jinertia1)
    parameters["jinertia"] = jinertia
    parameters["jinertia1"] = jinertia1
    Exparameter(parameters)
    if not impact_flag:
        np.savetxt(file,shape,delimiter=",") 
    else:
        return shape
   # print ("The best fit value of r is ", r)
  

