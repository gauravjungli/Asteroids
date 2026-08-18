#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 18 11:13:39 2023

@author: kumargaurav
"""

import math
import numpy as np
from collisions import G, qstarf
import time
from YORP import read_f_g_spline, shape_gen
from Diffusion_spherical import find_roots_parallel, parallel_root_computation, compute_energy



#%%


class Target:
    
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

    def __init__(self, parameters):
        
        self.name = parameters['Asteroid']
        self.d = float(parameters["Diameter"])
        self.atype = parameters["atype"]
        self.delta = float(parameters['Static Friction angle'])
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
            
        self.M = (math.pi / 6) * self.dens * self.d**3
        self.jinertia = [2/5 * self.M * (self.d/2)**2]*3
        
        self.omega = [0,0,2*np.pi/(float(parameters['Rotation period'])*3600)]
        self.kvg = 0.3
        self.crater_coeff = 0.62
        self.K = float(parameters["K"])
        self.grav = G*self.M/(self.d/2)**2
        self.obliq = float(parameters["Obliquity"])
        velave = float(parameters["Impactor velocity"])
        self.dstarave = qstarf(self, math.pi / 4, velave)[2]
        # Set the seed for NumPy's random number generator
        np.random.seed(int(time.time()/float(parameters['run'])))
        self.coeff_f,self.coeff_g = (1,1)#shape_gen(self.K) #change make it (1,1) if removing stochasticity
        self.sma = float(parameters['Semi major axis'])
        self.f_spline, self.g_spline = read_f_g_spline(parameters)
        self.K_0 = float(parameters["Bulk Modulus"])
        self.mu_0 = float(parameters["Shear Modulus"])
        self.wave_speed_P = np.sqrt(self.K_0+4/3*self.mu_0)*(((self.d/2)**2/15/self.dens)*(4*np.pi*G*self.dens-2*self.omega[2]**2))**(1/4)
        self.wave_speed_S = np.sqrt(self.mu_0)*(((self.d/2)**2/15/self.dens)*(4*np.pi*G*self.dens-2*self.omega[2]**2))**(1/4)
        self.failure_mode = parameters['Failure wave']
        if self.failure_mode == "P-wave":
            self.wave_speed =self.wave_speed_P
        else:
            self.wave_speed = self.wave_speed_S 
        self.k_s = float(parameters["Seismic Diffusivity"])
        self.efficiency = float(parameters["Seismic efficiency"])
        self.beta = float(parameters["Beta"])
        self.f = ((2*self.efficiency)/(np.pi*self.beta**2))**(1/3)*2*self.wave_speed_P/0.1 #change_P, we need to fix on this
        self.Q = float(parameters["Q"])
        self.k_d = 2*np.pi*self.f/self.Q
        self.N = 50 
        self.roots = parallel_root_computation(300,50,self.d) 
        self.theta = np.linspace(np.pi/12, np.pi,self.N) 
        self.energy, self.t_max = compute_energy(self)
        self.cohesion_cons = float(parameters["Cohesion constant"])
        self.cohesion_linear = float(parameters["Cohesion linear"])
        self.omegaLimit = (G*4/3 * math.pi*self.dens)**0.5
        self.rgrav = None 
        self.tgrav = None
        self.t_lan = None
        self.fast_rotation_flag = False
        self.epsilon = float(parameters["epsilon"])
        self.res = int(parameters["Resolution"])
        self.x_res = int(parameters['X Resolution'])
        self.y_res = int(parameters['Y Resolution'])
        self.initial_mass = None
        self.offset = float(parameters["offset"])

        self.dx = None
        self.dy = None
        self.dim = parameters['Dimension']
        self.shed_mass = 0
        self.min_epsilon = float(parameters["Minimum epsilon"])
        self.max_epsilon = float(parameters["Maximum epsilon"])
        self.folder = parameters['Output folder']
        self.number = int(parameters['run'])
        self.failure = parameters["Failure profile"]
        self.slides = 0
        
    
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
      



    
