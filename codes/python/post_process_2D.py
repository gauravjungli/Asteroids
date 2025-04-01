#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 28 10:19:21 2025

@author: g
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import glob
import time
import re
import matplotlib.animation as animation

# Load CSV file
def load_data(file_path):
    data = pd.read_csv(file_path, header=None)
    theta, phi = data[0], data[1]  # Convert to radians
    base_height, flow_height = data[2], data[3]
    vel_x, vel_y = data[4], data[5]
    return theta, phi, base_height, flow_height, vel_x, vel_y

def get_color_limits_vel(folder_path):
    file_paths = sorted(glob.glob(f"{folder_path}/field_*.csv"))
    min_val, max_val = float('inf'), float('-inf')
    for file_path in file_paths:
        _, _, _, flow_height, vel_x, vel_y = load_data(file_path)
        min_val = min(min_val, vel_x.min())
        max_val = max(max_val, vel_x.max())
    return min_val, max_val


def get_color_limits_height(folder_path):
    file_paths = sorted(glob.glob(f"{folder_path}/field_*.csv"))
    min_val, max_val = float('inf'), float('-inf')
    for file_path in file_paths:
        _, _, _, flow_height, vel_x, _ = load_data(file_path)
        min_val = min(min_val, flow_height.min())
        max_val = max(max_val, flow_height.max())
    return min_val, max_val

# Plot height as contour
def plot_height(theta, phi, flow_height, time_step,vmin,vmax):
    plt.rcParams.update({'font.size' : 14})
    sc = plt.tricontourf(phi*180/np.pi, theta*180/np.pi, flow_height, cmap='viridis')#,vmin=vmin,vmax=vmax)
    plt.colorbar(sc, label='Height')
    plt.xlabel('Phi (degree)')
    plt.ylabel('Theta (degree)')
    #plt.title(f'Flow Height Contour (Time Step {time_step})')
    plt.show()
    plt.pause(0.5)
    #plt.axis('equal')
    plt.minorticks_on()
    plt.tick_params(direction='in',right=True, top=True, left=True, bottom=True)
    plt.tick_params(labelsize=14)
    plt.tick_params(labelbottom=True, labeltop=False, labelright=False, labelleft=True)
    #plt.tick_params(direction='in',which='minor', length=5, bottom=True, top=True, left=True, right=True)
    #plt.tick_params(direction='in',which='major', length=10, bottom=True, top=True, left=True, right=True)
    #plt.legend(loc='best')
    #plt.legend(fontsize=14) 
    plt.savefig(f'/home/g/Asteroids/plots_2D/time_evolution/{time_step}.svg',dpi = 300)

# Plot velocity as quiver
def plot_velocity(theta, phi, vel_x, vel_y, time_step):
    plt.figure(figsize=(8, 6))
    plt.quiver(theta, phi, vel_x, vel_y, scale=10, color='red')
    plt.xlabel('Theta (radians)')
    plt.ylabel('Phi (radians)')
    plt.title(f'Velocity Field (Time Step {time_step})')
    plt.show()
    
# Initialize figure
def init_plot():
    fig = plt.figure(figsize=(8, 6))
    ax = fig.add_axes([0.1, 0.1, 0.70, 0.8])  # Manually set plot area (left, bottom, width, height)
    
    cbar_ax = fig.add_axes([0.85, 0.15, 0.03, 0.7])  # Colorbar axis
    return fig, ax, cbar_ax

def update(frame, file_paths, ax, cbar_ax, vmin, vmax):
    ax.clear()
    file_path = file_paths[frame]
    theta, phi, _, flow_height, vel_x, vel_y = load_data(file_path)
    sc = ax.tricontourf(theta*180/np.pi, phi*180/np.pi, flow_height, cmap='viridis', vmin=vmin, vmax=vmax)
    #ax.quiver(theta, phi, vel_x, vel_y, scale=10, color='red')
    ax.set_xlabel('Theta (degree)')
    ax.set_ylabel('Phi (degree)')
    ax.set_title(f'Flow Visualization (Time Step {frame})')
    cbar_ax.clear()
    plt.colorbar(sc, cax=cbar_ax, label='Height')
    return sc

# Main function to execute visualization and save animation
def main(folder_path, output_file="animation.mp4"):
    def extract_number(text):
        match = re.search(r'\d+', text)  # Find first occurrence of a number
        return int(match.group()) if match else float('inf')  # Default to large number if no match
    
    vmin, vmax = get_color_limits_height(folder_path)  # Get global color limits
    file_paths = sorted(glob.glob(f"{folder_path}/field_*.csv"),key=extract_number)
    fig, ax, cbar_ax = init_plot()
    ani = animation.FuncAnimation(fig, update, frames=len(file_paths), fargs=(file_paths, ax, cbar_ax, vmin, vmax), interval=500)
    ani.save(output_file, writer='ffmpeg', fps=2)
    plt.show()


# Main function to execute visualization
def plot_height_main(folder_path,epsilon,Gamma):
    
    def extract_number(text):
        #match = re.search(r'\d+', text)  # Find first occurrence of a number
        #return int(match.group()) if match else float('inf')  # Default to large number if no match
        numbers = [int(num) for num in re.findall(r'\d+', text)]
        return numbers  # Sorting will compare these tuples
    
    
    vmin, vmax = get_color_limits_height(folder_path)  # Get global color limits
    file_paths = sorted(glob.glob(f"{folder_path}/field_*.csv"),key=extract_number)
   # print(file_paths)
    plt.figure(figsize=(8, 6))
    plt.ion()  # Turn on interactive mode
    for i, file_path in enumerate(file_paths):
        print(f"Processing: {file_path}")
        theta, phi, base_height, flow_height, vel_x, vel_y = load_data(file_path)
        
        plt.clf()

        plot_height(theta, phi, (flow_height*epsilon)*250, i,vmin,vmax)
        #plot_velocity(theta, phi, vel_x, vel_y, i)
        time.sleep(1)  # Add delay to observe each step
        print(vmin,vmax)


# Main function to execute visualization
def plot_surface_main(folder_path):
    
    def extract_number(text):
        #match = re.search(r'\d+', text)  # Find first occurrence of a number
        #return int(match.group()) if match else float('inf')  # Default to large number if no match
        numbers = [int(num) for num in re.findall(r'\d+', text)]
        return numbers  # Sorting will compare these tuples
    
    file1 = "/home/g/Asteroids/output/Crater/run1/data/dia.txt"
    epsilon=np.loadtxt(file1,dtype=float)[:,2]
    Gamma = np.loadtxt(file1,dtype=float)[:,3]
    vmin, vmax = get_color_limits_height(folder_path)  # Get global color limits
    file_paths = sorted(glob.glob(f"{folder_path}/field_*.csv"),key=extract_number)
   # print(file_paths)
    plt.figure(figsize=(8, 6))
    plt.ion()  # Turn on interactive mode
    for i, file_path in enumerate(file_paths):
        print(f"Processing: {file_path}")
        theta, phi, base_height, flow_height, vel_x, vel_y = load_data(file_path)
        
        plt.clf()

        plot_height(theta, phi, (Gamma[i]*base_height+epsilon[i]*flow_height[i])*250, i,vmin,vmax)
        #plot_velocity(theta, phi, vel_x, vel_y, i)
        time.sleep(1)  # Add delay to observe each step
        print(vmin,vmax)

        
def plot_vel_main(folder_path):
    
    def extract_number(text):
        match = re.search(r'\d+', text)  # Find first occurrence of a number
        return int(match.group()) if match else float('inf')  # Default to large number if no match
    
    vmin, vmax = get_color_limits_vel(folder_path)  # Get global color limits
    file_paths = sorted(glob.glob(f"{folder_path}/field_*.csv"),key=extract_number)
    plt.figure(figsize=(8, 6))
    plt.ion()  # Turn on interactive mode
    for i, file_path in enumerate(file_paths):
        print(f"Processing: {file_path}")
        theta, phi, base_height, flow_height, vel_x, vel_y = load_data(file_path)
        
        plt.clf()
        plot_height(theta, phi, vel_x, i,vmin,vmax)

        #plot_velocity(theta, phi, vel_x, vel_y, i)
        time.sleep(1)  # Add delay to observe each step

        
        
# Example usage
file1 = "/home/g/Asteroids/output/Crater/run1/data/dia.txt"
epsilon=np.loadtxt(file1,dtype=float)[:,2]
Gamma = np.loadtxt(file1,dtype=float)[:,3]
for i in range(1,2):
  #  plot_height_main(f"/home/g/Asteroids/output/Crater/run1/data/landslides_{i}",epsilon[i],Gamma[i])  # Uncomment and replace with your folder path
    plt.close()

plot_surface_main("/home/g/Asteroids/output/crater_3/run1/data")  # Uncomment and replace with your folder path
