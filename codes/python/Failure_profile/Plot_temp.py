#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul  3 17:39:42 2026

@author: g
"""

import numpy as np
import matplotlib.pyplot as plt

# Define base parameters
k = 1.5
x = np.linspace(0, 3, 400)
y_kx = k * x

# Define the four different second plots
# 1. y = const.
y1 = np.full_like(x, 1.0)

# 2. y = c_0 + c_1*x with c_1 < k (c_0 = 1.5, c_1 = 0.5)
c0_2, c1_2 = 1.0, 0.5
y2 = c0_2 + c1_2 * x

# 3. y = c_0 + c_1*x with c_1 > k (c_0 = -1.0, c_1 = 2.5)
c0_3, c1_3 = 1, 2.0
y3 = c0_3 + c1_3 * x

# 4. y = a*(exp(b*x) - 1) (a = 0.5, b = 0.9)
a, b = 0.5, 0.9
y4 = a * (np.exp(b * x) - 1)

# List of subplots data for iteration
plots_data = [
    (y1, r"$c = 2.0$"),
    (y2, r"$c = 1.0 + 0.5h$"),
    (y3, r"$c = 1.0 + 2.0h$"),
    (y4, r"$c = 0.5(e^{0.9h} - 1)$")
]

# Initialize a 2x2 subplot figure
fig, axs = plt.subplots(2, 2, figsize=(14, 11))
axs = axs.ravel()  # Flatten the 2D array for easy 1D indexing
plt.rcParams.update({'font.size' : 20})
for i, (y_second, title) in enumerate(plots_data):
    ax = axs[i]
    
    # Plot the lines
    ax.plot(x, y_kx, color='blue', linewidth=2.5)
    ax.plot(x, y_second, label=title, color='black', linewidth=2.5,)
    
    # Set plot boundaries
    ax.set_xlim(0, 3)
    ax.set_ylim(0, 5)
    
    # Fill background based on condition:
    # Light green where kx < second plot, Light red where kx >= second plot
    ax.fill_between(x, -1, 5, where=(y_kx < y_second), facecolor='lightgreen', alpha=0.35)
    ax.fill_between(x, -1, 5, where=(y_kx >= y_second), facecolor='lightcoral', alpha=0.35)
    
    # Create masks to find the midpoint of the regions for text placement
    green_mask = y_kx < y_second
    red_mask = y_kx >= y_second
    
    # Add text label for the Stable Region (Green)
    if np.any(green_mask):
        x_green_center = np.mean(x[green_mask])
        ax.text(x_green_center, 2.5, "stable \n region", color='darkgreen', 
                fontsize=18, fontweight='bold', ha='center', va='center',
                bbox=dict(boxstyle='round,pad=0.3', facecolor='white', alpha=0.8, edgecolor='none'))
        
    # Add text label for the Failed Region (Red)
    if np.any(red_mask):
        x_red_center = np.mean(x[red_mask])
        ax.text(x_red_center +0.25, 2.5, "failed \n region", color='darkred', 
                fontsize=18, fontweight='bold', ha='center', va='center',
                bbox=dict(boxstyle='round,pad=0.3', facecolor='white', alpha=0.8, edgecolor='none'))
        
    # Formatting details
    #ax.set_title(title, fontsize=13, pad=10)
    ax.set_xlabel('$h$', fontsize=24)
    ax.set_ylabel('$c$', fontsize=24)
    ax.legend(loc='upper left', framealpha=0.9)
    ax.grid(True, linestyle=':', alpha=0.6)
    ax.tick_params(axis='both', which='major', labelsize=22)
# Adjust spacing and save the figure
plt.tight_layout()
plt.savefig('/home/g/Thesis/subplots_regions.svg', dpi=300)