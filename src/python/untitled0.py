#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jul  6 01:33:23 2024

@author: g
"""
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
from matplotlib.colors import LinearSegmentedColormap, Normalize
from matplotlib.cm import ScalarMappable

# Create some data
x = np.linspace(-5, 5, 100)
y = np.linspace(-5, 5, 100)
x, y = np.meshgrid(x, y)
z = np.sin(np.sqrt(x**2 + y**2))

# Calculate the sum of x, y, and z values at each grid point
sum_xyz = x + y + z

# Define a custom colormap
colors = [(0, 0, 0.5), (0, 0.5, 0), (1, 1, 0)]  # Dark blue to green to yellow
n_bins = 100  # Number of bins in the colormap

# Create the colormap
cmap_name = 'my_custom_cmap'
cm = LinearSegmentedColormap.from_list(cmap_name, colors, N=n_bins)

# Normalize the sum values to 0-1 range
norm = Normalize(vmin=sum_xyz.min(), vmax=sum_xyz.max())
normalized_values = norm(sum_xyz)

# Map the normalized values to colors
facecolors = cm(normalized_values)

# Create a figure with a specific size
fig = plt.figure(figsize=(10, 8))

# Add a 3D subplot with specific location and size
ax = fig.add_axes([0.1, 0.1, 0.75, 0.75], projection='3d')  # [left, bottom, width, height]

# Plot the surface with the custom colormap
surf = ax.plot_surface(x, y, z, facecolors=facecolors, rstride=1, cstride=1, linewidth=0, antialiased=False)

# Add a color bar
sm = ScalarMappable(cmap=cm, norm=norm)
sm.set_array(sum_xyz)
cbar_ax = fig.add_axes([0.86, 0.1, 0.03, 0.75])  # [left, bottom, width, height]
fig.colorbar(sm, cax=cbar_ax)

# Show the plot
plt.show()



