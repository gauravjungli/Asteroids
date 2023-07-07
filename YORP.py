#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 18 12:09:42 2023

@author: kumargaurav
"""

import numpy as np
import math
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick
from pathlib import Path
import glob

data1=np.loadtxt("/Users/kumargaurav/Documents/OrbFit/tests/gaurav/clo0.yorp",delimiter=" ")
#data2=np.loadtxt("/Users/kumargaurav/Documents/OrbFit/tests/gaurav/clo1.yorp",delimiter=" ")

plt.plot(data1[:,0],data1[:,1])
#plt.plot(data2[:,0],data2[:,1])

#plt.plot(data1[:,0],data1[:,2])
#plt.plot(data2[:,0],data2[:,2])