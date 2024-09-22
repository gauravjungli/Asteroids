#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jul  6 01:33:23 2024

@author: g
"""
from plots import post_process, show_omega
from gaurav import Parameter

parameters={'run': 1, 'Output folder': 'check' }
Parameter(parameters, filetype="output")

#post_process(parameters)
show_omega(parameters)


