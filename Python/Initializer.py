#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 27 15:37:17 2020

This file is to create a set of initial conditions and write them as json elements to the 
file Initial_Conditions.json  


@author: jonaszbinden
"""
import numpy as np
import json
from classes import constants as const

# x_init = np.linspace(-100,100,5)
# y_init = np.linspace(1000,1300,5)
# # y_init = [0]

# vy_init = np.linspace(-0.1,+0.1,3)
# # vy_init = 
# vx_init = -y_init*const.omega

R = 1200
phi = np.linspace(0,np.pi*2,10)
x_init = R*np.cos(phi)
y_init = R*np.sin(phi)

vy_init = np.linspace(-0.1,+0.1,3)
# vy_init = 
vx_init = -y_init*const.omega

# VX_init, VY_init = np.meshgrid(x_init,y_init)

with open('Initial_Conditions.json', 'w') as Init:
    for xy in (x_init,y_init):
        # for y in y_init:
        x = xy[0]
        y = xy[1] 
        for vx in vx_init:
            for vy in vy_init:
                json.dump({"x0": x, "y0": y, "z0": 0, "vx0": vx, "vy0": vy, "vz0": 0, "Sunr0": 1, "Sun": 0, "Shadow": 0, "plotting": 1}, Init)
                Init.write('\n')
