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
phi = np.linspace(0,np.pi*2,4)
x_init = R*np.cos(phi)
y_init = R*np.sin(phi)

vy_init = -const.omega*R*np.sin(phi)
# vy_init = 
vx_init = const.omega*R*np.cos(phi)

# VX_init, VY_init = np.meshgrid(x_init,y_init)

with open('Initial_Conditions.json', 'w') as Init:
    for x,y,vx,vy in zip(x_init,y_init,vx_init,vy_init):
        # for y in y_init:
        # x = xy[0]
        # y = xy[1] 
        # # for vy in vy_init:
        # vx = vxy[0]
        # vy = vxy[1]
        json.dump({"x0": x, "y0": y, "z0": 0, "vx0": vx, "vy0": vy, "vz0": 0, "Sunr0": 1, "Sun": 0, "Shadow": 0, "plotting": 1}, Init)
        Init.write('\n')
