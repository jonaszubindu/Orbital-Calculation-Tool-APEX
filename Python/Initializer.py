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
from matplotlib import pyplot as plt
from classes.methods_2 import cart_to_sph
import pandas as pd
# x_init = np.linspace(-100,100,5)
# y_init = np.linspace(1000,1300,5)
# # y_init = [0]

# vy_init = np.linspace(-0.1,+0.1,3)
# # vy_init = 
# vx_init = -y_init*const.omega

def end_to_init(orbits, Sun, Shadow):
    """
    

    Parameters
    ----------
    orbits : TYPE
        DESCRIPTION.

    Returns
    -------
    kwargs_list : TYPE
        DESCRIPTION.

    """
    kwargs_list = []
    if type(orbits) != list:
        orbits = [orbits]
    for orbit in orbits:
        x = orbit.rAPEX_x[-1]
        y = orbit.rAPEX_y[-1]
        z = orbit.rAPEX_z[-1]
        
        vx = orbit.vAPEX_x[-1]
        vy = orbit.vAPEX_y[-1]
        vz = orbit.vAPEX_z[-1]
        orbit_number = orbit.orbit_number
        kwargs_list.append({"x0": x, "y0": y, "z0": z, "vx0": vx, "vy0": vy, "vz0": vz, "Sunr0": 1, "Sun": Sun, "Shadow": Shadow, "plotting": 1, "orbit_number": orbit_number})
    
    return kwargs_list


def corot_speed(x,y):
    if np.shape(x) == ():
        z = 0
    else:
        z = np.zeros_like(x)
        x = pd.DataFrame(x)
        y = pd.DataFrame(y)
        z = pd.DataFrame(z)
    r, _, phi = cart_to_sph(x,y,z)
    vx = -const.omega*r*np.sin(phi)
    vy = const.omega*r*np.cos(phi)
    try:
        vx = vx.to_numpy()
        vy = vy.to_numpy()
    except Exception:
        pass
    return vx, vy

def create_initials(Sun, Shadow):
    """
    
    
    Parameters
    ----------
    R : TYPE
        DESCRIPTION.
    Sun : TYPE
        DESCRIPTION.
    Shadow : TYPE
        DESCRIPTION.
    
    Returns
    -------
    kwargs_list : TYPE
        DESCRIPTION.
    
    """
    # Sun = Shadow = 0
    kwargs_list = []
    # phi1 = np.linspace(.5,np.pi-.5,20)
    # phi2 = np.linspace(-np.pi+.5,-.5,20)
    # phi3 = np.linspace(-.7,.7,5)
    # phi4 = np.linspace(np.pi-.7,np.pi+.7,20)
    # R = np.linspace(700,1400,20)
    # R1 = np.linspace(600,1000,20)
    # R2 = np.linspace(1200,1500,20)
    # R3 = np.linspace(900,1300,20)
    # R4 = np.linspace(900,1100,20)
    
    # xf = np.linspace(-1600, 1600, 80)
    # yf = np.linspace(-1600, 1600, 80)
    # zf = np.linspace(-1600, 1600, 80)
    
    # XF, YF = np.meshgrid(xf, yf)
    # ZF = np.zeros_like(XF)
    
    # R=1200
    # phi = np.linspace(1,1.6,3)
    
    # x_init = R*np.cos(phi)
    # y_init = R*np.sin(phi)
    # vx_init = -const.omega*R*np.sin(phi)
    # vy_init = const.omega*R*np.cos(phi)
    
    # L4
    # x_init11 = np.outer(R,np.cos(phi1))
    # y_init11 = np.outer(R,np.sin(phi1))
    
    
    # # L5
    # x_init12 = np.outer(R1,np.cos(phi2))
    # y_init12 = np.outer(R1,np.sin(phi2))
    
    
    # # L2
    # x_init23 = np.outer(R2,np.cos(phi3))
    # y_init23 = np.outer(R2,np.sin(phi3))
    
    
    # # L1
    # x_init43 = np.outer(R4,np.cos(phi3))
    # y_init43 = np.outer(R4,np.sin(phi3))
    
    
    # # L3
    # x_init34 = np.outer(R3,np.cos(phi4))
    # y_init34 = np.outer(R3,np.sin(phi4))
    
    
    # x_init = [x_init11, x_init12, x_init23, x_init43, x_init34]  
    # y_init = [y_init11, y_init12, y_init23, y_init43, y_init34] 
    
    
    # vx_init = []
    # vy_init = []
    # for x_in, y_in in zip(x_init, y_init):
    #     vx_in, vy_in = corot_speed(x_in, y_in) 
    #     vx_init.append(vx_in)
    #     vy_init.append(vy_in)
    
    
    
    # if type(x_init) != np.ndarray:
    #     x_init = [x_init]
    #     y_init = [y_init]
    #     vx_init = [vx_init]
    #     vy_init = [vy_init]
    # # VX_init, VY_init = np.meshgrid(x_init,y_init)
    
    # R = np.linspace(375,1500,20)
    # phi = np.linspace(0,2*np.pi,36)
    
    # x_init = np.outer(R,np.cos(phi))
    # y_init = np.outer(R,np.sin(phi))
    
    # x_init = np.linspace(-1500,1500,100)
    # y_init = np.linspace(-1500,1500,100)
    # X_init, Y_init = np.meshgrid(x_init,y_init)
    # vx_init, vy_init = corot_speed(X_init, Y_init)
    
    x0 = 500
    y0 = 1000
    # z0 = 0
    vxc, vyc = corot_speed(x0, y0)
    vx_init = vxc + np.linspace(-1, 1, 100)
    vy_init = vyc + np.linspace(-1, 1, 100)
    
    VX_init, VY_init = np.meshgrid(vx_init, vy_init)
        
    # plt.scatter(x_init, y_init)
    # plt.show()
    n = 0
    with open('Initial_Conditions.json', 'w') as Init:
        for vx,vy in zip(VX_init, VY_init):#,vx_init,vy_init):
            for vx,vy in zip(np.nditer(vx),np.nditer(vy)):
                    
                # x = float(x)
                # y = float(y)
                vx = float(vx)
                vy = float(vy)
                """ Writing all the generated arguments to dicts that get stored as sets of initial conditions in kwargs_list """
                json.dump({"x0": x0, "y0": y0, "z0": 0, "vx0": vx, "vy0": vy, "vz0": 0, "Sunr0": 1, "Sun": Sun, "Shadow": Shadow, "plotting": 1}, Init)
                Init.write('\n')
                kwargs_list.append({"x0": x0, "y0": y0, "z0": 0, "vx0": vx, "vy0": vy, "vz0": 0, "Sunr0": 1, "Sun": Sun, "Shadow": Shadow, "plotting": 1, "orbit_number": n})
                n += 1
    return kwargs_list


