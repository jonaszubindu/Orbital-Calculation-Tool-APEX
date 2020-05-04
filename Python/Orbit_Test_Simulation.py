#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb  5 10:30:24 2020

This is a full physics orbital simulation software for a double comet system with a small probe 
The simulation checks automatically if there are data available to use for the Kepler Orbit of the 
secondary body. If you change something in the file 'constants.py' to try to investigate different scenarios,
you must compute the new Kepler orbit. This can be achieved by simply deleting the 'Kepler.json' file and
reload the program.


@author: jonaszbinden
"""
#import matplotlib.pyplot as plt
#from matplotlib import animation
# import os

import time
# import sys
from scipy.integrate import ode
from matplotlib import pyplot as plt
# script_dir = "~/Desktop/Studium/Uppsala/Space\ Mission\ Project\ Work/python/"
# sys.path.append(os.path.abspath(script_dir))

import classes
import timeit
from classes import methods_2 as meth
from classes import variables as var
from classes import constants as const
# from matplotlib.animation import FuncAnimation
# from mpl_toolkits.mplot3d import Axes3D

# from mpl_toolkits.mplot3d.axes3d import Axes3D
# from mpl_toolkits.mplot3d import proj3d

# import matplotlib as mpl
import numpy as np

import logging

from matplotlib import patches, ticker
import json
# import concurrent.futures
logging.basicConfig(filename='Orbit_Test_Simulation.log', level= logging.DEBUG)
kwargs_list = []

with open('Initial_Conditions.json', 'r') as inits: #Writing initial conditions from JSON file to a list
    for line in inits:
        kwargs_list.append(json.loads(line))
        
# with open('Stable_Orbits.json', 'r') as inits: #Writing initial conditions from JSON file to a list
#     for line in inits:
#         kwargs_list.append(json.loads(line))



"Main od"
def Orbital_Calculation(initials, variables, kwarg):
        
        v = ode(meth.fun).set_integrator('dop853') #fun is the righthandside of the differential equation, defined in class_ods
        vy0 = initials.vy0
        t0 = initials.t0
        v.set_initial_value(vy0, t0)
        v.set_f_params(variables)
        print('ODE solver successfully initialized...')
        time.sleep(1)
        unstable = 0
        percentage = 0
        "TIME LOOPs"
    
        
        "Integration of differential equation" 
        while variables.tn<(const.T_max-1) and v.successful():
            
            ######################################
            #                                    #
            #   Define breaking arguments here   #
            #                                    #
            ######################################
            
            v.get_return_code()
            variables.tn += 1
            v.integrate(v.t+const.dt)
            
            if variables.tn % int(const.T_max/10) == 0 and variables.tn > 0:
                percentage += 1
                print('Progress: {} %'.format(percentage*10))
            
            variables.vAPEX_x[variables.tn] = v.y[0]
            variables.vAPEX_y[variables.tn] = v.y[1]
            variables.vAPEX_z[variables.tn] = v.y[2]
            variables.rAPEX_x[variables.tn] = v.y[3]
            variables.rAPEX_y[variables.tn] = v.y[4]
            variables.rAPEX_z[variables.tn] = v.y[5]
            
            variables.r_Didymoon_x[variables.tn] = variables.r_Didymoon_i[0][variables.tn]
            variables.r_Didymoon_y[variables.tn] = variables.r_Didymoon_i[1][variables.tn]
            variables.r_Didymoon_z[variables.tn] = variables.r_Didymoon_i[2][variables.tn]
            
            r = np.sqrt(variables.rAPEX_x[variables.tn]**2+variables.rAPEX_y[variables.tn]**2+variables.rAPEX_z[variables.tn]**2)
            r_sec = 0
            r_sec_i = np.zeros(3, dtype = np.float128)
            r_sec_i[0] = variables.rAPEX_x[variables.tn]-variables.r_Didymoon_x[variables.tn]
            r_sec_i[1] = variables.rAPEX_y[variables.tn]-variables.r_Didymoon_y[variables.tn]
            r_sec_i[2] = variables.rAPEX_z[variables.tn]-variables.r_Didymoon_z[variables.tn]

            
            r_sec = np.sqrt(r_sec_i[0]**2+r_sec_i[1]**2+r_sec_i[2]**2)
            variables.t_real[variables.tn] = v.t+const.dt
            variables.rAPEX_r_corot[variables.tn], variables.rAPEX_theta_corot[variables.tn], variables.rAPEX_phi_corot[variables.tn] = meth.corot_frame(variables.rAPEX_x[variables.tn], variables.rAPEX_y[variables.tn], variables.rAPEX_z[variables.tn], variables.r_Didymoon_phi[variables.tn]) 
            variables.rAPEX_x_corot[variables.tn], variables.rAPEX_y_corot[variables.tn], variables.rAPEX_z_corot[variables.tn] = meth.sph_to_cart_coord(variables.rAPEX_r_corot[variables.tn], variables.rAPEX_theta_corot[variables.tn], variables.rAPEX_phi_corot[variables.tn])
            
            # "Lagrange points stable orbits criteria"
            # r_L_x = variables.rAPEX_x_corot[0]-variables.rAPEX_x_corot[variables.tn]
            # r_L_y = variables.rAPEX_y_corot[0]-variables.rAPEX_y_corot[variables.tn]
            # r_L_z = variables.rAPEX_z[0]-variables.rAPEX_z[variables.tn]
            # r_L = np.sqrt(r_L_x**2 + r_L_y**2 + r_L_z**2) 
            
            if r>10*const.a: #or r_L>1000:
                print("APEX left the Lagrange point, unstable", variables.t_real[variables.tn])
                unstable = 1
                break
            if r < const.RadDm:
                print("APEX crashed on Didymos.", variables.t_real[variables.tn])
                unstable = 1
                break
            if r_sec < const.RadDM: 
                print("APEX crashed on Didymoon", variables.t_real[variables.tn])
                unstable = 1
                break
            
        rAPEX_x = variables.rAPEX_x
        rAPEX_y = variables.rAPEX_y
        rAPEX_z = variables.rAPEX_z
        r_Didymoon_x = variables.r_Didymoon_x
        r_Didymoon_y = variables.r_Didymoon_y
        r_Didymoon_z = variables.r_Didymoon_z
        Tn = variables.tn
        print("Test successfully completed, trying next orbit...")
        print('successfully executed run for ')
        str_elem = ""
        for elem in kwarg.keys():
            str_elem = str_elem + elem + ':' + str(kwarg[elem]) + ' '
        print(str_elem)
                
        return unstable, str_elem, rAPEX_x, rAPEX_y, rAPEX_z, r_Didymoon_x, r_Didymoon_y, r_Didymoon_z, Tn


"2D-plot of the collective trajectories from a single run"

plt.clf()
fig1 = plt.figure(figsize = (15,14))
plt.style.use("fivethirtyeight")
ax1 = fig1.add_subplot(111)
ax1.set_xlim(-2000,2000)
ax1.set_ylim(-2000,2000)
ax1.set_xlabel("x-coordinate")
ax1.set_ylabel("y-coordinate")
ax1.set_title("Orbit of APEX in corotating frame")


"Loop that executes the simulation "

Args = []
Orbit_Times = []
Stable_Orbits = []
Stable_Orbits_Initialconditions = []
itern = 0
stable = 0
for kwargs in kwargs_list:
    print("Number of sets of initial conditions: ", itern, ":",len(kwargs_list))
    args = np.array([kwargs], dtype = object)
    for kwarg in args:
        Var = var.Variables(const.N)
        Initials = meth.Initial_Conditions(Var, **kwarg)
        Initials.initialize_initial_conditions(Var)
        Initials.wrap(Var)
        Args.append((Initials, Var, kwarg))
        "Initializing Kepler Orbit and choosing initial conditions"
        # if Var.Kep == 0:
        #     r_Didymoon_x, r_Didymoon_y, r_Didymoon_z, Tk = meth.Kepler(Initials.E, const.T_max, Var)
        # time.sleep(1)
        itern +=1
        Var.E_end = Initials.E
        start = timeit.timeit()
        try:    
            unstable, str_elem, rAPEX_x, rAPEX_y, rAPEX_z, r_Didymoon_x, r_Didymoon_y, r_Didymoon_z, Tn = Orbital_Calculation(Initials, Var, kwarg)
        except Exception:
            str_elem = ""
            for elem in kwarg.keys():
                str_elem = str_elem + elem + ':' + str(kwarg[elem]) + ' '
            print('Orbit not computable for ' + str_elem + 'continuing with next orbit...')
            
        stop = timeit.timeit()
        time_orbit_calc = stop - start
# with concurrent.futures.ProcessPoolExecutor() as executer:
#     rAPEX_x, rAPEX_y, rAPEX_z, r_Didymoon_x, r_Didymoon_y, r_Didymoon_z, Tn = executer.map(Orbital_Calculation, Args) 
        Orbit_Times.append(time_orbit_calc)
        if unstable == 0:
            stable += 1
            Stable_Orbits.append(Var)
            Stable_Orbits_Initialconditions.append(kwarg)
        ax1.plot(Var.rAPEX_x_corot[0:Var.tn], Var.rAPEX_y_corot[0:Var.tn], linewidth=1, color = 'k')
            
# with open('Stable_Orbits.json', 'w')  as Orbits:
#     for Init_elem in Stable_Orbits_Initialconditions:
#         json.dump(Init_elem, Orbits)
#         Orbits.write('\n')

Didymos_shape2 = patches.Ellipse(xy = (0,0), width = 2*const.a_Did,height = 2*const.b_Did)
Didymoon_shape2 = patches.Circle(xy = (r_Didymoon_x[0],0), radius = const.RadDM)
ax1.add_patch(Didymos_shape2)
ax1.add_patch(Didymoon_shape2)
unstable = len(kwargs_list) - stable

textstr1 = ""
file_name1 = ""

distr = 0
Text_box1 = []
if const.e == 0:
    Text_box1.append('Circular Orbit Didymoon')
else:
    Text_box1.append('Kepler Orbit Didymoon')
if len(kwargs_list)>1:
    distr = 3
    Text_box1.append('Number of stable orbits out of total: ' + str(stable) + ':' + str(len(kwargs_list)))
    Text_box1.append('Number of unstable orbits out of total: ' + str(unstable) + ':' + str(len(kwargs_list)))
if kwarg['Sun'] == 1:
    Text_box1.append('Solar radiation pressure on')
else:
    Text_box1.append('Solar radiation pressure off')


for name in Text_box1:
    textstr1 = textstr1 + name + '\n'


file_name1 = Text_box1[distr][0:5] + Text_box1[distr][-3:] + 'stable_orbits_' + str(stable)
props1 = dict(boxstyle='square', facecolor='white', alpha=1)
ax1.text(0.05, 0.95, textstr1, transform=ax1.transAxes, fontsize=14, verticalalignment='top', bbox = props1, zorder=10)

file_name1 = file_name1 + '.eps'
plt.tight_layout()
plt.show()
fig1.savefig(file_name1, format='eps')
    
"3D-Plot"
    
# def get_proj_scale(self):
#     """                                                                                                                                                                                                                                    
#     Create the projection matrix from the current viewing position.                                                                                                                                                                        

#     elev stores the elevation angle in the z plane                                                                                                                                                                                         
#     azim stores the azimuth angle in the x,y plane                                                                                                                                                                                         

#     dist is the distance of the eye viewing point from the object                                                                                                                                                                          
#     point.                                                                                                                                                                                                                                 

#     """
#     relev, razim = np.pi * self.elev/180, np.pi * self.azim/180

#     xmin, xmax = self.get_xlim3d()
#     ymin, ymax = self.get_ylim3d()
#     zmin, zmax = self.get_zlim3d()

#     # transform to uniform world coordinates 0-1.0,0-1.0,0-1.0                                                                                                                                                                             
#     worldM = proj3d.world_transformation(
#         xmin, xmax,
#         ymin, ymax,
#         zmin, zmax)

#     # look into the middle of the new coordinates                                                                                                                                                                                          
#     R = np.array([0.5, 0.5, 0.5])

#     xp = R[0] + np.cos(razim) * np.cos(relev) * self.dist
#     yp = R[1] + np.sin(razim) * np.cos(relev) * self.dist
#     zp = R[2] + np.sin(relev) * self.dist
#     E = np.array((xp, yp, zp))

#     self.eye = E
#     self.vvec = R - E
#     self.vvec = self.vvec / proj3d.mod(self.vvec)

#     if abs(relev) > np.pi/2:
#     # upside down                                                                                                                                                                                                                          
#       V = np.array((0, 0, -1))
#     else:
#       V = np.array((0, 0, 1))
#     zfront, zback = -self.dist, self.dist

#     viewM = proj3d.view_transformation(E, R, V)
#     perspM = proj3d.persp_transformation(zfront, zback)
#     M0 = np.dot(viewM, worldM)
#     M = np.dot(perspM, M0)

#     return np.dot(M, scale);

# #Make sure these are floating point values:                                                                                                                                                                                              
# scale_x = 12.0
# scale_y = 12.0
# scale_z = 5.0

# #Axes are scaled down to fit in scene                                                                                                                                                                                                    
# max_scale=max(scale_x, scale_y, scale_z)

# scale_x=scale_x/max_scale
# scale_y=scale_y/max_scale
# scale_z=scale_z/max_scale

# #Create scaling matrix                                                                                                                                                                                                                   
# scale = np.array([[scale_x,0,0,0],
#                   [0,scale_y,0,0],
#                   [0,0,scale_z,0],
#                   [0,0,0,1]])
# print(scale)



# Axes3D.get_proj=get_proj_scale

# mpl.rcParams['legend.fontsize'] = 10
# mpl.rcParams['lines.linewidth'] = 2
# fig1 = plt.figure(figsize=(5,5))
# # ax1 = fig.gca(projection='3d')
# # theta = np.linspace(-4 * np.pi, 4 * np.pi, 100)
# # z = np.linspace(-2, 2, 100)
# # r = z**2 + 1
# # x = r * np.sin(theta)
# # y = r * np.cos(theta)
# # ax.plot(x, y, z, label='parametric curve')
# # ax.legend()

# # plt.show()


    
# Didymos_shape1 = patches.Ellipse(xy = (0,0), width = 2*const.a_Did,height = 2*const.b_Did)
# Didymoon_shape1 = patches.Circle(xy = (r_Didymoon_x[0],0), radius = const.RadDM)


# lim_x = max(r_Didymoon_x)
# lim_y = max(r_Didymoon_y)
# lim = max(lim_x,lim_y)
# lim = (-int(lim*1.3),int(lim*1.3))
# # fig1 = plt.figure(figsize=(10,10,10))
# ax1 = fig1.add_subplot(111, projection='3d')
# # ax2 = fig1.add_subplot(212)
# ax1.plot(Var.rAPEX_x_corot[1:Tn], Var.rAPEX_y_corot[1:Tn], Var.rAPEX_z_corot[1:Tn])

# # ax2.plot(Var.rAPEX_x[1:Tn], Var.rAPEX_y[1:Tn], r_Didymoon_x[1:Tn], r_Didymoon_y[1:Tn])
# # ax2.add_patch(Didymoon_shape1)
# # ax2.add_patch(Didymos_shape1)
# ax1.set_xlim(lim)
# ax1.set_ylim(lim)
# ax1.set_zlim(-500,500)
# # ax2.set_xlim(lim)
# # ax2.set_ylim(lim)



# ax1.set_xlabel("x-coordinate")
# ax1.set_ylabel("y-coordinate")
# ax1.set_zlabel("z-coordinate")

# u = np.linspace(0,2*np.pi,100)
# v = np.linspace(0,np.pi,100)

# rDid_x = const.a_Did*np.outer(np.cos(u),np.sin(v))
# rDid_y = const.b_Did*np.outer(np.sin(u),np.sin(v))
# rDid_z = const.c_Did*np.outer(np.ones(np.size(u)),np.cos(v))

# rDidM_x = r_Didymoon_x[0] + const.RadDM*np.outer(np.cos(u),np.sin(v))
# rDidM_y = const.RadDM*np.outer(np.sin(u),np.sin(v))
# rDidM_z = const.RadDM*np.outer(np.ones(np.size(u)),np.cos(v))

# ax1.plot_surface(rDid_x,rDid_y,rDid_z, rstride=4, cstride=4, color='b')
# ax1.plot_surface(rDidM_x,rDidM_y,rDidM_z, rstride=4, cstride=4, color='b')

# fig1.savefig('Test_fig.png', format='png')
# plt.show()
    
    
"Contour Plots"
    
# textstr = ""
# file_name = ""

# cmap = plt.cm.get_cmap("winter")
# cmap.set_under("magenta")
# cmap.set_over("yellow")
    
# f = plt.figure(figsize = (15,14))
# plt.style.use("fivethirtyeight")
# ax = f.add_subplot(111)
# ax.set_xlim(-1600,1600)
# ax.set_ylim(-1600,1600)
# ax.set_xlabel("x-coordinate")
# ax.set_ylabel("y-coordinate")
# ax.set_title("Forcefield and Gravitational Potential")
# xs = np.linspace(-1600, 1600, 2000)
# ys = np.linspace(-1600, 1600, 2000)
# zs = np.linspace(-1600, 1600, 2000)

# xf = np.linspace(-1600, 1600, 80)
# yf = np.linspace(-1600, 1600, 80)
# zf = np.linspace(-1600, 1600, 80)

# XS, YS = np.meshgrid(xs, ys)
# ZS = np.zeros_like(XS)

# XF, YF = np.meshgrid(xf, yf)
# ZF = np.zeros_like(XF)

# corot = 1
# Text_box, F = meth.gradPhi(XF, YF, ZF, Var, grad2 = 1, g_main = 1, g_sec = 1, solar = 0, corot=1)
# # cf = ax.quiver(XF, YF, F[0], F[1])
# Big_Phi = meth.Phi(XS, YS, ZS, Var, grad2 = 1, g_main = 1, g_sec = 1, corot=corot)
# cs = ax.contour(XS, YS, Big_Phi, levels = 100, cmap=cmap, linewidths=0.7) #/np.mean(np.mean(Big_Phi))
# ax.clabel(cs, inline=1, fontsize=10)
# Didymos_shape2 = patches.Ellipse(xy = (0,0), width = 2*const.a_Did,height = 2*const.b_Did)
# Didymoon_shape2 = patches.Circle(xy = (r_Didymoon_x[0],0), radius = const.RadDM)
# ax.add_patch(Didymos_shape2)
# ax.add_patch(Didymoon_shape2)
# for name in Text_box:
#     textstr = textstr + name + '\n'   
#     file_name = file_name + str.split(name, sep = ' ')[0] + '_' + str.split(name, sep = ' ')[1] + '_'

# textstr = textstr[0:-1]
# file_name = file_name[0:-1]
# props = dict(boxstyle='square', facecolor='white', alpha=1)
# ax.text(0.05, 0.95, textstr, transform=ax.transAxes, fontsize=14, verticalalignment='top', bbox = props, zorder=10)
# # ax.plot(Var.rAPEX_x_corot[0:Var.tn], Var.rAPEX_y_corot[0:Var.tn], linewidth=2, color = 'k')

# file_name = file_name + '.eps'
# plt.tight_layout()
# plt.show()
# f.savefig(file_name, format='eps')
        
         
        
        
        
"Forcefield plots"    
# textstr = ""
# file_name = ""

# f = plt.figure(figsize = (10,8))
# plt.style.use("fivethirtyeight")
# ax = f.add_subplot(111)
# ax.set_xlim(-1400,1400)
# ax.set_ylim(-1400,1400)
# ax.set_xlabel("x-coordinate")
# ax.set_ylabel("y-coordinate")
# ax.set_title("Gravitational Potential")
# xs = np.linspace(-1400, 1400, 2000)
# ys = np.linspace(-1400, 1400, 2000)
# zs = np.linspace(-1400, 1400, 2000)

# cmap = plt.cm.get_cmap("winter")
# cmap.set_under("magenta")
# cmap.set_over("yellow")


# X, Y = np.meshgrid(xs, ys)
# Z = np.zeros_like(X)
    
# Text_box, Big_Phi = meth.Phi(X, Y, Z, Var, grad2 = 1, g_main = 1, g_sec = 1)
# cs = ax.contourf(X, Y, Big_Phi, levels = 20, cmap=cmap) #/np.mean(np.mean(Big_Phi))

# Didymos_shape2 = patches.Ellipse(xy = (0,0), width = 2*const.a_Did,height = 2*const.b_Did)
# Didymoon_shape2 = patches.Circle(xy = (r_Didymoon_x[0],0), radius = const.RadDM)
# f.colorbar(cs, ax=ax, shrink=0.9)

# # ax.add_patch(Didymos_shape2)
# # ax.add_patch(Didymoon_shape2)
# for name in Text_box:
#     textstr = textstr + name + '\n'   
#     file_name = file_name + str.split(name, sep = ' ')[0] + '_' + str.split(name, sep = ' ')[1] + '_'

# textstr = textstr[0:-1]
# file_name = file_name[0:-1]
# props = dict(boxstyle='square', facecolor='white', alpha=1)
# ax.text(0.05, 0.95, textstr, transform=ax.transAxes, fontsize=14, verticalalignment='top', bbox = props)
# file_name = file_name + '.eps'
# plt.tight_layout()
# plt.show()
# f.savefig(file_name, format='eps')
    
        
        