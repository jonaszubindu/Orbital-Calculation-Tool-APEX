#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb  6 15:53:38 2020

This file contains all the methods and mathematical functions used in the file 'Orbit_Test_Simulation.py',
including the initialization of the used and stored variables and the assignment of initial conditions for each 
run of the simulation.


@author: jonaszbinden
"""

import os
import numpy as np
# import sys
import timeit
# from matplotlib import pyplot as plt
import pandas as pd
import json

#script_dir = "~/Desktop/Studium/Uppsala/Space\ Mission\ Project\ Work/python"
#sys.path.append(os.path.abspath(script_dir))

from classes import constants as const
from classes import variables as Var


"JSON Encoder and Decoder"
from json import JSONEncoder

#Use NumpyArrayEncoder to write arrays and variables into json files or to read from json files
#For more advanced usage of data to JSON or the other way around, use pandas to_json method
# Example how to use JSON NumpyArrayEncoder:

# numpyArrayOne = numpy.array([[11 ,22, 33], [44, 55, 66], [77, 88, 99]])
# numpyArrayTwo = numpy.array([[51, 61, 91], [121 ,118, 127]])

# Serialization
# numpyData = {"arrayOne": numpyArrayOne, "arrayTwo": numpyArrayTwo}
# print("serialize NumPy array into JSON and write into a file")
# with open("numpyData.json", "w") as write_file:
#     json.dump(numpyData, write_file, cls=NumpyArrayEncoder)
# print("Done writing serialized NumPy array into file")

# Deserialization
# print("Started Reading JSON file")
# with open("numpyData.json", "r") as read_file:
#     print("Converting JSON encoded data into Numpy array")
#     decodedArray = json.load(read_file)

#     finalNumpyArrayOne = numpy.asarray(decodedArray["arrayOne"])
#     print("NumPy Array One")
#     print(finalNumpyArrayOne)
#     finalNumpyArrayTwo = numpy.asarray(decodedArray["arrayTwo"])
#     print("NumPy Array Two")
#     print(finalNumpyArrayTwo)

class NumpyArrayEncoder(JSONEncoder):
    def default(self, obj):
        if isinstance(obj, np.integer):
            return int(obj)
        elif isinstance(obj, np.floating):
            return float(obj)
        elif isinstance(obj, np.ndarray):
            return obj.tolist()
        else:
            return super(NumpyArrayEncoder, self).default(obj)    


"Kepler"
def Kepler(E, T_max, variables):
    
    path = os.getcwd()
    
    try:
        
        os.path.exists(path + '/' + 'Kepler.json')
        with open('Kepler.json', 'r') as Kpl:
            Kepler = json.load(Kpl)
    
        print('Kepler data loaded with Eccentricity: e = {} and Semi-major axis: a = {}'.format(Kepler['Eccentricity e'],Kepler['semi-major axis a']))
        if Kepler['Eccentricity e'] == const.e and Kepler['semi-major axis a'] == const.a:
            variables.r_Didymoon_x = np.asarray(Kepler['r_Didymoon_x']) 
            variables.r_Didymoon_y = np.asarray(Kepler['r_Didymoon_y'])
            variables.r_Didymoon_z = np.asarray(Kepler['r_Didymoon_z'])
            del Kepler
        
            Tk = len(variables.r_Didymoon_x)-1
            for n in range(const.N):
                variables.r_Didymoon_r[n], variables.r_Didymoon_theta[n], variables.r_Didymoon_phi[n] = cart_to_sph(variables.r_Didymoon_x[n], variables.r_Didymoon_y[n], variables.r_Didymoon_z[n]) 
            print('Data for Kepler orbit retrieved from file: Kepler.json')
        else:
            print('Kepler data not matching new Kepler orbit values...')
            raise Exception
    except Exception:
        print('Beginning with calculating Kepler orbit of Didymoon around Didymos')
        tk = 0
        i = 0
        j = 0
        
        variables.E_end = E
        Ex = 10
        M0 = E - const.e*np.sin(E)
        percentage = 0
        while tk<(T_max-1):
            
            "Body of the Kepler orbit computation"
            if tk % int(T_max/10) == 0 and tk > 0:
                percentage += 1
                print('Progress: {} %'.format(percentage*10))
            tk += 1
            Mx = E - const.e*np.sin(E)
            M = const.n*((tk-const.t0)*const.dt) + M0
            Ex = M
            start = timeit.timeit()
            deltaT = 0
            j = 0
            while abs(Ex-E)>1E-12:
                if deltaT<1E-1:
                    Ex = E
                    E = Ex + (M + const.e*np.sin(Ex)-Ex)/(1-const.e*np.cos(Ex))
                    
                    j += 1
                    end = timeit.timeit()
                    deltaT = end - start
                else:
                    message = 'Time to compute Kepler too long! i = ' + str(i) + ' j = ' + str(j)
                    raise ValueError(message)
                    
                
            
                
            if E != 0:
                i += 1
                Var.deltaM[tk] = abs(M-Mx)
                Var.deltaE[tk] = abs(E-Ex)
                variables.E_end = E
                variables.Ex_end = Ex
                Var.J[tk] = j 
                
                
                "True anomaly tv"    
                Var.tv[tk] = 2*np.arctan(np.sqrt((1+const.e)/(1-const.e))*(np.sin(E/2)/(np.cos(E/2))))
                variables.r_Didymoon_r[tk] = const.a*(1-const.e*np.cos(E))    
                variables.r_Didymoon_x[tk] = variables.r_Didymoon_r[tk]*np.cos(Var.tv[tk])
                variables.r_Didymoon_y[tk] = variables.r_Didymoon_r[tk]*np.sin(Var.tv[tk])
                _, variables.r_Didymoon_phi[tk], variables.r_Didymoon_theta[tk] = cart_to_sph(variables.r_Didymoon_x[tk], variables.r_Didymoon_y[tk], variables.r_Didymoon_z[tk])
                    
        
        variables.r_Didymoon_i = [variables.r_Didymoon_x, variables.r_Didymoon_y, variables.r_Didymoon_z]
        variables.Kep = 1
        Tk = tk
        # plt.plot(variables.r_Didymoon_x[1:Tk],variables.r_Didymoon_y[1:Tk])
        # plt.show()
        for n in range(const.N):
            variables.r_Didymoon_r[n], variables.r_Didymoon_theta[n], variables.r_Didymoon_phi[n] = cart_to_sph(variables.r_Didymoon_x[n], variables.r_Didymoon_y[n], variables.r_Didymoon_z[n]) 
        try:
            del Kepler
        except UnboundLocalError:
            pass
        Kepler = {}
        
        
        Kepler['Eccentricity e'] = const.e
        Kepler['semi-major axis a'] = const.a
        Kepler['r_Didymoon_x'] = variables.r_Didymoon_x
        Kepler['r_Didymoon_y'] = variables.r_Didymoon_y
        Kepler['r_Didymoon_z'] = variables.r_Didymoon_z
        
        with open('Kepler.json', 'w') as out_file:
            json.dump(Kepler, out_file, cls=NumpyArrayEncoder, indent = 2)
        print('Kepler Orbit of Didymoon around Didymos successfully calculated, moving on to solving differential equation.')
    

    # fig = plt.figure()
    # ax = fig.add_subplot(111)

    variables.r_Didymoon_i[0] = r_Didymoon_x = variables.r_Didymoon_x
    variables.r_Didymoon_i[1] = r_Didymoon_y = variables.r_Didymoon_y
    variables.r_Didymoon_i[2] = r_Didymoon_z = variables.r_Didymoon_z
    
    
    return r_Didymoon_x, r_Didymoon_y, r_Didymoon_z, Tk

        ####################################
        #                                  #
        # Differential Equation APEX       #
        # and Time Loop                    #
        #                                  #
        ####################################



"Derivatives: dv/dt = f(v) and dr/dt = v"
def fun(tim, v_i, var):
    
    G = const.G
    m_Did = const.m_Did
    m_DidM = const.m_DidM
    RadDm = const.RadDm
    
    
    dist_Sun = var.R_Sun_r[var.tn]  
    r_hat, theta_hat, phi_hat,gradU_2, R_Sun, e_Sun, n_surf1, n_surf2, drdt, dvdt, r_i, r_sec_i, sol_pan_1, sol_pan_2 = [np.squeeze(np.zeros(3, dtype = np.float128)) for _ in range(14)] 
    
    R_Sun[0] = var.R_Sun_r[0] 
    for i in range(3):
        e_Sun[i] = R_Sun[i]/var.R_Sun_r[var.tn]
    
    for j in range(3):
        r_i[j] = v_i[j+3] #r_Didymoon_i append
        r_sec_i[j] = r_i[j]-var.r_Didymoon_i[j][var.tn] # The position of Didymoon does not get re-evaluated in one integration step. To include that, the position of Didymoon must be dynamically calculated as well. 
    r, theta, phi = cart_to_sph(r_i[0],r_i[1],r_i[2])
    r_sec = np.sqrt(r_sec_i[0]**2+r_sec_i[1]**2+r_sec_i[2]**2)
    
    R_Sun[0] = var.R_Sun_x[var.tn]
    R_Sun[1] = var.R_Sun_y[var.tn]
    R_Sun[2] = var.R_Sun_z[var.tn]
    
    n_surf1[0] = n_surf1_x = np.cos(30)*r_i[0]-np.sin(30)*r_i[1]
    n_surf1[1] = n_surf1_y = np.sin(30)*r_i[0]-np.sin(30)*r_i[1]
    n_surf1[2] = n_surf1_z = r_i[2]

    n_surf2[0] = n_surf2_x = np.cos(-30)*r_i[0]-np.sin(-30)*r_i[1]
    n_surf2[1] = n_surf2_y = np.sin(-30)*r_i[0]-np.sin(-30)*r_i[1]
    n_surf2[2] = n_surf2_z = r_i[2]
    
    cos_il1 = ilangle_cos(n_surf1_x,n_surf1_y,n_surf1_z,R_Sun[0],R_Sun[1],R_Sun[2])
    cos_il2 = ilangle_cos(n_surf2_x,n_surf2_y,n_surf2_z,R_Sun[0],R_Sun[1],R_Sun[2])
    
    G_main = G_first(m_Did,r,r_i)
    G_sec = G_first(m_DidM,r_sec,r_sec_i)
    vf = f_shadow(dist_Sun, r_i[0], r_i[1], r_i[2], r_sec_i[0], r_sec_i[1], r_sec_i[2], var.Shadow, var.Sun)
    
    
    
    "Solar radiation pressure"
    
    for i in range(3):
        sol_pan_1[i] = vf*f_sol_A(const.P_Sun, dist_Sun, const.S_1, const.m_SC, cos_il1, const.emiss_sol, e_Sun[i], n_surf1[i])
        sol_pan_2[i] = vf*f_sol_A(const.P_Sun, dist_Sun, const.S_2, const.m_SC, cos_il2, const.emiss_sol, e_Sun[i], n_surf2[i])
    
    
    "Secondary terms C20 and C22 Didymos"
    
    C_20 = -0.023 # calculation from general a and b for the main body? 
    C_22 = -0.0013 # calculation from general a and b for the main body?
    
    
    r_hat[0] = np.sin(theta)*np.cos(phi)
    r_hat[1] = np.sin(theta)*np.sin(phi)
    r_hat[2] = np.cos(theta)
    
    theta_hat[0] = np.cos(theta)*np.cos(phi)
    theta_hat[1] = np.cos(theta)*np.sin(phi)
    theta_hat[2] = -np.sin(theta) 
    
    phi_hat[0] = -np.sin(phi)
    phi_hat[1] = np.cos(phi)
    phi_hat[2] = 0
    
    gradU_2[0] = (3*G*m_Did*RadDm**2*(1/2*C_20*(-1 + 3*np.sin(np.pi/2 - theta)**2) + 3*C_22*np.cos(2*phi)*np.sin(theta)**2))/r**4 
    gradU_2[1] = -(G*m_Did*RadDm**2*(-3*C_20*np.sin(np.pi/2 - theta) + 6*C_22*np.cos(2*phi)*np.cos(theta)*np.sin(theta)))/r**4
    gradU_2[2] = (6*C_22*G*m_Did*RadDm**2*np.sin(2*phi)*np.sin(theta))/r**4
        
    gradU_2_r = gradU_2[0]*r_hat
    gradU_2_theta = gradU_2[1]*theta_hat 
    gradU_2_phi = gradU_2[2]*phi_hat
        
    GradU_2 = gradU_2_r + gradU_2_theta + gradU_2_phi
    
    gradU_2_x = GradU_2[0]
    gradU_2_y = GradU_2[1]
    gradU_2_z = GradU_2[2]
      
    dvdt[0] = sol_pan_1[0] + sol_pan_2[0] + G_main[0] + G_sec[0] + gradU_2_x
    dvdt[1] = sol_pan_1[1] + sol_pan_2[1] + G_main[1] + G_sec[1] + gradU_2_y
    dvdt[2] = sol_pan_1[2] + sol_pan_2[2] + G_main[2] + G_sec[2] + gradU_2_z
    
    
    for i in range(3):
        drdt[i] = v_i[i]
    
    return [dvdt[0], dvdt[1], dvdt[2], drdt[0], drdt[1], drdt[2]]


"Potential Functions"
def Phi(x, y, z, var, grad2 = 1, g_main = 1, g_sec = 1, corot = 1):
    
    G = const.G
    m_Did = const.m_Did
    m_DidM = const.m_DidM
    RadDM = const.RadDM
    RadDm = const.RadDm
    
    x = pd.DataFrame(x)
    y = pd.DataFrame(y)
    z = pd.DataFrame(z)
    r_i = [x, y, z]
    
    
    r_sec_i = []
    
    
    for j in range(3):
        r_sec_i.append(r_i[j]-var.r_Didymoon_i[j][0])
    r, theta, phi = cart_to_sph(r_i[0], r_i[1], r_i[2])
    r_sec = np.sqrt(r_sec_i[0]**2+r_sec_i[1]**2+r_sec_i[2]**2)
    
    
    G_main = -G*m_Did/r
    G_sec = -G*m_DidM/r_sec
    
    
    "Secondary terms C20 and C22 Didymos"
    
    C_20 = -0.023 # unnormalized (paper)
    C_22 = -0.0013 # unnormalized (paper)
    

    U_2 = -(G*m_Did/r)*RadDm**2/r**2*(1/2*(3*np.sin(np.pi/2-theta)**2-1)*C_20 + 3*np.cos(phi)**2*C_22*np.cos(2*phi))

    
    "Corot frame, zentripetal force term added"
    
    omega = 2*np.pi/const.T
    
    zentripet = -1/2*omega**2*r**2


    "Sum of all contributions of the RHS of the differential equation"
    
    # Text_box = []
    
    if grad2 == 0: # grad2 terms off
        U_2 = 0
    else:
        # Text_box.append('GravPotential secondary terms Didymos on')
        pass
        
    if g_main == 0: # Grav Didymos main off
        G_main = 0
    else:
        # Text_box.append('GravPotential main terms Didymos on')   
        pass
        
    if g_sec == 0: # Grav Didymoon off
        G_sec = 0
    else:
        # Text_box.append('GravPotential Didymoon on')
        pass
    
    if corot == 0: # Grav Didymoon off
        # print(zentripet)
        zentripet = 0
        # print(zentripet)
    else:
        # Text_box.append('GravPotential Didymoon on')
        pass
    
    
    Big_Phi = G_main + G_sec + U_2 + zentripet
    
    x_shift = x - var.r_Didymoon_x[0]
    y_shift = y - var.r_Didymoon_y[0]
    z_shift = z - var.r_Didymoon_z[0]
    r_shift = np.sqrt(x_shift**2 + y_shift**2 + z_shift**2)
    filt_DidM = r_shift < RadDM
        
    Big_Phi[r < RadDm] = 0
    Big_Phi[filt_DidM] = 0
    
    
    # dvdt[0] = G_main[0] + G_sec[0] + gradU_2[0]
    # dvdt[1] = G_main[1] + G_sec[1] + gradU_2[0]
    # dvdt[2] = G_main[2] + G_sec[2] + gradU_2[0]
    
    
    # for i in range(3):
    #     drdt[i] = v_i[i]
    
    # return [dvdt[0], dvdt[1], dvdt[2], drdt[0], drdt[1], drdt[2]]
    return Big_Phi


"Force Field around Didymos and Didymoon"
def gradPhi(x, y, z, var, grad2=1, g_main=1, g_sec=1, solar=1, corot=1):
    
    G = const.G
    m_Did = const.m_Did
    m_DidM = const.m_DidM
    RadDM = const.RadDM
    RadDm = const.RadDm
    
    # y = np.zeros_like(x)*400
    x = pd.DataFrame(x)
    y = pd.DataFrame(y)
    z = pd.DataFrame(z)
    r_i = [x, y, z]
    
    dist_Sun = var.R_Sun_r[0]  
    R_Sun, e_Sun = [np.squeeze(np.zeros(3, dtype = np.float128)) for _ in range(2)] 
    
    r_sec_i = []
    sol_pan_1 = []
    sol_pan_2 = []
    
    R_Sun[0] = var.R_Sun_x[0]
    R_Sun[1] = var.R_Sun_y[0]
    R_Sun[2] = var.R_Sun_z[0]
    
    for i in range(3):
        e_Sun[i] = R_Sun[i]/var.R_Sun_r[0]
    
    
    for j in range(3):
        r_sec_i.append(r_i[j]-var.r_Didymoon_i[j][0])
    r, theta, phi = cart_to_sph(r_i[0], r_i[1], r_i[2])
    r_sec = np.sqrt(r_sec_i[0]**2+r_sec_i[1]**2+r_sec_i[2]**2)
    
    n_surf1_x = np.cos(30)*r_i[0]-np.sin(30)*r_i[1]
    n_surf1_y = np.sin(30)*r_i[0]-np.sin(30)*r_i[1]
    n_surf1_z = r_i[2]
    n_surf_1 = [n_surf1_x, n_surf1_y, n_surf1_z]

    n_surf2_x = np.cos(-30)*r_i[0]-np.sin(-30)*r_i[1]
    n_surf2_y = np.sin(-30)*r_i[0]-np.sin(-30)*r_i[1]
    n_surf2_z = r_i[2]
    n_surf_2 = [n_surf2_x, n_surf2_y, n_surf2_z]
    # n_surf_2 = n_surf_1 = -e_Sun
    
    cos_il1 = ilangle_cos(n_surf1_x,n_surf1_y,n_surf1_z,R_Sun[0],R_Sun[1],R_Sun[2])
    cos_il2 = ilangle_cos(n_surf2_x,n_surf2_y,n_surf2_z,R_Sun[0],R_Sun[1],R_Sun[2])
    # cos_il1 = cos_il2 = 1 # APEX in no specific shape
    
    G_main = G_first(m_Did,r,r_i)
    G_sec = G_first(m_DidM,r_sec,r_sec_i)
    vf = f_shadow(dist_Sun, r_i[0], r_i[1], r_i[2], r_sec_i[0], r_sec_i[1], r_sec_i[2], var.Shadow, var.Sun)
    
    
    "Solar radiation pressure"
    
    for i in range(3):
        sol_pan_1.append(vf*f_sol_A(const.P_Sun, dist_Sun, const.S_1, const.m_SC, cos_il1, const.emiss_sol, e_Sun[i], n_surf_1[i]))
        sol_pan_2.append(vf*f_sol_A(const.P_Sun, dist_Sun, const.S_2, const.m_SC, cos_il2, const.emiss_sol, e_Sun[i], n_surf_2[i]))
    
    
    
    "Secondary terms C20 and C22 Didymos"
    
    C_20 = -0.023 # calculation from general a and b for the main body? 
    C_22 = -0.0013 # calculation from general a and b for the main body?
    
    
    r_hat_0 = np.sin(theta)*np.cos(phi)
    r_hat_1 = np.sin(theta)*np.sin(phi)
    r_hat_2 = np.cos(theta)
    
    r_hat = [r_hat_0, r_hat_1, r_hat_2]
    
    theta_hat_0 = np.cos(theta)*np.cos(phi)
    theta_hat_1 = np.cos(theta)*np.sin(phi)
    theta_hat_2 = -np.sin(theta) 
    
    theta_hat = [theta_hat_0, theta_hat_1, theta_hat_2]
    
    phi_hat_0 = -np.sin(phi)
    phi_hat_1 = np.cos(phi)
    phi_hat_2 = 0
    
    phi_hat = [phi_hat_0, phi_hat_1, phi_hat_2]
    
    gradU_2_0 = (3*G*m_Did*RadDm**2*(1/2*C_20*(-1 + 3*np.sin(np.pi/2 - theta)**2) + 3*C_22*np.cos(2*phi)*np.sin(theta)**2))/r**4 
    gradU_2_1 = -(G*m_Did*RadDm**2*(-3*C_20*np.sin(np.pi/2 - theta) + 6*C_22*np.cos(2*phi)*np.cos(theta)*np.sin(theta)))/r**4
    gradU_2_2 = (6*C_22*G*m_Did*RadDm**2*np.sin(2*phi)*np.sin(theta))/r**4
    
    gradU_2_r = [gradU_2_0*rh for rh in r_hat]
    gradU_2_theta = [gradU_2_1*th for th in theta_hat]
    gradU_2_phi = [gradU_2_2*ph for ph in phi_hat]
        
    GradU_2 = [gradU_2_r[i] + gradU_2_theta[i] + gradU_2_phi[i] for i in range(3)]
    
    gradU_2_x = -GradU_2[0]
    gradU_2_y = -GradU_2[1]
    gradU_2_z = -GradU_2[2]
    
    
    "Corot frame, zentripetal force term added"
    
    omega = const.omega
    
    zentripet_x = omega**2*r[0]
    zentripet_y = omega**2*r[1]
    zentripet_z = omega**2*r[2]
    
    
    "Sum of all contributions of the RHS of the differential equation"
    
    Text_box = []
    
    if grad2 == 0: # grad2 terms off
        gradU_2_x = gradU_2_y = gradU_2_z = 0
    else:
        Text_box.append('Gravitation secondary terms Didymos on')
        
    if g_main == 0: # Grav Didymos main off
        G_main[0] = G_main[1] = G_main[2] = 0
    else:
        Text_box.append('Gravitation main terms Didymos on')   
        
    if g_sec == 0: # Grav Didymoon off
        G_sec[0] = G_sec[1] = G_sec[2] = 0
    else:
        Text_box.append('Gravitation Didymoon on')
        
    if solar == 0: # solar radiation pressure off
        sol_pan_1[0] = sol_pan_2[0] = sol_pan_1[1] = sol_pan_2[1] = sol_pan_1[2] = sol_pan_2[2] = 0
    else:
        Text_box.append('Solarradiation pressure on')
    
    if corot == 0: # evaluate in corot frame, add zentripetal force
        zentripet_x = zentripet_y = zentripet_z = 0
    else:
        Text_box.append('Forcefield in Corotating frame')
    
    dvdt_x =  gradU_2_x + G_main[0] + G_sec[0] + sol_pan_1[0] + sol_pan_2[0] + zentripet_x
    dvdt_y =  gradU_2_y + G_main[1] + G_sec[1] + sol_pan_1[1] + sol_pan_2[1] + zentripet_y
    dvdt_z =  gradU_2_z + G_main[2] + G_sec[2] + sol_pan_1[2] + sol_pan_2[2] + zentripet_z
    
    dvdt = [dvdt_x, dvdt_y, dvdt_z] 
    
    
    x_shift = x - var.r_Didymoon_x[0]
    y_shift = y - var.r_Didymoon_y[0]
    z_shift = z - var.r_Didymoon_z[0]
    r_shift = np.sqrt(x_shift**2 + y_shift**2 + z_shift**2)
    filt_DidM = r_shift < RadDM+200
    for elem in dvdt:    
        elem[r < RadDm] = 0
        elem[filt_DidM] = 0
    
    return Text_box, dvdt



"From Spherical Coordinates to Cartesian Coordinates" 
    
def sph_to_cart_coord(r,theta,phi):
    return [r*np.sin(theta)*np.cos(phi), r*np.sin(theta)*np.sin(phi), r*np.cos(theta)]
    
def cart_to_sph(x,y,z):
    r = np.sqrt(x**2+y**2+z**2)
    if np.shape(r) == ():
        if r == 0:
            raise ValueError("vector length = " + str(r))
            #print("vector length = " + str(r))
        if x == 0:
            phi = np.sign(y)*np.pi/2
        elif x<0 and y<0:
            phi = np.arctan(y/x)-np.pi
        elif x<0 and y>=0:
            phi = np.arctan(y/x)+np.pi
        else:
            phi = np.arctan(y/x)
        theta = np.arccos(z/r)
    else:
        phi, theta = [np.zeros_like(r) for _ in range(2)]
        phi = pd.DataFrame(phi)
        theta = pd.DataFrame(theta)
        # r = pd.DataFrame(r)
        # print(type(x),type(y),type(z))
        if np.any(r==0):
            raise ValueError("vector length = " + str(r))
            #print("vector length = " + str(r))
        if np.any(x == 0):
            phi[x == 0] = np.sign(y[x == 0])*np.pi/2
        if np.any(x < 0) and np.any(y < 0):
            filt1 = x < 0
            filt2 = y < 0
            filt = filt1 & filt2
            phi[filt] = np.arctan(y[filt]/x[filt])-np.pi
        if np.any(x < 0) and np.any(y >= 0):
            filt1 = x < 0
            filt2 = y >= 0
            filt = filt1 & filt2
            phi[filt] = np.arctan(y[filt]/x[filt])+np.pi
        if np.any(x > 0):
            filt = x > 0
            phi[filt] = np.arctan(y[filt]/x[filt])
        theta = np.arccos(z/r)
        
    return [r, theta, phi]

"Ilumination Angle"
def ilangle_cos(x1,y1,z1,x2,y2,z2): #illumination angle on a surface
    r1 = np.sqrt(x1**2+y1**2+z1**2)
    r2 = np.sqrt(x2**2+y2**2+z2**2)
    if np.any(r1 == 0) or np.any(r2 == 0):
        raise ValueError("Vector can not be zero!")
        #print("Vector can not be zero!")
    cos_ill = (x1*x2+y1*y2+z1*z2)/(r1*r2)
    return cos_ill

"Shadow Function"
def f_shadow(d,x,y,z,xM,yM,zM,shadow, Sun): #x,y,z position of APEX, xM,yM,zM position of APEX with respect to Didymoon
    if shadow == 1 and Sun == 1:
        
        dM = np.sqrt((d+xM)**2+yM**2+zM**2)
        if np.shape(x) == ():
            
            "Didymos part"
            frac = const.RadSun/const.RadDm
            s0 = d/(frac - 1)
            sin_a = const.RadSun/(d + s0)
            a = np.arcsin(sin_a)
            tan_a = np.tan(a)
            y_sh = (s0-x)*tan_a
            c = np.sqrt(y**2 + z**2)
            
            "Didymoon part"
            fracM = const.RadSun/const.RadDM
            s0M = dM/(fracM - 1)
            sin_aM = const.RadSun/(dM + s0M)
            aM = np.arcsin(sin_aM)
            tan_aM = np.tan(aM)
            y_shM = (s0M-xM)*tan_aM
            cM = np.sqrt(yM**2 + zM**2)
            
            if (np.abs(y_shM) > np.abs(cM) and x < 0 and x > -s0) or (np.abs(y_shM) > np.abs(cM) and xM < 0 and x > -s0M):
                vf = 0
            else:
                vf = 1
            
        else:
            
            "Didymos part"
            frac = const.RadSun/const.RadDm
            s0 = d/(frac - 1)
            sin_a = const.RadSun/(d + s0)
            a = np.arcsin(sin_a)
            tan_a = np.tan(a)
            y_sh = (s0-np.abs(x))*tan_a #cone radius of shadow
            c = np.sqrt(y**2 + z**2)
            
            vf = np.ones_like(x)
            filt1 = np.abs(c) < np.abs(y_sh)
            filt2 = x < 0
            filt3 = x > -s0
            filt = filt2 & filt3 & filt1
            
            "Didymoon part"
            fracM = const.RadSun/const.RadDM
            s0M = dM/(fracM - 1)
            sin_aM = const.RadSun/(dM + s0M)
            aM = np.arcsin(sin_aM)
            tan_aM = np.tan(aM)
            y_shM = (s0M-np.abs(xM))*tan_aM #cone radius of shadow
            cM = np.sqrt(yM**2 + zM**2)
            
            filt1M = np.abs(cM) < np.abs(y_shM)
            filt2M = xM < 0
            filt3M = xM > -s0M
            filtM = filt2M & filt3M & filt1M 
            FILT = filt | filtM
            vf[FILT] = 0
            
            
    elif Sun == 1:
        if np.shape(x) == ():
            vf = np.ones_like(x)
        else:
            vf = 1
    else:
        if np.shape(x) == ():
            vf = np.zeros_like(x)
        else:
            vf = 0
    return vf

"Solar radiation pressure term, componentwise"
def f_sol_A(P, d, A, m, COSilAngle, emiss, eSun, nsurf): #A = surface area, m = mass s/c
    return -P*(1.496E+11/d)**2*(A/m)*COSilAngle*((1-emiss)*eSun+2*emiss*COSilAngle*nsurf)    

   
def G_first(mf, rf,r_fi):
    "Gravitational Field first term:"
    if np.any(rf == 0):
        raise ValueError("r is = " + str(rf)) 
        #print("r is = " + str(rf)) 
    G_main_x = -(const.G*mf/rf**3)*r_fi[0]
    G_main_y = -(const.G*mf/rf**3)*r_fi[1]
    G_main_z = -(const.G*mf/rf**3)*r_fi[2]
    return [G_main_x, G_main_y, G_main_z]
    
    
    
"Initial Conditions"    
class Initial_Conditions:
    
    def __init__(self, variables, x0, y0, z0, vx0, vy0, vz0, Sunr0, Sun = 1, Shadow = 1, plotting = 0):
        "Initial Conditions, self belongs to Initial_Conditions here"
        self.rAPEX_x0 = x0 # m
        self.rAPEX_y0 = y0 # m
        self.rAPEX_z0 = z0 # m
        
        self.vAPEX_x0 = vx0 # m
        self.vAPEX_y0 = vy0 # m
        self.vAPEX_z0 = vz0 # m
        self.Sunr0 = Sunr0
        self.t0 = const.t0
        self.plotting = plotting
    
        self.r_Didymoon_x0 = variables.r_Didymoon_x[0] = const.a-const.a*const.e
        self.r_Didymoon_y0 = variables.r_Didymoon_y[0] = 0
        self.r_Didymoon_z0 = variables.r_Didymoon_z[0] = 0
        self.Shadow = Shadow
        self.Sun = Sun
        
        
    def initialize_initial_conditions(self, variables):
        variables.rAPEX_y[0] = self.rAPEX_y0
        variables.rAPEX_x[0:const.N] = self.rAPEX_x0
        variables.rAPEX_z[0] = self.rAPEX_z0
        
        variables.vAPEX_x[0] = self.vAPEX_x0
        variables.vAPEX_y[0] = self.vAPEX_y0
        variables.vAPEX_z[0] = self.vAPEX_z0
        variables.R_Sun_r[:] = const.a + self.Sunr0*const.AU
        if const.i == 0:
            variables.R_Sun_x[:] = variables.R_Sun_r[:]
        else: 
            variables.R_Sun_x[:] = variables.R_Sun_r[:]*np.cos(np.pi-const.i) #inclination dependent ilumination from the sun
            variables.R_Sun_z[:] = variables.R_Sun_r[:]*np.sin(np.pi-const.i) #inclination dependent ilumination from the sun
        variables.Shadow = self.Shadow
        variables.Sun = self.Sun
        
        
    def wrap(self, variables):
        self.vy0 = np.array([self.vAPEX_x0, self.vAPEX_y0, self.vAPEX_z0, self.rAPEX_x0, self.rAPEX_y0, self.rAPEX_z0], dtype = float)
        variables.r_Didymoon_r[0],_ ,_ = cart_to_sph(variables.r_Didymoon_x[0], variables.r_Didymoon_y[0], variables.r_Didymoon_z[0])
        if const.e == 0 and variables.r_Didymoon_r[0] == const.a:
            arg = 0
        else:
            arg = (1-variables.r_Didymoon_r[0]/const.a)/const.e
        if arg <= 1 and arg >= -1:
            arg = arg
        elif arg > 1 and arg < 1.1:
            arg = 1
        elif arg < -1 and arg > -1.1:
            arg = -1
        else:
            raise TypeError("Error: argument for E invalid")
            #print("Error: argument for E invalid")
        E0 = np.arccos(arg)
        if self.r_Didymoon_x0>=0 and self.r_Didymoon_y0>=0:
            self.E = E0
        elif self.r_Didymoon_x0>0 and self.r_Didymoon_y0<0:
            self.E = 2*np.Pi - E0
        elif self.r_Didymoon_x0<0 and self.r_Didymoon_y0<0:
            self.E = E0 + np.Pi
        elif self.r_Didymoon_x0<0 and self.r_Didymoon_y0>=0:
            self.E = np.pi - E0
        else:
            raise ValueError("Error: something went wrong with the initial conditions of Didymoon")

def corot_frame(x, y, z, rDidymoon_phi):
    r, theta, phi = cart_to_sph(x, y, z)
    r_new = r
    theta_new = theta
    phi_new = phi-rDidymoon_phi
    return r_new, theta_new, phi_new
    
    

