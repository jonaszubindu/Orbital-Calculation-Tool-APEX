# -*- coding: utf-8 -*-
"""
Created on Thu Feb  6 15:53:38 2020

This file contains the variables class for Orbit_Test_Simulation

@author: jonaszbinden
""" 


import numpy as np
from classes.constants import N  

deltaM = (np.zeros(N+1, dtype = np.float128))
deltaE = (np.zeros(N+1, dtype = np.float128))
J = (np.zeros(N+1, dtype = np.float128))
tv, Ex_end = [np.zeros(N+1, dtype = np.float128) for _ in range(2)]

class Variables:
    
    def __init__(self, N):  
            
        "Time"
        self.t_real = (np.zeros(N+1, dtype = np.float128))
        
        "APEX coordinates in cartesian"
        self.rAPEX_y, self.rAPEX_z, self.rAPEX_x = [np.zeros(N+1, dtype = np.float128) for _ in range(3)]
        self.rAPEX_i = np.vstack((self.rAPEX_x, self.rAPEX_y, self.rAPEX_z)) 
        self.vAPEX_x, self.vAPEX_y, self.vAPEX_z = [(np.zeros(N+1, dtype = np.float128)) for _ in range(3)]
        "Spherical coordinates"
        self.rAPEX_theta, self.rAPEX_phi, self.rAPEX_r = [np.zeros(N+1, dtype = np.float128) for _ in range(3)]
        
        "Didymoon, Sun and solar panel surface vectors in cartesian coordinates"
        self.r_Didymoon_x, self.r_Didymoon_y, self.r_Didymoon_z = [np.zeros(N+1, dtype = np.float128) for _ in range(3)]
        self.R_Sun_theta, self.R_Sun_phi, self.R_Sun_r = [np.zeros(N+1, dtype = np.float128) for _ in range(3)]
        self.R_Sun_x, self.R_Sun_y, self.R_Sun_z = [np.zeros(N+1, dtype = np.float128) for _ in range(3)]  
        
        
        "ORIGIN OF COORDINATE SYSTEM IS AT CENTER OF MASS OF DIDYMOS AND"
        "THETA = 0 IS AT ECLIPTIC OF DIDYMOS AND DIDYMOON"
        
        "Didymoon position from Didymos"
        self.r_Didymoon_theta, self.r_Didymoon_phi, self.r_Didymoon_r = [np.zeros(N+1, dtype = np.float128) for _ in range(3)]
        self.r_Didymoon_i = np.vstack((self.r_Didymoon_x, self.r_Didymoon_y, self.r_Didymoon_z))
        self.rAPEX_r_corot, self.rAPEX_theta_corot, self.rAPEX_phi_corot = [np.zeros(N+1, dtype = np.float128) for _ in range(3)]
        self.rAPEX_x_corot, self.rAPEX_y_corot, self.rAPEX_z_corot = [np.zeros(N+1, dtype = np.float128) for _ in range(3)]
        
        self.Shadow = 0
        self.Sun = 0
        self.Kep = 0
        self.tn = 0
        self.E_end = 0
        self.Ex_end = 0
        self.t_end = 0
        
        
