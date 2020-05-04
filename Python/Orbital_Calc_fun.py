# -*- coding: utf-8 -*-
    

"Student Project Work Jonas Zbinden, APEX orbital calculations"

import numpy as np
import timeit

class const():
    "Variables and const" 
    G = 6.67430E-11 #N*m^2/kg^2
    AU = 1.496E+11 # 1 AU astronomical unit in m
    
    "Mass APEX, mass Didymos mDid, mass Didymoon mDidM"
    m_SC = 12
    m_Did = 5.24E+11
    m_DidM = 3.45E+9
    
    
    "mean motion n, period T, semi-major aX_is a, eccentricity e"
    "mean anomaly M, true anomaly E"
    T = 11.92*3600 #s
    a = 1.178E+3 #m
    n = np.sqrt(G*(m_DidM+m_Did)/a**3)
    e = 0.05
    b = np.sqrt(a**2*(1-e**2))
    p = a*(1-e**2)
    
    "Time steps"
    dt = 10 #s
    rev = 2 #number of revolutions of Didymoon
    T_max = int(rev*T/dt)
    t0 = 0
    N = T_max
    
    
    "Solarpanel Surface vector 1 + 2"
    S_1 = 0.0828258 #m^2, surface area solar panel 1
    S_2 = S_1 #surface area solar panel 2
    SC_1 = S_2 #surface area s/c 1
    SC_2 = SC_1 #surface area s/c 2
    SC_back = 0.0366 #m^2, surface area s/c back
    SC_top = 0.02263 #m^2, surface area s/c top
    SC_bottom = SC_top #surface area s/c bottom
    emiss_sol = 0.21
    emiss_al = 0.88
    RadSun = 696340000 #km
    RadDm = 385 #m
    RadDM = 75 #m
    P_Sun = 4.56E-6 #N/m^2


class init_Variables:
    def __init__(self,N):
        self.N = N
        
    "APEX coordinates in cartesian"
    rAPEX_y, rAPEX_z, rAPEX_x = [np.zeros([1,self.N+1], dtype = np.float128) for _ in range(3)]
    rAPEX_i = np.concatenate((rAPEX_x, rAPEX_y, rAPEX_z), axis = 0) 
    vAPEX_x, vAPEX_y, vAPEX_z = [(np.zeros([1,self.N+1], dtype = np.float128)) for _ in range(3)]
    "Spherical coordinates"
    rAPEX_theta, rAPEX_phi, rAPEX_r = [np.zeros([1,self.N+1], dtype = np.float128) for _ in range(3)]
    
    "Didymoon, Sun and solar panel surface vectors in cartesian coordinates"
    r_Didymoon_x, r_Didymoon_y, r_Didymoon_z = [np.zeros([1,self.N+1], dtype = np.float128) for _ in range(3)]
    R_Sun_theta, R_Sun_phi, R_Sun_r = [np.zeros([1,self.N+1], dtype = np.float128) for _ in range(3)]
    R_Sun_x, R_Sun_y, R_Sun_z = [np.zeros([1,self.N+1], dtype = np.float128) for _ in range(3)]  
    
    
    "ORIGIN OF COORDINATE SYSTEM IS AT CENTER OF MASS OF DIDYMOS AND"
    "THETA = 0 IS AT ECLIPTIC OF DIDYMOS AND DIDYMOON"
    
    "Didymoon position from Didymos"
    r_Didymoon_theta, r_Didymoon_phi, r_Didymoon_r = [np.zeros([1,self.N+1], dtype = np.float128) for _ in range(3)]
    
    "First order gravitation terms"
    G_main_y, G_main_z, G_main_x = [np.zeros([1,self.N+1], dtype = np.float128) for _ in range(3)]
    
    G_sec_y, G_sec_z, G_sec_x = [np.zeros([1,self.N+1], dtype = np.float128) for _ in range(3)]
    
            
            
    
    r_mu_r, r_mu_phi, r_mu_theta, tv = [np.zeros([1,self.N+1], dtype = np.float128) for _ in range(4)]
    deltaM = (np.zeros([1,self.N+1], dtype = np.float128))
    deltaE = (np.zeros([1,self.N+1], dtype = np.float128))
    J = (np.zeros([1,self.N+1], dtype = np.float128))
    t_real = (np.zeros([1,self.N+1], dtype = np.float128))
    E_end, Ex_end = [np.zeros([1,self.N+1], dtype = np.float128) for _ in range(2)]
    
    
     


           

            ####################################
            #                                  #
            # Motion of Didymoon and Didymos   #
            # described through Kepler's Laws  #
            #                                  #
            ####################################






