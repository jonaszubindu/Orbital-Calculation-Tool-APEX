# -*- coding: utf-8 -*-
    
"""
Created on Thu Feb  6 15:53:38 2020

This file contains all the necessary constants for Orbit_Test_Simulation

@author: jonaszbinden
"""


import numpy as np
    
  
"Constants" 
G = 6.67430E-11 #N*m^2/kg^2
AU = 1.496E+11 # 1 AU astronomical unit in m

"Mass APEX, mass Didymos mDid, mass Didymoon mDidM"
m_SC = 12
m_Did = 5.24E+11
m_DidM = 3.45E+9
# m_Did = 1E+12
# m_DidM = 8E+11


"mean motion n, period T, semi-major aX_is a, eccentricity e, inclination i, orbital pole A and B orb_pole_A, orb_pole_b"
T = 11.92*3600 #s
a = 1.178E+3 #m
n = np.sqrt(G*(m_DidM+m_Did)/a**3)
e = 0.05 # 0.05 in literature
b = np.sqrt(a**2*(1-e**2))
p = a*(1-e**2)
orb_pole_A = (300,-60)
orb_pole_B = (310,-84)
i = 0 #none given in literature

a_Did = 395 # m
b_Did = 390 # m
c_Did = 370 # m
omega = 2*np.pi/T

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
    
    
    
