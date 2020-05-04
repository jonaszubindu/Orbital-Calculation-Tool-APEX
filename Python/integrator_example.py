#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec 18 13:40:16 2019

@author: jonaszbinden
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import ode
    
G = 6.672*(10**-11)
M = 5.972*(10**24)
mu = G*M

def f(t, y):
    "Returns dy/dt"
    r = np.sqrt(y[0]**2 + y[1]**2 + y[2]**2)

    dy0 = y[3]
    dy1 = y[4]
    dy2 = y[5]
    dy3 = -(mu / (r**3)) * y[0]
    dy4 = -(mu / (r**3)) * y[1]
    dy5 = -(mu / (r**3)) * y[2]
    return [dy0, dy1, dy2, dy3, dy4, dy5]


t = 0.
y0 = [10.e6, 0., 0., 0., np.sqrt(G * M / 10e6), 0.]

integrator = ode(f)
integrator.set_integrator('dop853')
integrator.set_initial_value(y0, t)
plt.close('all')

norm = plt.Normalize(0., 8640.)
cmap = plt.get_cmap('viridis')
#map(norm,integrator.t)

while integrator.successful():
    if integrator.t > 8640.: 
        break
    integrator.integrate(integrator.t + 10.)
    plt.plot(integrator.y[0]/1e6, integrator.y[1]/1e6, 'k.', color=cmap(norm(integrator.t)))


plt.show() 