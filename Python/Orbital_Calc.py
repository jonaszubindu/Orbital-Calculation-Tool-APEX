# -*- coding: utf-8 -*-
"""
Created on Thu Feb  6 15:53:38 2020

This This is the first dummy setup of the orbital simulation software for single run of orbit test cases.

@author: jonaszbinden
"""


import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import ode
import timeit
from matplotlib import animation

plt.clf()

"Variables and Constants" 
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
tk = 0
tn = 0

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


"APEX coordinates in cartesian"
rAPEX_y, rAPEX_z, rAPEX_x = [np.zeros([1,N+1], dtype = np.float128) for _ in range(3)]
rAPEX_i = np.concatenate((rAPEX_x, rAPEX_y, rAPEX_z), axis = 0) 
vAPEX_x, vAPEX_y, vAPEX_z = [(np.zeros([1,N+1], dtype = np.float128)) for _ in range(3)]
"Spherical coordinates"
rAPEX_theta, rAPEX_phi, rAPEX_r = [np.zeros([1,N+1], dtype = np.float128) for _ in range(3)]

"Didymoon, Sun and solar panel surface vectors in cartesian coordinates"
r_Didymoon_x, r_Didymoon_y, r_Didymoon_z = [np.zeros([1,N+1], dtype = np.float128) for _ in range(3)]
R_Sun_theta, R_Sun_phi, R_Sun_r = [np.zeros([1,N+1], dtype = np.float128) for _ in range(3)]
R_Sun_x, R_Sun_y, R_Sun_z = [np.zeros([1,N+1], dtype = np.float128) for _ in range(3)]  


"ORIGIN OF COORDINATE SYSTEM IS AT CENTER OF MASS OF DIDYMOS AND"
"THETA = 0 IS AT ECLIPTIC OF DIDYMOS AND DIDYMOON"

"Didymoon position from Didymos"
r_Didymoon_theta, r_Didymoon_phi, r_Didymoon_r = [np.zeros([1,N+1], dtype = np.float128) for _ in range(3)]

"First order gravitation terms"
G_main_y, G_main_z, G_main_x = [np.zeros([1,N+1], dtype = np.float128) for _ in range(3)]

G_sec_y, G_sec_z, G_sec_x = [np.zeros([1,N+1], dtype = np.float128) for _ in range(3)]

            
            ####################################
            #                                  #
            # Motion of Didymoon and Didymos   #
            # described through Kepler's Laws  #
            #                                  #
            ####################################

r_mu_r, r_mu_phi, r_mu_theta, tv = [np.zeros([1,N+1], dtype = np.float128) for _ in range(4)]
deltaM = (np.zeros([1,N+1], dtype = np.float128))
deltaE = (np.zeros([1,N+1], dtype = np.float128))
J = (np.zeros([1,N+1], dtype = np.float128))
t_real = (np.zeros([1,N+1], dtype = np.float128))
E_end, Ex_end = [np.zeros([1,N+1], dtype = np.float128) for _ in range(2)]

"From Spherical Coordinates to Cartesian Coordinates" 

def sph_to_cart_coord(r,theta,phi):
    return [r*np.sin(theta)*np.cos(phi), r*np.sin(theta)*np.sin(phi), r*np.cos(theta)]
    
def cart_to_sph(x,y,z):
    r = np.sqrt(x**2+y**2+z**2)
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
    return [r, theta, phi]

"Functions"
P_Sun = 4.56E-6 #N/m^2

"Ilumination Angle"
def ilangle_cos(x1,y1,z1,x2,y2,z2): #illumination angle on a surface
    r1 = np.sqrt(x1**2+y1**2+z1**2)
    r2 = np.sqrt(x2**2+y2**2+z2**2)
    if r1 == 0 or r2 == 0:
        raise ValueError("Vector can not be zero!")
        #print("Vector can not be zero!")
    cos_ill = np.arccos((x1*x2+y1*y2+z1*z2)/(r1*r2))
    return cos_ill

"Shadow Function"
def f_shadow(d,x,y,z,xM,yM,zM): #x,y,z position of APEX, xM,yM,zM position of APEX with respect to Didymoon
    s, theta, phi = cart_to_sph(x,y,z)
    sM, thetaM, phiM = cart_to_sph(xM,yM,zM)
    dM = np.sqrt((d+xM)**2+yM**2+zM**2)
    cosxi = np.cos(np.pi/2-theta)*np.cos(phi)
    cosxiM = np.cos(np.pi/2-thetaM)*np.cos(phiM)
    s0m = s*cosxi
    s0M = sM*cosxiM
    sinf2m = (RadSun - RadDm)/d
    sinf2M = (RadSun - RadDM)/dM
    c2m = RadDm/sinf2m - s0m
    c2M = RadDM/sinf2M - s0M
    f2m = np.arcsin(sinf2m)
    f2M = np.arcsin(sinf2M)
    l2m = c2m*np.tan(f2m)
    l2M = c2M*np.tan(f2M)
    lDm = np.sqrt(s**2 - s0m**2)
    lDM = np.sqrt(sM**2 - s0M**2)
    if lDm > l2m or lDM > l2M:
        vf = 1
    else:
        vf = 0
    #print(vf)
    return vf

"Solar radiation pressure term, componentwise"
def f_sol_A(P, d, A, m, COSilAngle, emiss, eSun, nsurf): #A = surface area, m = mass s/c
    return -P*(1.496E+11/d)**2*(A/m)*COSilAngle*((1-emiss)*eSun+2*emiss*COSilAngle*nsurf)    

   
def G_first(mf,rf,r_fi):
    "Gravitational Field first term:"
    if rf == 0:
        raise ValueError("r is = " + str(rf)) 
        #print("r is = " + str(rf)) 
    G_main_x = -(G*mf/rf**3)*r_fi[0]
    G_main_y = -(G*mf/rf**3)*r_fi[1]
    G_main_z = -(G*mf/rf**3)*r_fi[2]
    return [G_main_x, G_main_y, G_main_z]


tn = 0

"Derivatives: dv/dt = f(v) and dr/dt = v"
def f(tim,v_i):
    gradU_2, R_Sun, e_Sun, n_surf1, n_surf2, drdt, dvdt, r_i,r_sec_i,sol_pan_1,sol_pan_2 = [np.squeeze(np.zeros([1,3], dtype = np.float128)) for _ in range(11)] 
    
    for i in range(3):
        e_Sun[i] = R_Sun[i]/R_Sun_r[0,tn]
    
    for j in range(3):
        r_i[j] = v_i[j+3]
        r_sec_i[j] = r_i[j]-r_Didymoon_i[j][0][tn]
    r, theta, phi = cart_to_sph(r_i[0],r_i[1],r_i[2])
    r_sec = np.sqrt(r_sec_i[0]**2+r_sec_i[1]**2+r_sec_i[2]**2)
    
    R_Sun[0] = R_Sun_x[0,tn]
    R_Sun[1] = R_Sun_y[0,tn]
    R_Sun[2] = R_Sun_z[0,tn]
    
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
    vf = f_shadow(dist_Sun, r_i[0], r_i[1], r_i[2], r_sec_i[0], r_sec_i[1], r_sec_i[2])
    
    "Solar radiation pressure"
    
    for i in range(3):
        sol_pan_1[i] = vf*f_sol_A(P_Sun, dist_Sun, S_1, m_SC, cos_il1, emiss_sol, e_Sun[i], n_surf1[i])
        sol_pan_2[i] = vf*f_sol_A(P_Sun, dist_Sun, S_2, m_SC, cos_il2, emiss_sol, e_Sun[i], n_surf2[i])
    
    "Secondary terms C20 and C22 Didymos"
    
    C_20 = -0.023
    C_22 = -0.0013
    
    gradU_2[0] = 3*G*m_Did*RadDM**2/r**4*(C_20*1/2*3*(np.sin(np.pi/2-theta)**2-1)+C_22*3*np.cos(np.pi/2-theta)**2*np.cos(2*phi))*np.sin(theta)*np.cos(phi) + G*m_Did*RadDM**2/r**4*(C_20*3*np.sin(np.pi/2-theta)*np.cos(np.pi/2-theta) + C_22*6*np.cos(np.pi/2-theta)*np.sin(np.pi/2-theta))*np.cos(2*phi)*np.cos(theta)*np.cos(phi) + 1/np.sin(theta)*(G*m_Did*RadDM**2/r**4*C_22*3*np.cos(np.pi/2-theta)**2*2*np.sin(2*phi))*(-np.sin(phi)) 
    gradU_2[1] = 3*G*m_Did*RadDM**2/r**4*(C_20*1/2*3*(np.sin(np.pi/2-theta)**2-1)+C_22*3*np.cos(np.pi/2-theta)**2*np.cos(2*phi))*np.sin(theta)*np.sin(phi) + G*m_Did*RadDM**2/r**4*(C_20*3*np.sin(np.pi/2-theta)*np.cos(np.pi/2-theta) + C_22*6*np.cos(np.pi/2-theta)*np.sin(np.pi/2-theta))*np.cos(2*phi)*np.cos(theta)*np.sin(phi) + 1/np.sin(theta)*(G*m_Did*RadDM**2/r**4*C_22*3*np.cos(np.pi/2-theta)**2*2*np.sin(2*phi))*np.cos(phi)
    gradU_2[2] = 3*G*m_Did*RadDM**2/r**4*(C_20*1/2*3*(np.sin(np.pi/2-theta)**2-1)+C_22*3*np.cos(np.pi/2-theta)**2*np.cos(2*phi))*np.cos(theta) + G*m_Did*RadDM**2/r**4*(C_20*3*np.sin(np.pi/2-theta)*np.cos(np.pi/2-theta) + C_22*6*np.cos(np.pi/2-theta)*np.sin(np.pi/2-theta)*np.cos(2*phi))*(-np.sin(theta))
    
    
    
    dvdt[0] = sol_pan_1[0] + sol_pan_2[0] + G_main[0] + G_sec[0] + gradU_2[0]
    dvdt[1] = sol_pan_1[1] + sol_pan_2[1] + G_main[1] + G_sec[1] + gradU_2[1]
    dvdt[2] = sol_pan_1[2] + sol_pan_2[2] + G_main[2] + G_sec[2] + gradU_2[2]
    
    
    for i in range(3):
        drdt[i] = v_i[i]
     
    return [dvdt[0], dvdt[1], dvdt[2], drdt[0], drdt[1], drdt[2]]

    
    

"Initial Conditions"
rAPEX_x0 = rAPEX_x[0,0:N] = -100 # m
rAPEX_y0 = rAPEX_y[0,0] = 1000 # m
rAPEX_z0 = rAPEX_z[0,0] = 0 # m

vAPEX_x0 = vAPEX_x[0,0] = -0.2 # m
vAPEX_y0 = vAPEX_y[0,0] = 0 # m
vAPEX_z0 = vAPEX_z[0,0] = 0 # m

vAPEX_i0 = [vAPEX_x0, vAPEX_y0, vAPEX_z0]
vAPEX_i = [vAPEX_x, vAPEX_y, vAPEX_z]

r_Didymoon_x0 = r_Didymoon_x[0,0] = a-a*e
r_Didymoon_y0 = r_Didymoon_y[0,0] = 0
r_Didymoon_z0 = r_Didymoon_z[0,0] = 0

R_Sun_r[0,0:N] = 1.6*AU # 1.6 AU
R_Sun_x[0,0:N] =R_Sun_r[0,0]

vy0 = [vAPEX_x0, vAPEX_y0, vAPEX_z0, rAPEX_x0, rAPEX_y0, rAPEX_z0]        


            ####################################
            #                                  #
            # Differential Equation APEX       #
            # and Time Loop                    #
            #                                  #
            ####################################



Integrator = 1
Kepler = 1


    
v = ode(f).set_integrator('dop853')
v.set_initial_value(vy0, t0)
v.set_integrator('dop853')

        
i = 0
j = 0

r_Didymoon_r[0,0],_ ,_ = cart_to_sph(r_Didymoon_x[0,0],r_Didymoon_y[0,0],r_Didymoon_z[0,0])
arg = (1-r_Didymoon_r[0,0]/a)/e
if arg <= 1 and arg >= -1:
    arg = arg
if arg > 1 and arg < 1.1:
    arg = 1
elif arg < -1 and arg > -1.1:
    arg = -1
else:
    raise TypeError("Error: argument for E invalid")
    #print("Error: argument for E invalid")
E0 = np.arccos(arg)
if r_Didymoon_x[0,0]>=0 and r_Didymoon_y[0,0]>=0:
    E = E0
elif r_Didymoon_x[0,0]>0 and r_Didymoon_y[0,0]<0:
    E = 2*np.Pi - E0
elif r_Didymoon_x[0,0]<0 and r_Didymoon_y[0,0]<0:
    E = E0 + np.Pi
elif r_Didymoon_x[0,0]<0 and r_Didymoon_y[0,0]>=0:
    E = np.pi - E0
else:
    raise ValueError("Error: something went wrong with the initial conditions of Didymoon") 
    #print("Error: something went wrong with the initial conditions of Didymoon")

E_end[0,0] = E
Ex = 10
M0 = E - e*np.sin(E)



"TIME LOOPs"

"Kepler"
if Kepler == 1:
    while tk<(T_max-1):

        tk += 1
        Mx = E - e*np.sin(E)
        M = n*((tk-t0)*dt) + M0
        Ex = M
        start = timeit.timeit()
        deltaT = 0
        j = 0
        while abs(Ex-E)>1E-12:
            if deltaT<1E-1:
                Ex = E
                E = Ex + (M + e*np.sin(Ex)-Ex)/(1-e*np.cos(Ex))
                #print(E)
                j += 1
                end = timeit.timeit()
                deltaT = end - start
            else:
                message = 'Time to compute Kepler too long! i = ' + str(i) + ' j = ' + str(j)
                raise ValueError(message)
                #print(message)
            
        
            
        if E != 0:
            i += 1
            deltaM[0,tk] = abs(M-Mx)
            deltaE[0,tk] = abs(E-Ex)
            E_end[0,tk] = E
            Ex_end[0,tk] = Ex
            J[0,tk] = j 
        
        
        "True anomaly tv"    
        tv[0,tk] = 2*np.arctan(np.sqrt((1+e)/(1-e))*(np.sin(E/2)/(np.cos(E/2))))
        r_Didymoon_r[0,tk] = a*(1-e*np.cos(E))    
        r_Didymoon_x[0,tk] = r_Didymoon_r[0,tk]*np.cos(tv[0,tk])
        r_Didymoon_y[0,tk] = r_Didymoon_r[0,tk]*np.sin(tv[0,tk])
        _, r_Didymoon_phi[0,tk], r_Didymoon_theta[0,tk] = cart_to_sph(r_Didymoon_x[0,tk], r_Didymoon_y[0,tk], r_Didymoon_z[0,tk])
        
    r_Didymoon_R = r_Didymoon_r
    r_Didymoon_X = r_Didymoon_x
    r_Didymoon_Y = r_Didymoon_y
    r_Didymoon_Z = r_Didymoon_z
    r_Didymoon_i = [r_Didymoon_x, r_Didymoon_y, r_Didymoon_z]
    Tk = tk
        
    
"Integration of differential equation" 
if Integrator == 1 :
    tn = 0
    fig = plt.figure()
    ax1 = plt.axes(xlim=(-1500, 1500), ylim=(-2000,2000))
    line, = ax1.plot([], [], lw=2)
    plt.xlabel("x-coordinate")
    plt.ylabel("y-coordinate")
    
    plotlays, plotcols = [2], ["black","red"]
    lines = []
    for index in range(2):
        lobj = ax1.plot([],[],lw=2,color=plotcols[index])[0]
        lines.append(lobj)
    
    
    def init():
        for line in lines:
            line.set_data([],[])
        return lines
    
    
    #plt.plot(r_Didymoon_X[0,0:Tk],r_Didymoon_Y[0,0:Tk]) 
    while tn<(T_max-1) and v.successful():
        
        ######################################
        #                                    #
        #   Define breaking arguments here   #
        #                                    #
        ######################################
        
        
        dist_Sun = R_Sun_r[0,tn]
        
        tn += 1
        v.integrate(v.t+dt)
        
        vAPEX_x[0,tn] = v.y[0]
        vAPEX_y[0,tn] = v.y[1]
        vAPEX_z[0,tn] = v.y[2]
        rAPEX_x[0,tn] = v.y[3]
        rAPEX_y[0,tn] = v.y[4]
        rAPEX_z[0,tn] = v.y[5]
        
        r = np.sqrt(rAPEX_x[0,tn]**2+rAPEX_y[0,tn]**2+rAPEX_z[0,tn]**2)
        r_sec = 0
        r_sec_i = np.zeros([1,3], dtype = np.float128)
        r_sec_i[0,0] = rAPEX_x[0,tn]-r_Didymoon_x[0,tn]
        r_sec_i[0,1] = rAPEX_y[0,tn]-r_Didymoon_y[0,tn]
        r_sec_i[0,2] = rAPEX_z[0,tn]-r_Didymoon_z[0,tn]
        for j in range(3):
            r_sec = r_sec + r_sec_i[0,j]**2
        r_sec = np.sqrt(r_sec)
        t_real[0,tn] = v.t+dt
        if r>10*a:
            print(r, tn)
            plt.plot(rAPEX_x[0,1:tn],rAPEX_y[0,1:tn],r_Didymoon_X[0,1:tn],r_Didymoon_Y[0,1:tn])
            raise Warning("APEX left the system, unstable")
            plt.show()
            
        elif r < RadDm:
            print(r, RadDm, tn)
            plt.plot(rAPEX_x[0,1:tn],rAPEX_y[0,1:tn],r_Didymoon_X[0,1:tn],r_Didymoon_Y[0,1:tn])
            raise Warning("APEX crashed on Didymos")
            plt.show()
            
        elif r_sec < RadDM:
            print(r_sec, RadDM, tn)
            plt.plot(rAPEX_x[0,1:tn],rAPEX_y[0,1:tn],r_Didymoon_X[0,1:tn],r_Didymoon_Y[0,1:tn])
            raise Warning("APEX crashed on Didymoon")
            plt.show()
            

x1,y1 = [],[]
x2,y2 = [],[]

frame_num = tn

def animate(i):

    x = r_Didymoon_X[0, i]
    y = r_Didymoon_Y[0, i]
    x1.append(x)
    y1.append(y)

    a = rAPEX_x[0,i]
    b = rAPEX_y[0,i]
    x2.append(a)
    y2.append(b)

    xlist = [x1, x2]
    ylist = [y1, y2]

    #for index in range(0,1):
    for lnum,line in enumerate(lines):
        line.set_data(xlist[lnum], ylist[lnum]) # set data for each line separately. 

    return lines

# call the animator.  blit=True means only re-draw the parts that have changed.
anim = animation.FuncAnimation(fig, animate, init_func=init,
                               frames=frame_num, interval=10, blit=True)

plt.plot(rAPEX_x[0,1:tn],rAPEX_y[0,1:tn],r_Didymoon_X[0,1:tn],r_Didymoon_Y[0,1:tn])
plt.show()
print(tn)
