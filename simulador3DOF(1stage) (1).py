#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 21 12:11:24 2019

@author: leonardo
"""

#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 21 00:20:26 2019

@author: leonardo
"""
import numpy as np
from scipy import integrate, interpolate 
from scipy.integrate import ode
from gravity import gravity
from atmosphere import atmosphere
from conticap import conticap
import matplotlib.pyplot as plt

def rocket(t, o):
    global dtr; global mu; global omega; global S; global c; global rm; 
    global tb1; global tb2; global fT1; global fT2; global m01; 
    global m02; global mL; global mp1; global mp2; global f8; global Rocket; global deriv   

    lo = o[0]; la = o[1]
    clo = np.cos(lo); slo = np.sin(lo); cla = np.cos(la); sla = np.sin(la)
    fpa = o[4]; chi = o[5]
    cfpa = np.cos(fpa); sfpa = np.sin(fpa); cchi = np.cos(chi); schi = np.sin(chi) 
    
    # Acceleration due to gravity (nonespherical earth)
    g0 = gravity(o[2],o[1])
    g = g0[0]; gn = g0[1]
    
    # Atmospheric properties
    if o[2] < rm:
        o[2] = rm

    alt = o[2] - rm;
    v  = o[3];
    if v < 0.0:
        v = 0.0

    if alt <= 2000.0e3:
        atmosp = atmosphere(alt,v,c)
        rho = atmosp[1]
        mach = atmosp[2] 
        Kn = atmosp[3]
        s = atmosp[4]
        iflow = atmosp[5]
        Qinf = 0.5*rho*(v**2)
        CDC = conticap(mach) 
        CDFM = 1.75 + np.sqrt(np.pi)/(2.0*s) 
        
        
        if iflow == 2:   # Continum flow
            CD = CDC
        elif iflow == 1: # Free-molecule flow 
            CD = CDFM
        else:            # Transition flow
        	CD = CDC + (CDFM - CDC)*(0.333*np.log10(Kn/np.sin(np.pi/6)) + 0.5113)
    else:
        rho = 0.0
        Qinf = 0.0
        CD = 0.0
        mach = 0.0

    if t <= tb1: 
        fT = fT1
        m = m01 - mp1*t/tb1
   
    else:
        fT = 0.0
        m = mL
    
    # Rocket ascent trajectories 
    Rocket = np.array([t, alt, m, mach])
    print Rocket 
    Rocket = np.append(Rocket, [t, alt, m, mach])
    
    D = Qinf*S*CD;
    Xfo = fT - D; Yfo = 0.0; Zfo = 0.0
    
    # Trajectory equations follow
    longidot = o[3]*cfpa*schi/(o[2]*cla) 
    latidot =  o[3]*cfpa*cchi/o[2] 
    raddot = o[3]*sfpa  
    veldot = -g*sfpa + gn*cchi*cfpa + Xfo/m + omega*omega*o[2]*cla*(sfpa*cla - cfpa*cchi*sla)
    if t <= 10.0:
        headdot = 0.0; gammadot = 0.0
    else:
        gammadot = (o[3]/o[2] - g/o[3])*cfpa - gn*cchi*sfpa/o[3] + Zfo/(o[3]*m) + 2.0*omega*schi*cla + omega*omega*o[2]*cla*(cfpa*cla + sfpa*cchi*sla)/o[3]
        if abs(cfpa) > 1e-6:
            headdot = o[3]*schi*np.tan(o[1])*cfpa/o[2] - gn*schi/o[3] - Yfo/(o[3]*cfpa*m) - 2.0*omega*(np.tan(o[4])*cchi*cla - sla) + omega*omega*o[2]*schi*sla*cla/(o[3]*cfpa)
        else:
            headdot = 0.0
            
    # Time derivates    
    deriv = np.array([longidot, latidot, raddot, veldot, gammadot, headdot]) 
    
    if alt <= 120.0e3:
        Qdot = Qinf*v*S*CD/20.0

    #print deriv[3]
    
    return deriv

def main():
    global dtr; dtr = np.pi/180;
    global S; S = 0.114;
    global c; c = 0.25;
    global rm; rm = 6378.135e3;
    global omega; omega = 7.292115090e-5 # Earth rotation 
    global f8; f8 = open('data8.dat', 'w');
    global tb1; tb1 = 22.0
    global fT1; fT1 = 15.0e3 
    global m01; m01 = 254.57 
    global mL;  mL = 10.0
    global mp1; mp1 = 137.72
    global Rocket; global deriv
   
    
    lo = -80.55*dtr; # longitude            
    lat = 28.5*dtr;  # latitude
    rad = rm         # Radius Earth
    fpa = 86.0*dtr;  # flight path angle 
    chi = 170.0*dtr; # Azimuth 
    vel = 0.0        #omega*rm*np.cos(lat)*np.cos(chi); # Velocidade inicial 
    orbinit = np.array([lo, lat, rad, vel, fpa, chi])
    t0 = 0.0; t_final = 180.0
    
    state_ode_f = ode(rocket)
    state_ode_f.set_integrator('dopri5', rtol=5e-9, atol=5e-10, nsteps=50000, max_step=1e-2)
    state2 = orbinit  # initial conditions
    state_ode_f.set_initial_value(state2, 0)
    
    simulation = np.array([t0, state2[0], state2[1], state2[2], state2[3], state2[4], state2[5]], dtype=float)
    while state_ode_f.successful() and state_ode_f.t < (t_final):
        state_ode_f.integrate(t_final, step=True)
        simulation = np.append(simulation, [state_ode_f.t, state_ode_f.y[0], state_ode_f.y[1], state_ode_f.y[2], state_ode_f.y[3], state_ode_f.y[4], state_ode_f.y[5]], axis=0)
        #print("{0:0.8f}\t {1:0.4e} \t{2:10.3e}\t {3:0.3e}\t {4:0.3e}\t {5:0.3e}\t {6:0.3e}\t {7:0.3e}\t {8:0.3e}\t {9}".format(state_ode_f.t, sol[-1, 0]- sol[-2, 0], state_ode_f.y[0], state_ode_f.y[1],  state_ode_f.y[2],  state_ode_f.y[3],  state_ode_f.y[4],  state_ode_f.y[5], state_ode_f.successful()))
    
        plt.figure()
        plt.subplot(3,2,1)
        plt.plot(simulation[:,0], simulation[:,2])
   
    return  

if __name__ == "__main__":
    main()
