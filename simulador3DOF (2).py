#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 21 00:20:26 2019

@author: leonardo
"""
import numpy as np
from scipy.integrate import ode, odeint
import myFun as mf
import matplotlib.pyplot as plt

def rocket(t,o):
    global dtr; global mu; global omega; global S; global c; global rm; 
    global tb1; global tb2; global fT1; global fT2; global m01; 
    global m02; global mL; global mp1; global mp2; global f8   

    lo = o[0]; la = o[1]
    clo = np.cos(lo); slo = np.sin(lo); cla = np.cos(la); sla = np.sin(la)
    fpa = o[4]; chi = o[5]
    cfpa = np.cos(fpa); sfpa = np.sin(fpa); cchi = np.cos(chi); schi = np.sin(chi) 
    
    # Acceleration due to gravity (nonespherical earth)
    g0 = mf.gravity(o[2],o[1])
    g = g0[0]; gn = g0[1]
    
    # Atmospheric properties
    if o[2] < rm:
        o[2] = rm

    alt = o[2] - rm;
    v  = o[3];
    if v < 0.0:
        v = 0.0

    if alt <= 2000.0e3:
        atmosp = mf.atmosphere(alt,v,c)
        rho = atmosp[1]
        mach = atmosp[2] 
        Kn = atmosp[3]
        s = atmosp[4]
        iflow = atmosp[5]
        Qinf = 0.5*rho*(v**2)
        CDC = mf.conticap(mach) 
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
        CD = 8.0*CD # Cd first stage = 8*Cd_p(Cd of payload)
    elif t <= (tb1+tb2):
        fT = fT2 
        m = m02 - mp2*(t-tb1)/tb2
        CD = 3.0*CD # Cd second stage = 3*Cd_p(Cd of payload)
    else:
        fT = 0.0
        m = mL
    
    # Rocket ascent trajectories    
    Rocket = np.array([t, alt/1e3, m, mach])
    print (Rocket)
    D = Qinf*S*CD;
    Xfo = fT - D; Yfo = 0.0; Zfo = 0.0
    
    # Trajectory equations follow
    longidot = o[3]*cfpa*schi/(o[2]*cla) 
    latidot =  o[3]*cfpa*cchi/o[2] 
    raddot = o[3]*sfpa  
    veldot = -g*sfpa + gn*cchi*cfpa + Xfo/m + omega*omega*o[2]*cla*(sfpa*cla - cfpa*cchi*sla)
    if t <= 10.0:
        headdot = 0.0; 
        gammadot = 0.0
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
    
    return deriv 

def main():
    global dtr; dtr = np.pi/180;
    global S; S = 0.114;
    global c; c = 0.25;
    global rm; rm = 6378.135e3;
    global omega; omega = 7.292115090e-5 # Earth rotation 
    global f8; f8 = open('data8.dat', 'w');
    global tb1; tb1 = 50.0
    global tb2; tb2 = 87.5
    global fT1; fT1 = 346017.7314
    global fT2; fT2 = 249640.1712
    global m01; m01 = 16528.42
    global m02; m02 = 7046.715
    global mL;  mL = 350.0
    global mp1; mp1 = 8817.985
    global mp2; mp2 = 6361.88
   
    
    lo = -80.55*dtr; # longitude            
    lat = 28.5*dtr;  # latitude
    rad = rm         # Radius Earth
    fpa = 90.0*dtr;  # flight path angle 
    chi = 170.0*dtr; # Azimuth 
    vel = 0.0#np.sqrt((omega*rm*np.cos(lat)*np.cos(chi))**2) # Velocidade inicial 
    orbinit = np.array([lo, lat, rad, vel, fpa, chi])    # initial conditions

    t0 = 0.0; t_final = 100
    
    state_ode_f = ode(rocket)
    state_ode_f.set_integrator('dopri5', rtol=5e-9, atol=5e-10, nsteps=1e6, max_step=1e-2) 
    state_ode_f.set_initial_value(orbinit, 0)
    
    sol = np.array([t0, orbinit[0], orbinit[1], orbinit[2], orbinit[3], orbinit[4], orbinit[5]], dtype=float)
    while state_ode_f.successful() and state_ode_f.t <= (t_final):
        state_ode_f.integrate(t_final, step=True)
        sol = np.append(sol, [state_ode_f.t, state_ode_f.y[0], state_ode_f.y[1], state_ode_f.y[2], state_ode_f.y[3], state_ode_f.y[4], state_ode_f.y[5]], axis=0)
        #print("{0:0.8f}\t {1:0.4e} \t{2:10.3e}\t {3:0.3e}\t {4}".format(state_ode_f.t, sol[-1, 0]- sol[-2, 0], state_ode_f.y[0], state_ode_f.y[1], state_ode_f.successful()))

    plt.figure(1)
    plt.subplot(311)
    plt.xlabel('Time [s]')
    plt.ylabel(r'$\lambda$ [deg]')
    plt.plot(sol[:,0], sol[:,1]/dtr)
    plt.grid(True)
    
    plt.subplot(312)
    plt.xlabel('Time [s]')
    plt.ylabel(r'$\delta$ [deg]')
    plt.plot(sol[:,0], sol[:,2]/dtr)
    plt.grid(True)
    
    plt.subplot(313)
    plt.xlabel('Time [s]')
    plt.ylabel('Altitude [km]')
    plt.plot(sol[:,0], (sol[:,3] - rm)/1e3)
    plt.grid(True)

    plt.figure(2)
    plt.subplot(311)
    plt.xlabel('Time [s]')
    plt.ylabel('Velocity [m/s]')
    plt.plot(sol[:,0], sol[:,4])
    plt.grid(True)
    
    plt.subplot(312)
    plt.xlabel('Time [s]')
    plt.ylabel(r'$\phi$ [deg]')
    plt.plot(sol[:,0], sol[:,5]/dtr)
    plt.grid(True)

    plt.subplot(313)
    plt.xlabel('Time [s]')
    plt.ylabel(r'$A$ [deg]')
    plt.plot(sol[:,0], sol[:,6]/dtr)
    plt.grid(True)

    return sol 
    """
    t = np.linspace(0.0, 1000.0, 10000.0)
    simulation = odeint(rocket, orbinit, t, rtol=5e-10, atol=5e-15)    
    
    plt.figure(1)
    plt.subplot(311)
    plt.xlabel('Time [s]')
    plt.ylabel(r'$\lambda$ [deg]')
    plt.plot(t, simulation[:,0]/dtr)
    plt.grid(True)
    
    plt.subplot(312)
    plt.xlabel('Time [s]')
    plt.ylabel(r'$\delta$ [deg]')
    plt.plot(t, simulation[:,1]/dtr)
    plt.grid(True)
    
    plt.subplot(313)
    plt.xlabel('Time [s]')
    plt.ylabel('Altitude [km]')
    plt.plot(t, (simulation[:,2] - rm)/1e3)
    plt.grid(True)

    plt.figure(2)
    plt.subplot(311)
    plt.xlabel('Time [s]')
    plt.ylabel('Velocity [m/s]')
    plt.plot(t, simulation[:,3])
    plt.grid(True)
    
    plt.subplot(312)
    plt.xlabel('Time [s]')
    plt.ylabel(r'$\phi$ [deg]')
    plt.plot(t, simulation[:,4]/dtr)
    plt.grid(True)

    plt.subplot(313)
    plt.xlabel('Time [s]')
    plt.ylabel(r'$A$ [deg]')
    plt.plot(t, simulation[:,5]/dtr)
    plt.grid(True)
       
    return simulation
    """
if __name__ == "__main__":
    main()
