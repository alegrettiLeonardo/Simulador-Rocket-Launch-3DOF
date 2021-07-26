import numpy as np
from scipy import integrate
from gravity import gravity
from atmosphere import atmosphere
from conticap import conticap
import matplotlib.pyplot as plt

def rocket(o, t):
    global dtr; global mu; global omega; global S; global c; global rm; 
    global tb1; global tb2; global fT1; global fT2; global m01; global me;
    global m02; global mL; global mp1; global mp2; global f8; global Rocket; 

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
        CD = CD # Cd first stage = 8*Cd_p(Cd of payload)
    elif t <= (tb1+tb2):
        fT = fT2 
        m = m02 - mp2*(t-tb1)/tb2
        CD = CD # Cd second stage = 3*Cd_p(Cd of payload)
    else:
        fT = 0.0
        m = mL
    """    
    if t <= tb1: 
        fT = fT1
        m = m01 - mp1*t/tb1
        CD = CD 
    
    elif t > tb1 and t<tb1 + 12.0:
        fT = 0.0
        m = m02
        CD = CD
    elif t <= (tb1+tb2):
        fT = fT2 
        m = m02 - mp2*(t-tb1)/tb2
        CD = CD 
    else:
        fT = 0.0
        m = mL
    """
    # Rocket ascent trajectories 
    Rocket = np.array([t, alt, m, mach])
    print(Rocket) 


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
        if abs(cfpa) > 1e-8:
            headdot = o[3]*schi*np.tan(o[1])*cfpa/o[2] - gn*schi/o[3] - Yfo/(o[3]*cfpa*m) - 2.0*omega*(np.tan(o[4])*cchi*cla - sla) + omega*omega*o[2]*schi*sla*cla/(o[3]*cfpa)
        else:
            headdot = 0.0
            
    # Time derivates    
    deriv = np.array([longidot, latidot, raddot, veldot, gammadot, headdot]) 
    Rocket = np.append(Rocket, [t, alt, m, mach])
    
    if alt <= 120.0e3:
        Qdot = Qinf*v*S*CD/20.0

    #print deriv[3]
    
    return deriv

def main():
    global dtr; dtr = np.pi/180;
    global S; S = np.pi*(3.7**2)/4;
    global c; c = 0.25;
    global rm; rm = 6378.135e3;
    global omega; omega = 7.292115090e-5 # Earth rotation 
    global f8; f8 = open('data8.dat', 'w');
    global tb1; tb1 = 162.0
    global tb2; tb2 = 397.0
    global fT1; fT1 = 7607.0e3
    global fT2; fT2 = 934.0e3 
    global m01; m01 = 511170.0
    global m02; m02 = 115470.0
    global mL;  mL = 13150.0
    global mp1; mp1 = 395700.0
    global mp2; mp2 = 92670.0
   
   
    lo = -44.3676*dtr; # longitude            
    lat = -2.3173*dtr; # latitude
    rad = rm           # Radius Earth
    fpa = 90.0*dtr;    # flight path angle 
    chi = 170.0*dtr;   # Azimuth 
    vel = 0.0          # Velocidade inicial 
    orbinit = np.array([lo, lat, rad, vel, fpa, chi])
    
    t = np.linspace(0.0, 1500.0, 15000000.0)
    simulation = integrate.odeint(rocket, orbinit, t, rtol=5e-8, atol=5e-10)    

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

    #plt.figure(3)
    #plt.plot(Rocket[:,0], Rocket[:,2])
    
    return simulation 
    
if __name__ == "__main__":
    main()