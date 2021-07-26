import numpy as np
from scipy import integrate
from gravity import gravity
from atmosphere import atmosphere
from conticap import conticap
import matplotlib.pyplot as plt
from matplotlib.pyplot import  xticks, gca, yticks

def rocket(o, t):
    global dtr; global mu; global omega; global S; global c; global rm; 
    global tb1; global tb2; global fT1; global fT2; global m01; global me;
    global m02; global mL; global mp1; global mp2; global Rocket; global deriv   

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
        m = mL + me
    
    # Rocket ascent trajectories
    f = open('rocket.txt','w')
    Rocket = np.array([t, alt/1e3, m, mach])
    Rocket = np.append(Rocket, [t, alt/1e3, m, mach])
    f.write('%.2f   %.2f    %.2f    %.2f '%(Rocket[0], Rocket[1], Rocket[2], Rocket[3],))
    
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
    
    return deriv

def main():
    global dtr; dtr = np.pi/180;
    global S; S = 0.114;
    global c; c = 0.5;
    global rm; rm = 6378.135e3;
    global omega; omega = 7.292115090e-5 # Earth rotation 
    global tb1; tb1 = 22.0
    global fT1; fT1 = 15.0e3 
    global m01; m01 = 254.57
    global mL;  mL = 10.0
    global mp1; mp1 = 137.72
    global me; me = 106.85
   
   
    lo = -44.3676*dtr; # longitude            
    lat = -2.3173*dtr;  # latitude
    rad = rm         # Radius Earth
    fpa = 86.0*dtr;  # flight path angle 
    chi = 45.0*dtr; # Azimuth 
    vel = 0.0        # Velocidade inicial 
    orbinit = np.array([lo, lat, rad, vel, fpa, chi])
    
    t = np.linspace(0.0, 205.0, 205000.0)
    simulation = integrate.odeint(rocket, orbinit, t, rtol=5e-10, atol=5e-15)    
    """        
    plt.figure(1, figsize = (60,60))
    #gca().set_aspect('equal', adjustable='box')
    plt.tick_params(axis='x', labelsize=15)
    plt.tick_params(axis='y', labelsize=15)
    plt.xlabel('Time [s]')
    plt.ylabel(r'$\lambda$ [deg]', fontsize = 20)
    plt.plot(t, simulation[:,0]/dtr)
    plt.grid(True)
    
    plt.figure(2, figsize = (60,60))
    #gca().set_aspect('equal', adjustable='box')
    plt.tick_params(axis='x', labelsize=15)
    plt.tick_params(axis='y', labelsize=15)
    plt.xlabel('Tempo [s]', fontsize = 20)
    plt.ylabel(r'$\delta$ [deg]', fontsize = 20)
    plt.plot(t, simulation[:,1]/dtr)
    plt.grid(True)
    
    plt.figure(3, figsize = (60,60))
    #gca().set_aspect('equal', adjustable='box')
    plt.tick_params(axis='x', labelsize=15)
    plt.tick_params(axis='y', labelsize=15)
    plt.xlabel('Time [s]')
    plt.ylabel(r'$\phi$ [deg]', fontsize = 20)
    plt.plot(t, simulation[:,4]/dtr)
    plt.grid(True)

    plt.figure(4, figsize = (60,60))
    #gca().set_aspect('equal', adjustable='box')
    plt.tick_params(axis='x', labelsize=15)
    plt.tick_params(axis='y', labelsize=15)
    plt.xlabel('Tempo [s]', fontsize = 20)
    plt.ylabel(r'$A$ [deg]', fontsize = 20)
    plt.plot(t, simulation[:,5]/dtr)
    plt.grid(True)
    
    plt.figure(5, figsize = (60,60))
    #gca().set_aspect('equal', adjustable='box')
    plt.tick_params(axis='x', labelsize=15)
    plt.tick_params(axis='y', labelsize=15)
    plt.xlabel('Time [s]', fontsize = 18)
    plt.ylabel('Altitude [km]', fontsize = 20)
    plt.plot(t, (simulation[:,2] - rm)/1e3)
    plt.grid(True)
    
    plt.figure(6, figsize = (60,60))
    #gca().set_aspect('equal', adjustable='box')
    plt.tick_params(axis='x', labelsize=15)
    plt.tick_params(axis='y', labelsize=15)
    plt.xlabel('Tempo [s]', fontsize = 20)
    plt.ylabel('Velocidade [m/s]', fontsize = 20)
    plt.plot(t, simulation[:,3])
    plt.grid(True)
    """
    plt.figure(7, figsize = (60,60))
    #gca().set_aspect('equal', adjustable='box')
    plt.tick_params(axis='x', labelsize=15)
    plt.tick_params(axis='y', labelsize=15)
    plt.xlabel(r'$\lambda$ [deg]', fontsize = 20)
    plt.ylabel(r'$\delta$ [deg]', fontsize = 20)
    plt.plot(simulation[:,0]/dtr, simulation[:,1]/dtr)
    plt.grid(True)
    
    range_x = -((rm*np.cos(simulation[:,1]/dtr)*np.cos(simulation[:,0]/dtr)) - rm)/1e3
    range_y = -((rm*np.cos(simulation[:,1]/dtr)*np.sin(simulation[:,0]/dtr)) - rm)/1e3
    
    
    plt.figure(8, figsize = (60,60))
    #gca().set_aspect('equal', adjustable='box')
    plt.tick_params(axis='x', labelsize=15)
    plt.tick_params(axis='y', labelsize=15)
    plt.title(r'Alcance em Coordenadas Cartesianas',fontsize = 20)
    plt.xlabel(r'Eixo x [km]', fontsize = 20)
    plt.ylabel(r'Eixo y [km]', fontsize = 20)
    plt.plot(range_y, range_x)
    plt.grid(True)

    return simulation 

if __name__ == "__main__":
    main()