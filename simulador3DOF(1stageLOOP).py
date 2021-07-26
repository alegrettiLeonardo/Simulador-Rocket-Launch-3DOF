import numpy as np
from scipy import integrate
from gravity import gravity
from atmosphere import atmosphere
from conticap import conticap
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from mpl_toolkits.mplot3d import Axes3D
import matplotlib as mpl

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
    #print Rocket
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
   
    global dtr;
    global S;
    global c; 
    global rm; 
    global omega; 
    global tb1; 
    global fT1; 
    global m01; m0 = [305.49, 277.72, 254.57]
    global mL; 
    global mp1; 
    global me; ms = [157.77, 129.99, 106.85] 
   
    # Storage matrices initialization 
    varZ = [210000, 3]
    lo = np.zeros(varZ)
    lat = np.zeros(varZ)
    rad = np.zeros(varZ)
    fpa = np.zeros(varZ)
    chi = np.zeros(varZ)
    vel = np.zeros(varZ)
    range_x = np.zeros(varZ)
    range_y = np.zeros(varZ)
    range_z = np.zeros(varZ)
    Range = np.zeros(varZ)
    t = np.linspace(0, 210, 210000)
    
    for i in np.arange(0,3):
        dtr = np.pi/180
        S = 0.114
        c = 0.5
        rm = 6378.135e3
        omega = 7.292115090e-5 # Earth rotation
        tb1 = 22.0
        fT1 = 15.0e3
        mL = 10.0
        mp1 = 137.72
        me = ms[i]
        m01 = m0[i]

        loi = -44.3676*dtr # longitude            
        lati = -2.3173*dtr # latitude
        radi = rm         # Radius Earth
        fpai = 86.0*dtr   # flight path angle 
        chii = 45.0*dtr   # Azimuth 
        veli = 0.0        # Velocidade inicial 

        orbinit = np.array([loi, lati, radi, veli, fpai, chii])    
        simulation = integrate.odeint(rocket, orbinit, t, rtol=5e-10, atol=5e-15)
        
        range_x[:,i] = (rm*np.cos(simulation[:,1])*np.cos(simulation[:,0]))/1e3
        range_y[:,i] = (rm*np.cos(simulation[:,1])*np.sin(simulation[:,0]))/1e3
        range_z[:,i] = (rm*np.sin(simulation[:,1]))
        Range[:,i] = (rm*np.cos(simulation[:,1])*np.cos(simulation[:,0]))+ (rm*np.cos(simulation[:,1])*np.sin(simulation[:,0])) + (rm*np.sin(simulation[:,1]))
 
        lo[:,i] = simulation[:,0]
        lat[:,i] = simulation[:,1]
        rad[:,i] = simulation[:,2]
        vel[:,i] = simulation[:,3]
        fpa[:,i] = simulation[:,4]
        chi[:,i] = simulation[:,5]
    
    #img = mpimg.imread('C:\Users\orsur040\Pictures\map_ftbl.png')
    """   
    plt.figure(1, figsize = (60,60))
    #gca().set_aspect('equal', adjustable='box')
    plt.rcParams.update({'font.size': 15})
    plt.title('Longitude', fontsize = 25)
    plt.xlabel('Tempo [s]', fontsize = 20)
    plt.ylabel(r'$\lambda$ [deg]', fontsize = 20)
    plt.plot(t, lo[:,0]/dtr, label=r'$a_0 = 4.0g$')
    plt.plot(t, lo[:,1]/dtr, label=r'$a_0 = 4.5g$')
    plt.plot(t, lo[:,2]/dtr, label=r'$a_0 = 5.0g$')
    plt.legend()
    plt.grid(True)
    
    plt.figure(2, figsize = (60,60))
    #gca().set_aspect('equal', adjustable='box')
    plt.rcParams.update({'font.size': 15})
    plt.title('Latitude', fontsize = 25)
    plt.xlabel('Tempo [s]', fontsize = 20)
    plt.ylabel(r'$\delta$ [deg]', fontsize = 20)
    plt.plot(t, lat[:,0]/dtr,  label=r'$a_0 = 4.0g$')
    plt.plot(t, lat[:,1]/dtr,  label=r'$a_0 = 4.5g$')
    plt.plot(t, lat[:,2]/dtr,  label=r'$a_0 = 5.0g$')
    plt.legend()
    plt.grid(True)
    """
    plt.figure(3, figsize = (60,60))
    #plt.gca().set_aspect('equal', adjustable='box')
    plt.rcParams.update({'font.size': 15})
    plt.title('Altitude', fontsize = 25)
    plt.xlim(0, max(t) + 5)
    plt.ylim(0, max((rad[:,2] - rm)/1e3) + 5)
    plt.xlabel('Tempo [s]', fontsize = 20)
    plt.ylabel('Altitude [km]', fontsize = 20)
    plt.plot(t, (rad[:,0] - rm)/1e3,  label=r'$a_0 = 4.0g$',color='k',linestyle = ':')
    plt.plot(t, (rad[:,1] - rm)/1e3,  label=r'$a_0 = 4.5g$',color='k',linestyle = '-.')
    plt.plot(t, (rad[:,2] - rm)/1e3,  label=r'$a_0 = 5.0g$',color='k',linestyle = '--')
    plt.legend()
    plt.grid(True)
    
    plt.figure(4, figsize = (60,60))
    #gca().set_aspect('equal', adjustable='box')
    plt.rcParams.update({'font.size': 15})
    plt.title('Velocidade Resultante', fontsize = 25)
    plt.xlabel('Tempo [s]', fontsize = 20)
    plt.ylabel('Velocidade [m/s]', fontsize = 20)
    plt.plot(t, vel[:,0],  label=r'$a_0 = 4.0g$',color='k',linestyle = ':')
    plt.plot(t, vel[:,1],  label=r'$a_0 = 4.5g$',color='k',linestyle = '-.')
    plt.plot(t, vel[:,2],  label=r'$a_0 = 5.0g$',color='k',linestyle = '--')
    plt.legend()
    plt.grid(True)
    """
    plt.figure(5, figsize = (60,60))
    #gca().set_aspect('equal', adjustable='box')
    plt.rcParams.update({'font.size': 15})
    plt.title('Angulo de Voo', fontsize = 25)
    plt.xlabel('Tempo [s]',fontsize = 20)
    plt.ylabel(r'$\phi$ [deg]', fontsize = 20)
    plt.plot(t, fpa[:,0]/dtr,  label=r'$a_0 = 4.0g$')
    plt.plot(t, fpa[:,1]/dtr,  label=r'$a_0 = 4.5g$')
    plt.plot(t, fpa[:,2]/dtr,  label=r'$a_0 = 5.0g$')
    plt.legend()
    plt.grid(True)   
    
    plt.figure(6, figsize = (60,60))
    #gca().set_aspect('equal', adjustable='box')
    plt.rcParams.update({'font.size': 15})
    plt.title('Azimute', fontsize = 25)
    plt.xlabel('Tempo [s]', fontsize = 20)
    plt.ylabel(r'$A$ [deg]', fontsize = 20)
    plt.plot(t, chi[:,0]/dtr,  label=r'$a_0 = 4.0g$')
    plt.plot(t, chi[:,1]/dtr,  label=r'$a_0 = 4.5g$')
    plt.plot(t, chi[:,2]/dtr,  label=r'$a_0 = 5.0g$')
    plt.legend()
    plt.grid(True)
    """
    plt.figure(7, figsize = (60,60))
#    plt.imshow(img, imlim = ([-2.32, -2.25], [-44.37, -44.32]))
    plt.rcParams.update({'font.size': 15})
    plt.title('Ground Track', fontsize = 25)
    plt.xlabel(r'$\lambda$ [deg]', fontsize = 20)
    plt.ylabel(r'$\delta$ [deg]', fontsize = 20)
    plt.plot(lo[:,0]/dtr, lat[:,0]/dtr,  label=r'$a_0 = 4.0g$',color='k',linestyle = ':')
    plt.plot(lo[:,1]/dtr, lat[:,1]/dtr,  label=r'$a_0 = 4.5g$',color='k',linestyle = '-.')
    plt.plot(lo[:,2]/dtr, lat[:,2]/dtr,  label=r'$a_0 = 5.0g$',color='k',linestyle = '--')
    plt.legend()
    plt.grid(True)
    
    
    plt.figure(8, figsize = (60,60))
    #plt.gca().set_aspect('equal', adjustable='box')
    plt.rcParams.update({'font.size': 15})
    plt.title('Trajetoria', fontsize = 25)
    plt.xlabel('Alcance [km]', fontsize = 20)
    plt.ylabel('Altitude [km]', fontsize = 20)
    plt.ylim(0, max((rad[:,2] - rm)/1e3) + 5)
    plt.xlim(0, max(((Range[:,2] - min(Range[:,2])))/1e3) + 0.5)
    plt.plot((Range[:,0] - min(Range[:,0]))/1e3, (rad[:,0] - rm)/1e3,  label=r'$a_0 = 4.0g$',color='k',linestyle = ':')
    plt.plot((Range[:,1] - min(Range[:,1]))/1e3, (rad[:,1] - rm)/1e3,  label=r'$a_0 = 4.5g$',color='k',linestyle = '-.')
    plt.plot((Range[:,2] - min(Range[:,2]))/1e3, (rad[:,2] - rm)/1e3,  label=r'$a_0 = 5.0g$',color='k',linestyle = '--')
    plt.legend()
    plt.grid(True)
    
    mpl.rcParams['legend.fontsize'] = 15
    fig = plt.figure(9)
    ax = fig.gca(projection='3d')
    ax.plot((range_y[:,0] - min(range_y[:,0])), (range_x[:,0] - min(range_x[:,0])), (rad[:,0] - rm)/1e3, label='Trajetoria a = 4.0g',color='k',linestyle = ':')
    ax.plot((range_y[:,1] - min(range_y[:,1])), (range_x[:,1] - min(range_x[:,1])), (rad[:,1] - rm)/1e3, label='Trajetoria a = 4.5g',color='k',linestyle = '-.')
    ax.plot((range_y[:,2] - min(range_y[:,2])), (range_x[:,2] - min(range_x[:,2])), (rad[:,2] - rm)/1e3, label='Trajetoria a = 5.0g',color='k',linestyle = '--')
    ax.legend()

    return simulation 

if __name__ == "__main__":
    main()