# -*- coding: utf-8 -*-
"""
Created on Mon Mar 18 11:56:20 2019

@author: orsur040
"""

import numpy as np 
from scipy.integrate import odeint

def atmosphere(h, vel, CL):
    R = 287; 
    go = 9.806; 
    Na = 6.0220978e23; 
    sigma = 3.65e-10; 
    S = 110.4;  
    Mo = 28.964; 
    To = 288.15; 
    Po = 1.01325e5; 
    re = 6378.14e3; 
    Beta = 1.458e-6; 
    gamma = 1.405; 
    B = 2/re;
    layers = 21;
    P = np.zeros([]); rho = np.zeros([])
    Z = 1e3*[0.00, 11.0191, 20.0631, 32.1619, 47.3501, 51.4125, 71.8020, 86.00, 100.00, 110.00, 120.00, 150.00, 160.00, 170.00, 190.00, 230.00, 300.00, 400.00, 500.00, 600.00, 700.00, 2000.00];
    T = [To, 216.65, 216.65, 228.65, 270.65, 270.65, 214.65, 186.946, 210.65, 260.65, 360.65, 960.65, 1110.60, 1210.65, 1350.65, 1550.65, 1830.65, 2160.65, 2420.65, 2590.65, 2700.00, 2700.0];
    M = [Mo, 28.964, 28.964, 28.964, 28.964, 28.964, 28.962, 28.962, 28.880, 28.560, 28.070, 26.920, 26.660, 26.500, 25.850, 24.690, 22.660, 19.940, 17.940, 16.840, 16.170, 16.17];
    LR = [-6.5e-3, 0.0, 1e-3, 2.8e-3, 0.0, -2.8e-3, -2e-3, 1.693e-3, 5.00e-3, 1.0e-2, 2.0e-2, 1.5e-2, 1.0e-2, 7.0e-3, 5.0e-3, 4.0e-3, 3.3e-3, 2.6e-3, 1.7e-3, 1.1e-3, 0.0];
    
    rho0 = Po/(R*To);
    P[0] = Po;
    T[0] = To;
    rho[0] = rho0;
    for i in range(0, layers):
        if ~(LR[i] == 0):
            C1 = 1 + B*( T[i]/LR[i] - Z[i]);
            C2 = C1*go/(R*LR[i]);
            C3 = T[i+1]/T[i];
            C4 = C3**(-C2);
            C5 = np.exp(go*B*(Z[i+1] - Z[i])/(R*LR[i]));
            P[i+1] = P[i]*C4*C5;
            C7 = C2 + 1;
            rho[i+1] = rho[i]*C5*(C3**(-C7));
        else:
            C8 = -go*(Z[i+1] - Z[i])*(1 - B*(Z[i+1] + Z[i])/2.0)/(R*T[i]);
            P[i+1] = P[i]*np.exp(C8);
            rho[i+1] = rho[i]*np.exp(C8);
        for i in range (0,21):
            if h < Z[i+1]:
                if ~(LR(i) == 0):
                    C1 = 1 + B*( T[i]/LR[i] - Z[i] );
                    TM = T[i] + LR[i]*(h - Z[i]);
                    C2 = C1*go/(R*LR[i]);
                    C3 = TM/T[i];
                    C4 = C3**(-C2);
                    C5 = np.exp(B*go*(h - Z[i])/(R*LR[i]));
                    PR = P[i]*C4*C5;
                    C7 = C2 + 1;
                    rhoE = C5*np.rho[i]*(C3**(-C7));
                else:
                    TM = T[i];
                    C8 = -go*(h - Z[i])*(1 - (h + Z[i])*B/2.0)/(R*T[i]);
                    PR = P[i]*np.exp(C8);
                    rhoE = rho[i]*np.exp(C8);     
                MOL = M(i) + ( M(i+1)-M(i) )*( h - Z(i) )/( Z(i+1) - Z(i) );
                TM = MOL*TM/Mo;
                asound = np.sqrt(gamma*R*TM); 
                MU = Beta*TM^1.5/(TM + S); 
                KT = 2.64638e-3*TM^1.5/(TM + 245.4*10^(-12/TM));
                Vm = np.sqrt(8*R*TM/np.pi);
                m = MOL*1e-3/Na;
                n = rhoE/m;
                F = np.sqrt(2)*np.pi*n*sigma^2*Vm;
                L = Vm/F; 
                Mach = vel/asound; 
                T0 = TM*(1 + (gamma - 1)*Mach^2/2);
                MU0 = Beta*T0^1.5/(T0 + S);
                RE0 = rhoE*vel*CL/MU0;
                RE = rhoE*vel*CL/MU; 
                Kn = L/CL; 
                Kno = 1.25*np.sqrt(gamma)*Mach/RE0;
                if Kn >= 10.0:
                    d = 1;
                elif Kn <= 0.01: 
                    d = 2;
                else:         
                    d = 3;
                Y = [TM, rhoE, Mach, Kn, asound, d, PR, MU, RE];
        return Y

def conticap(mach):
    machr = [0.0, 0.2, 0.3, 0.4, 0.5, 0.6, 0.8, 0.9, 0.95, 1.05, 1.1, 1.2, 1.6, 2.0, 2.5, 3.0, 3.8, 5.0, 10.0, 99.0];
    Cdr = [.475475, .475475, .47576, .48336, .488965, .508345, .56563, .618165, .668135,\
       1.031795, 1.01707, .990565, .815955, .69236, .60971, .54606, .513, .494, .48317, .48317];
   
    Cd = scipy.interp1(machr, Cdr, mach);
    return Cd

def gravity(r,lat):
    phi = np.pi/2.0 - lat;
    mu = 3.986004e14;
    Re = 6378.135e3;
    J2 = 1.08263e-3;
    J3 = 2.532153e-7;
    J4 = 1.6109876e-7;
    
    gc = mu*(1.0 - 1.5*J2*(3.0*np.cos(phi)**2 - 1)*(Re/r)**2 - 2.0*J3*np.cos(phi)\
    *(5.0*np.cos(phi)**2 - 3.0)*(Re/r)**3 - (5.0/8.0)*J4*(35.0*np.cos(phi)**4 - 30.0*np.cos(phi)**2 + 3.0)*(Re/r)**4)/r**2
   
    gnorth = -3.0*mu*np.sin(phi)*np.cos(phi)*(Re/r)*(Re/r)*(J2 + 0.5*J3*(5.0*np.cos(phi)**2 - 1)\
    *(Re/r)/np.cos(phi) + (5.0/6.0)*J4*(7.0*np.cos(phi)**2 - 1.0)*(Re/r)**2)/r**2
    
    g = [gc, gnorth]
    return g

def rocket(t,o):
    global dtr; global mu; global omega; global S; global c; global rm; 
    global tb1; global tb2; global fT1; global fT2; global m01; 
    global m02; global mL; global mp1; global mp2; global Gamma; 
    global f8; 
    
    # Acceleration due to gravity (nonspherical earth)
    [g, gn] = gravity(o[2], o[1]); 
    lo = o[0];la = o[1];
    clo = np.cos(lo); slo = np.sin(lo); cla = cos(la); sla = np.sin(la);
    fpa = o[4]; chi = o[5];
    cfpa = np.cos(fpa); sfpa = np.sin(fpa); cchi = np.cos(chi); schi = np.sin(chi); 
    if o[2] < rm:
        o[2] = rm;

    alt = o[2] - rm;
    v  = o[3];
    if v < 0.:
        v = 0.0;

    if alt <= 2000.0e3:
        atmosp = atmosphere(alt,v,c);
        rho = atmosp[1];     
        Qinf = 0.5*rho*(v**2);
        mach = atmosp[2]; 
        Kn = atmosp[3];
        CDC = conticap(mach); 
        s = mach*np.sqrt(Gamma/2.0);
        CDFM = 1.75 + np.sqrt(np.pi)/(2.0*s); 
        iflow = atmosp[5];
        if iflow == 2.0:
            CD = CDC;
        elif iflow == 1.0:
            CD = CDFM;
        else:
            CD = CDC + (CDFM - CDC)*(0.333*np.log10(Kn/np.sin(np.pi/6.0)) + 0.5113);
    else:
        rho = 0.0;
        Qinf = 0.0;
        CD = 0.0;
        mach = 0.0;
    if t<= tb1: 
        fT = fT1; 
        m = m01 - mp1*t/tb1; 
        CD = 8.0*CD;
    elif t <= (tb1+tb2):
        fT = fT2; 
        m = m02 - mp2*(t - tb1)/tb2; 
        CD = 3.0*CD;
    else:
        fT = 0.0;
        m = mL;
#[t alt m mach]
    D = Qinf*S*CD;
    Xfo = fT - D; Yfo = 0.0; Zfo = 0.0;
    longidot = o[3]*cfpa*schi/(o[2]*cla); 
    latidot =  o[3]*cfpa*cchi/o[2]; 
    raddot = o[3]*sfpa;  
    veldot = -g*sfpa + gn*cchi*cfpa + Xfo/m + omega*omega*o[2]*cla*(sfpa*cla - cfpa*cchi*sla);
    if t <= 10.0:
        headdot = 0.0; 
        gammadot = 0.0;
    else:
        gammadot = (o[3]/o[2] -g/o[3])*cfpa - gn*cchi*sfpa/o[3] + Zfo/(o[3]*m) + 2.0*omega*schi*cla + omega*omega*o[2]*cla*(cfpa*cla + sfpa*cchi*sla)/o[3];
        if abs(cfpa) > 1.0e-6:
            headdot = o[3]*schi*np.tan(o[1])*cfpa/o[2] - gn*schi/o[3] - Yfo/(o[3]*cfpa*m) - 2.0*omega*(np.tan(o[4])*cchi*cla - sla) + omega*omega*o[2]*schi*sla*cla/(o[3]*cfpa);
        else:
            headdot = 0.0;

    deriv = [longidot, latidot, raddot, veldot, gammadot, headdot];
    if alt <= 120.0e3:
        Qdot = Qinf*v*S*CD/20.0;
        print(f8,'\t%1.5e\t%1.5e\t%1.5e\t%1.5e\t%1.5e\t%1.5e\t%1.5e\n',t,alt,m,mach,veldot,Qinf,Qdot);
    
    return deriv

def main():
    global dtr; dtr = np.pi/180;
    global mu; mu = 3.986004e14;
    global S; S = 4.0;
    global c; c = 0.5;
    global m; m = 350;
    global rm; rm = 6378140;
    global omega; omega = 2*np.pi/(23*3600 + 56*60 + 4.0905);
    global Gamma; Gamma = 1.41; 
    global f8; f8 = open('data8.dat', 'w');
    global tb1;
    global tb2; 
    global fT1; 
    global fT2; 
    global m01; 
    global m02; 
    global mL; 
    global mp1; 
    global mp2;
   
    
    longi = -10*dtr;            
    lat = -79.8489182889*dtr; 
    rad = 6579.89967e3;  
    vel = 7589.30433867; 
    fpa = 0.54681217*dtr;
    chi = 99.955734*dtr; 
    options = 1.0e-8;    
    orbinit = [longi, lat, rad, vel]; 
    [t, o] = odeint(rocket,[0, 1750.0],orbinit,args=(options,));
    
    return

if __name__ == "__main__":
    main()