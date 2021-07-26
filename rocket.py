# -*- coding: utf-8 -*-
"""
Created on Mon Mar 18 11:56:20 2019

@author: orsur040
"""

import numpy as np 







def rocket(t,o):
    global dtr; global mu; global omega; global S; global c; global rm; 
    global tb1; global tb2; global fT1; global fT2; global m01; 
    global m02; global mL; global mp1; global mp2; global Gamma; 
    global f8; 
    
    # Acceleration due to gravity (nonspherical earth)
    [g, gn] = gravity(o[2], o[1]); 
    lo = o[0];la = o[1];
    clo = np.cos(lo); slo = np.sin(lo); cla = cos(la); sla = sin(la);
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
        CD = 8*CD;
    elif t <= (tb1+tb2):
        fT = fT2; 
        m = m02 - mp2*(t - tb1)/tb2; 
        CD = 3*CD;
    else:
        fT = 0.0;
        m = mL;
#[t alt m mach]
    D = Qinf*S*CD;
    Xfo = fT-D; Yfo = 0; Zfo = 0;
    longidot = o(4)*cfpa*schi/(o(3)*cla); 
    latidot =  o(4)*cfpa*cchi/o(3); 
    raddot = o(4)*sfpa;  
    veldot = -g*sfpa + gn*cchi*cfpa + Xfo/m + omega*omega*o[2]*cla*(sfpa*cla - cfpa*cchi*sla);
    if t <= 10.0:
        headdot = 0.0; 
        gammadot = 0.0;
    else:
        gammadot = (o[3]/o[2] -g/o[3])*cfpa - gn*cchi*sfpa/o[3] + Zfo/(o[3]*m) + 2.0*omega*schi*cla + omega*omega*o[2]*cla*(cfpa*cla + sfpa*cchi*sla)/o[3];
        if abs(cfpa) > 1e-6:
            headdot = o[3]*schi*np.tan(o[1])*cfpa/o[2] - gn*schi/o[3] - Yfo/(o[3]*cfpa*m) - 2.0*omega*(np.tan(o[4])*cchi*cla - sla) + omega*omega*o[2]*schi*sla*cla/(o[3]*cfpa);
        else:
            headdot = 0.0;

    deriv = [longidot, latidot, raddot, veldot, gammadot, headdot];
    if alt <= 120.0e3:
        Qdot=Qinf*v*S*CD/20;
        print(f8,'\t%1.5e\t%1.5e\t%1.5e\t%1.5e\t%1.5e\t%1.5e\t%1.5e\n',t,alt,m,mach,veldot,Qinf,Qdot);
    
    return deriv