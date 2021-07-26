# -*- coding: utf-8 -*-
"""
Created on Mon Mar 18 14:42:24 2019

@author: orsur040
"""
import numpy as np

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
