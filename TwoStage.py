# -*- coding: utf-8 -*-
"""
Created on Mon Feb 25 09:40:51 2019

@author: orsur040
"""

import numpy as np
import matplotlib.pyplot as plt

def NStage(vf, beta, epsilon, alpha):
    
    # Program for solving the multi-stage (N) rocket equation for the first-satge
    # Payload ratio (p)
    # Beta = (Nx1) vector of ratios of specific impulses to that of the first stage
    # Beta(k) = Isp_k/Isp_1, first element of beta should be 1.0
    # Alpha = (Nx1) vector of ratios of payload ratios to that of the first stage
    # Alpha(k) = lambda_k/lambda_1, first element of alpha should be 1.0
    # Epsilon = (Nx1) vector of structural ratios of the stages
    # Vf = ratio of total velocity impulse to exhaust speed of first stage
    # vf = Delta_v/v_e1
    
    N = np.size(beta)
    p = 0.1
    f = vf
    tol = 1e-9
    for k in np.arange(0, N):
        f = f + beta[k]*np.log(epsilon[k] + alpha[k]*((1 - epsilon[k])*p))
    
    while np.abs(f) > tol:
        f = vf
        fp = 0.0
        for k in np.arange(0, N):
            f =  f + beta[k]*np.log(epsilon[k] + alpha[k]*((1 - epsilon[k])*p))
            fp = fp + alpha[k]*beta[k]/(epsilon[k] + alpha[k]*(1 - epsilon[k])*p)
        d = -f/fp
        p = p + d
    
    return p

def main():
    
    # Parameters
    Isp_1 = 200.0; Delta_v = 9.5e3 
    # Nondimensional speed change (delta_v/ve_1)
    vf = Delta_v/(9.81*Isp_1)  
    
    # Stage ratios 
    b = 1.5
    sigma = [0.07, 0.05]
    
    # Storage matrices initialization 
    P = np.zeros([2,1000])
    A = np.zeros([2,1000]); mp = np.zeros([1000,2]) 
    aopt = np.zeros([1000,2]); popt = np.zeros([1000,2]) 
    for j in np.arange(0, 2):
        a = 0.05    # Alpha initial value 
        for k in np.arange(0, 1000):  # Alpha interation 
            # Newton interation for first-satge payload ratio
            p = NStage(vf, [1.0, b], sigma, [1.0, a])
            P[j,k] = p
            A[j,k] = a
            a = a + 0.0005
              
        # Optimal total payload, alpha and fist-stage payload ratios 
        mp[:,j] = (P[j,:]*P[j,:]*A[j,:])
        kidx = np.argmax(mp[:,j])
        aopt[:,j] = A[j,kidx]
        popt[:,j] = P[j,kidx]
        
        #plt.fig1 = plt.figure(1, figsize = (10, 20))    
        plt.plot(A[j,:], P[j,:]*P[j,:]*A[j,:],'.', label = r'$\beta_2$')
        plt.title('Vatiation of total payload of a two-stage rocket', fontsize = 14)
        plt.xlabel(r'$\alpha_2$', fontsize = 20)
        plt.ylabel(r'$\lambda_T$', fontsize = 20)
        plt.legend()
        plt.grid(True)
        
        print '----------',
        print'\n Optimal Total Payload, Alpha and First-Stage Payload Ratios:',
        print'\n Beta2 = ', b, 
        print'\n mp = %.5f'%mp.max(),
        print'\n aopt = %.5f'%aopt[:,j].max(),
        print'\n popt = %.5f'%popt[:,j].max()
        b = b + 0.25      
    #print'\nmp = ', mp[:,0].max(), mp[:,1].max()   
    return 

if __name__ == "__main__":
    main()