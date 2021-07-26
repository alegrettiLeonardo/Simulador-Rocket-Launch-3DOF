# -*- coding: utf-8 -*-
"""
Created on Fri Mar 15 11:33:11 2019

@author: orsur040
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.pyplot import plot, grid, title, figure, xlabel, ylabel, xticks, legend, gca, yticks,  twinx
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
    g0 = 9.80665
    # Legend
    Legend = [str('UDMH|N2O4'), str('Ethanol|LOx'), str('Kerosene|LOx'), str('Methane|LOx')] 
    
    # Parameters
    Isp_1 = 200.0 # Specific Impulse First Stage Engine 
    Delta_v = 9.5e3 # Low Orit Earth 
    PL = 100.0 # Payload [kg]
    
    # Nondimensional speed change (delta_v/ve_1)
    vf = Delta_v/(9.81*Isp_1)  
    
    # Stage ratios 
    b = 1.5
    sigma = [0.07, 0.05]
    
    # Engine's
    fT_1 = 95.0e3
    fT_2 = 75.0e3
    
    # Storage matrices initialization 
    P = np.zeros([4,1000])
    A = np.zeros([4,1000]); mp = np.zeros([1000,4]) 
    aopt = np.zeros([1000,4]); popt = np.zeros([1000,4]);# m = np.zeros([3,1000]) 
    for j in np.arange(0, 4):
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
        
        # Rocket Mass
        m_1 = PL/mp.max()
        m_2 = (P.max())*m_1
        mp_1 = (m_1 - m_2)*(1 - sigma[0])
        mp_2 = (m_2 - PL)*(1 - sigma[1])
        
        ms_1 = m_1 - (PL + mp_1)
        ms_2 = m_2 - (PL + mp_2)
        
        # Rocket Engine
        m1_dot = fT_1/(g0*Isp_1)
        tb_1 = mp_1/m1_dot
        
        m2_dot = fT_2/(g0*Isp_1*b)
        tb_2 = mp_2/m2_dot
        
        '''
        #plt.fig1 = plt.figure(1, figsize = (10, 20))    
        plt.plot(A[j,:], P[j,:]*P[j,:]*A[j,:],'.', label = r'$\beta_2$')
        plt.title('Vatiation of total payload of a two-stage rocket', fontsize = 14)
        plt.xlabel(r'$\alpha_2$', fontsize = 20)
        plt.ylabel(r'$\lambda_T$', fontsize = 20)
        plt.legend()
        plt.grid(True)
        '''
        print'----------',
        print'\n Optimal Total Payload, Alpha and First-Stage Payload Ratios:',
        print'\n Propellants: ', Legend[j],
        print'\n Beta2 = ', b, 
        print'\n Lambda_T = %.5f'%mp.max(),
        print'\n Alpha_opt = %.5f'%aopt[:,j].max(),
        print'\n Lambda_opt = %.5f'%popt[:,j].max(),
        print'\n First Stage Mass = %.2f'%m_1, '[kg]',
        print'\n Second Stage Mass = %2.f'%m_2, '[kg]',
        print'\n First Stage Propellant Mass = %.2f'%mp_1, '[kg]',
        print'\n Second Stage Propellant Mass = %.2f'%mp_2, '[kg]',
        print'\n First Stage Dry Mass = %.2f'%ms_1, '[kg]',
        print'\n Second Stage Dry Mass = %.2f'%ms_2, '[kg]',
        print'\n Burn Time First Stage = %.2f'%tb_1, '[s]',
        print'\n Burn Time Second Stage = %.2f'%tb_2, '[kg]',
        print'\n Specific Impulse First Stage = %.2f'%Isp_1, '[s]',
        print'\n Specific Impulse Second Stage = %.2f'%(Isp_1*b), '[s]',
        print'\n----------',
        
        b = b + 0.15    
    
    plt.fig1 = plt.figure(1, figsize = (50, 50))    
    plt.rcParams.update({'font.size': 12})
    plt.plot(A[0,:], P[0,:]*P[0,:]*A[0,:],':' , color='k', label = 'Etanol|LOx')
    plt.plot(A[1,:], P[1,:]*P[1,:]*A[1,:],'-.', color='k', label = 'UDMH|N2O4')
    plt.plot(A[2,:], P[2,:]*P[2,:]*A[2,:],'--', color='k', label = 'Kerosene|LOx')
    plt.plot(A[3,:], P[3,:]*P[3,:]*A[3,:],'-' , color='k', label = 'Metano|LOx')
    plt.title('Foguete de Dois Estagios', fontsize = 25)
    plt.xlabel(r'$\alpha_2$', fontsize = 25)
    plt.ylabel(r'$\lambda_T$', fontsize = 25)
    plt.legend()
    legend(bbox_to_anchor=(1.005, 1), loc=2, borderaxespad=0.,fontsize=10)
    plt.grid(True)
    #print'\nmp = ', mp[:,0].max(), mp[:,1].max()   
    return 

if __name__ == "__main__":
    main()