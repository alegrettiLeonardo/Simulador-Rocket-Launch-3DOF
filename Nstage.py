# -*- coding: utf-8 -*-
"""
Created on Fri Feb 22 11:51:07 2019

@author: orsur040
"""
import numpy as np

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
    tol = 1.0e-12
    for k in np.arange(0, N):
        f = f + beta[k]*np.log(epsilon[k] + alpha[k]*((1 - epsilon[k])*p))
    print'Payload Ratio:'
    while np.abs(f) > tol:
        f = vf
        fp = 0.0
        for k in np.arange(0, N):
            f =  f + beta[k]*np.log(epsilon[k] + alpha[k]*((1 - epsilon[k])*p))
            fp = fp + alpha[k]*beta[k]/(epsilon[k] + alpha[k]*(1 - epsilon[k])*p)
        d = -f/fp
        p = p + d
        print 'p = ',p
    
    return p

def main():
    # Stage ratios 
    beta = [1.0, 340/220.0]
    sigma = [0.07, 0.05]
    alpha = [1.0, 1.2]
    N = np.size(beta)
    Delta_v = 9.5e3 
    PL = 100 # Payload [kg]
    # Nondimensional speed change (delta_v/ve_1)
    vf = Delta_v/(9.81*220.0)
    print 'Number of Stage:', N
    
    # Newton interation for first-satge payload ratio
    p = NStage(vf, beta, sigma, alpha)
    
    # Total payload ratio
    pT = ((p)**N)*(alpha[1])#*alpha[2]
    
    # Rocket Mass
    m_1 = PL/pT
    m_2 = p*m_1
    #m_3 = p*alpha[1]*m_2
    #mp_3 = (m_3 - PL)*(1 - sigma[2])
    mp_1 = (m_1 - m_2)*(1 - sigma[0])
    mp_2 = (m_2 - PL)*(1 - sigma[1])
    ms_1 = m_1 - (PL + mp_1)
    ms_2 = m_2 - (PL + mp_2)
    
    print '\n'
    print '--------------------------'
    print 'Total Mass_1º = %.4f'%m_1, '[kg]'
    print 'Total Mass_2º = %.4f'%m_2, '[kg]'
    print '--------------------------'
    #print 'Mass_3 = %.4f'%m_3, '[kg]'
    print 'Propellant_1º = %.4f'%mp_1, '[kg]'
    print 'Propellant_2º = %.4f'%mp_2, '[kg]'
    #print 'Propellant_3 = %.4f'%mp_3, '[kg]'
    print '--------------------------'
    print 'Dry Mass_1º = %.4f'%ms_1, '[kg]'
    print 'Dry Mass_2º = %.4f'%ms_2, '[kg]'
    print '--------------------------'
    return p

if __name__ == "__main__":
    main()