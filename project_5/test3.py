#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec  5 10:22:06 2018

@author: maritkollstuen
"""

import numpy as np
import matplotlib.pyplot as plt


def initialize(N, dx):
    init_psi = np.zeros(N)
    init_zeta = np.zeros(N)
#    init_psi2 = np.zeros(N)
#    init_zeta2 = np.zeros(N)
    
#    init_psi_gauss = np.zeros(N)
#    init_zeta_gauss = np.zeros(N)  
#    sigma = 0.1
    
    for i in range(0, N-1):
        x = i*dx
        init_psi[i] = np.sin(4.0*np.pi*x) 
        init_zeta[i] = -16.0*np.pi**2*np.sin(4.0*np.pi*x)
#        init_psi2[i] = np.sin(4.0*np.pi*x)
##        print(init_psi2[i], init_psi[i])          
#        init_zeta2[i] = -16.0*np.pi**2*np.sin(4.0*np.pi*x)
        
#        init_psi_gauss[i] = np.exp(-((x-0.5)/sigma)**2)
#        init_zeta_gauss[i] = (4*((x-0.5)/sigma)**2) - (2/sigma**2)*(np.exp(-((x-0.5)/sigma)**2))
#   
    return init_psi, init_zeta

def tridiagonal_solidb(psi, force, dx, N):
    e1 = np.ones(N)
    e2 = np.ones(N)
    d  = np.ones(N)*(-2)
    
    psi[-1] = 0.0
    psi[0] = 0.0
    
    for i in range(2, N):
        d[i] -=(e1[i]*e2[i-1]/d[i-1])
        force[i] -= ((e1[i])*(force[i-1]))/d[i-1]
    for i in range(N-2, 1, -1):
        psi[i] = (force[i] - (e2[i]*(psi[i+1])))/d[i]
    return psi


def forward_solidb(psi, vorticity, N, alpha):
    psi[0] = 0.0
    psi[-1] = 0.0
    
    for i in range(1, N-1):
        vorticity[i] = vorticity[i] + (psi[i-1] - psi[i +1])*alpha
        
    return vorticity

def center_solidb(psi, vorticity, first_vorticity, N, gamma, alpha, time):
    
    if time == 0:
        first_vorticity = vorticity
        
        for i in range(1, N-1):
            vorticity[i] = vorticity[i] + (psi[i-1] - psi[i +1])*alpha
    if time != 0:
        for i in range(1, N-1):
            temp_first_vorticity = first_vorticity[i]
            temp_psi1 = psi[i +1]
            temp_psi2 = psi[i-1]
            
            first_vorticity = vorticity
            
            vorticity[i] = temp_first_vorticity + (temp_psi2- temp_psi1)*gamma
            
    return vorticity


if __name__ == "__main__":
    T = 150
    dt = 0.001

    dx = 1.0/40
    L = 1.0
    N = int(L/dx + 1)
    
    alpha = dt/(2*dx)
    gamma = dt/dx
    
    
    init_psi, init_vorticity = initialize(N, dx)
    first_vorticity = init_vorticity
    time = 0.0
    
    force_fwd = np.zeros(N)
    force_center = np.zeros(N)
    
    while time < T:
        v_fwd = forward_solidb(init_psi, first_vorticity, N, alpha)
        #v_center = center_solidb(init_psi, init_vorticity, first_vorticity, N, gamma, alpha, time)
    
            
        for i in range(0, N-1):
            force_fwd[i] = v_fwd[i]*dx**2
        #    force_center[i] = v_center[i]*dx**2
        psi_fwd = tridiagonal_solidb(init_psi, force_fwd, dx, N)
        #psi_center = tridiagonal_solidb(init_psi, force_center, dx, N)
    
        time +=dt
    
    x = np.linspace(0,1,N)
    plt.figure()
    plt.plot(x, psi_fwd, 'b-.')
    
    
    
    
    
    
    
        