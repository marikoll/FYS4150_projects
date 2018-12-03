#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec  3 10:33:10 2018

@author: maritkollstuen
"""

import numpy as np
import matplotlib.pyplot as plt

def assertion(init_psi, init_zeta,N_x):
    epsilon = 1e-10
    bc_0 = 0.0
    bc_N = 0.0
    if abs(init_psi[0] - bc_0) > epsilon and abs(init_psi[N_x] - bc_N) > epsilon:
        print('psi_0: ', init_psi[0], 'psi_N:', init_psi[N_x])
        print('Error, initial condition does not satisfy BC')
            
    psi_0 = np.empty(len(init_psi))
    zeta_0 = np.empty(len(init_zeta))
    
    
    for i in range(0,N_x):
        psi_0[i] = init_psi[i]
        zeta_0[i] = init_zeta[i]
    return psi_0, zeta_0, bc_0, bc_N

def tridiag(b, y, N, soltn):
    b[0] = 2.0
    
    #forward substitution
    for i in range(1, N):
         b[i] = (i + 2)/(i + 1)
         y[i] = y[i] + (y[i-1]/b[i-1])
         
    #backward substitution
    soltn[N-1] = y[N-1]/b[N-1]
    
    for i in range(N-2, 0, -1):
        soltn[i] = (y[i] + soltn[i+1])/b[i]
    
    return soltn

def euler(init_psi, init_zeta, N_x, dx, T, dt):
    psi_0, zeta_0,  bc_0, bc_N = assertion(init_psi, init_zeta, N_x)
    alpha = dt/(2*dx)
    
    psi_prev = np.zeros(N_x)
    psi_curr = np.zeros(N_x)
    zeta_prev = np.zeros(N_x)
    zeta_curr = np.zeros(N_x)
    
    
    diag = np.zeros(N_x-2)
    rhs_diag = np.zeros(N_x-2)

    for i in range(0,N_x):
        psi_prev[i]  = psi_0[i]
        zeta_prev[i] = zeta_0[i]

    psi_curr[0] = bc_0; psi_curr[N_x-1] = bc_N 
    zeta_curr[0] = zeta_prev[0]; zeta_curr[N_x-1] = zeta_prev[N_x-1]

    outstuff = np.zeros((N_x, int(float(T)/dt)))
    t = 0.0
    n = 0

    while t < T:
        #forward Euler:
        for i in range(1, N_x-1):
            zeta_curr[i] = zeta_prev[i] - alpha*(psi_prev[i+1] - psi_prev[i-1])
        for i in range(1, N_x-1):
            rhs_diag[i-1] = -dx**2*zeta_curr[i]
        psi_curr = tridiag(diag, rhs_diag, N_x -2, psi_curr)
        for i in range(1, N_x-1):
            psi_prev[i] = psi_curr[i]
            zeta_prev[i] = zeta_curr[i]
        
        t += dt
        if (n % 20 == 0):
            outstuff[:, n] = psi_curr[:]
        n += 1

    return psi_curr, outstuff   


def leapfrog(init_psi, init_zeta, N_x, dx, T, dt):
    psi_0, zeta_0,  bc_0, bc_N = assertion(init_psi, init_zeta, N_x)
    alpha = dt/(2.0*dx)
    gamma = dt/dx

    psi_prev = np.zeros(N_x)
    psi_curr = np.zeros(N_x)
    zeta_prev = np.zeros(N_x)
    zeta_pp = np.zeros(N_x)
    zeta_curr = np.zeros(N_x)
    
    diag = np.zeros(N_x-2)
    rhs_diag = np.zeros(N_x-2)
    
    for i in range(0,N_x):
        psi_prev[i]  = psi_0[i]
        zeta_pp[i] = zeta_0[i]
    
    psi_curr[0] = bc_0; psi_curr[N_x-1] = bc_N 
    zeta_curr[0] = zeta_pp[0]; zeta_curr[N_x-1] = zeta_pp[N_x-1]


    #initial Euler:
    for i in range(1, N_x-1):
        zeta_prev[i] = zeta_0[i] - alpha*(psi_0[i+1] - psi_0[i-1])
    for i in range(1, N_x-1):
        rhs_diag[i-1] = -dx**2*zeta_prev[i]

   
    psi_prev = tridiag(diag, rhs_diag, N_x-2, psi_prev)

    outstuff = np.zeros((N_x, int(float(T)/dt)))
    t = 0.0
    n = 0

    while t < T:
        for i in range(1, N_x-1):
            zeta_curr[i] = zeta_pp[i] - gamma*(psi_prev[i+1] - psi_prev[i-1])
        for i in range(1,N_x-1):
            rhs_diag[i-1] = -dx**2*zeta_curr[i]
        psi_curr = tridiag(diag, rhs_diag, N_x -2, psi_curr)
    
        for i in range(1, N_x -1):
            psi_prev[i] = psi_curr[i]
            zeta_pp[i] = zeta_prev[i]
            zeta_prev[i] = zeta_curr[i]
        t += dt
        if (n % 20 == 0):
            outstuff[:, n] = psi_curr[:]
        n += 1

    return psi_curr, outstuff   


if __name__ == "__main__":
    N = 40
    T = 150
    
    dx = 1.0/40
    dt = 0.001
    
    init_psi = np.zeros(N)
    init_zeta = np.zeros(N)
    
    init_psi_gauss = np.zeros(N)
    init_zeta_gauss = np.zeros(N)  
    sigma = 0.1
    
    for i in range(0, N-1):
        x = i*dx
        init_psi[i] = np.sin(4.0*np.pi*x)
        init_zeta[i] = -16.0*np.pi**2*np.sin(4.0*np.pi*x)
        
#        init_psi_gauss[i] = np.exp(-((x-0.5)/sigma)**2)
#        init_zeta_gauss[i] = (4*((x-0.5)/sigma)**2) - (2/sigma**2)*(np.exp(-((x-0.5)/sigma)**2))
#    

    psiE, outstuff= euler(init_psi, init_zeta, N, dx, T, dt)
    psiLF, outstuff2 = leapfrog(init_psi, init_zeta, N, dx, T, dt)
    
    psiE_gauss = euler(init_psi_gauss, init_zeta_gauss, N, dx, T, dt)
    psiLF_gauss = leapfrog(init_psi_gauss, init_zeta_gauss, N, dx, T, dt)
    
    x = np.linspace(0, 1, N)
    
    plt.figure()
    plt.plot(x, outstuff[:,0], 'r-')
    plt.plot(x, outstuff2[:,0], 'b-.')

    
#    plt.figure()
#    plt.plot(x, psiE_gauss[1:N-3], 'r-')
#    plt.plot(x, psiLF_gauss[1:N-3], 'b-.')




