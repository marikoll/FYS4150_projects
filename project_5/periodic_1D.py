#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec  3 15:22:26 2018

@author: maritkollstuen
"""

import numpy as np
import matplotlib.pyplot as plt

def assertion(init_psi, init_zeta,N_x):
    epsilon = 1e-10
 
    if abs(init_psi[0] - init_psi[N_x -1]) > epsilon and abs(init_zeta[0] - init_zeta[N_x-1]) > epsilon:
        print('psi_0: ', init_psi[0], 'psi_N:', init_psi[N_x])
        print('Error, initial condition does not satisfy BC')
            
    psi_0 = np.empty(len(init_psi))
    zeta_0 = np.empty(len(init_zeta))
    

    psi_0 = init_psi
    zeta_0 = init_zeta
    return psi_0, zeta_0

def periodic_matrix(n_rows, n_cols):
    A = np.zeros((n_rows, n_cols))
    for i in range(0,n_rows):
        for j in range(0,n_cols):
            if i == j:
                A[i, j] = 2.0
            elif abs(i-j) ==1:
                A[i,j] = -1.0
    A[0, n_cols-1] = -1.0
    A[n_rows -1, 0] = -1.0
    
    return A


def euler(init_psi, init_zeta, N_x, dx, T, dt):
    psi_0, zeta_0= assertion(init_psi, init_zeta, N_x)
    alpha = dt/(2*dx)
    
    psi_prev = np.zeros(N_x)
    psi_curr = np.zeros(N_x)
    zeta_prev = np.zeros(N_x)
    zeta_curr = np.zeros(N_x)

    
    rhs_poisson = np.zeros(N_x-1)
    A = periodic_matrix(int(N_x-1),int(N_x-1))
    
    
    psi_prev  = psi_0
    zeta_prev = zeta_0
    
    outstuff = np.zeros((N_x-1, int(float(T)/dt)+1))
    t = 0.0
    n = 0

    while t < T:
        #forward Euler:
        for i in range(1, N_x-1):
            zeta_curr[i] = zeta_prev[i] - alpha*(psi_prev[i+1] - psi_prev[i-1])
            
        zeta_curr[0] = zeta_prev[0] - alpha*(psi_prev[1] - psi_prev[N_x -2])
        zeta_curr[N_x-1] = zeta_curr[0]

        for i in range(0, N_x-1):
            rhs_poisson[i-1] = -dx**2*zeta_curr[i]


        psi_curr = np.linalg.solve(A, rhs_poisson)
        
        psi_curr[N_x-2] = psi_curr[0]
        
        for i in range(0, N_x-1):
            psi_prev[i] = psi_curr[i]
            zeta_prev[i] = zeta_curr[i]
        
        t += dt
        if (n % 20 == 0):
            outstuff[:, n] = psi_curr[:]
 
        n += 1

    return outstuff   
        
def leapfrog(init_psi, init_zeta, N_x, dx, T, dt):
    psi_0, zeta_0 = assertion(init_psi, init_zeta, N_x)

    alpha = dt/(2*dx)
    gamma =  dt/dx
    
    psi_prev = np.zeros(N_x)
    psi_curr = np.zeros(N_x)
    zeta_prev = np.zeros(N_x)
    zeta_pp = np.zeros(N_x)
    zeta_curr = np.zeros(N_x)

    
    rhs_poisson = np.zeros(N_x-1)
    A = periodic_matrix(int(N_x-1),int(N_x-1))
    
    psi_prev  = psi_0
    zeta_prev = zeta_0
    
    #initial Euler:
    for i in range(1, N_x-1):
        zeta_prev[i] = zeta_0[i] - alpha*(psi_0[i+1] - psi_0[i-1])

    zeta_prev[0] = zeta_0[0] - alpha*(psi_0[1] - psi_0[N_x -2])
    zeta_prev[N_x-1] = zeta_prev[0]    
    
    for i in range(0, N_x-1):
        rhs_poisson[i-1] = -dx**2*zeta_prev[i]

    psi_prev = np.linalg.solve(A, rhs_poisson)
    psi_prev[N_x-2] = psi_prev[0]
    
    
    outstuff = np.zeros((N_x-1, int(float(T)/dt)+1))
    t = 0.0
    n = 0
    
    while t < T:
        #forward Euler:
        for i in range(1, N_x-2):
            zeta_curr[i] = zeta_pp[i] - gamma*(psi_prev[i+1] - psi_prev[i-1])
        
        zeta_curr[0] = zeta_pp[0] - gamma*(psi_prev[1] - psi_prev[N_x -2])
        zeta_curr[N_x-1] = zeta_curr[0]
        
        for i in range(1, N_x-1):
            rhs_poisson[i] = -dx**2*zeta_curr[i]
        
        psi_curr = np.linalg.solve(A, rhs_poisson)
        
        psi_curr[N_x-2] = psi_curr[0]
        
        
        for i in range(0, N_x-1):
            psi_prev[i] = psi_curr[i]
            zeta_pp[i] = zeta_prev[i]
            zeta_prev[i] = zeta_curr[i]
        
        t += dt
        if (n % 20 == 0):
            outstuff[:, n] = psi_curr[:]

        n += 1

    return outstuff   
        
        
if __name__ == "__main__":

    T = 150
    dt = 0.5
    
    dx = 1.0/40
    L = 1.0
    N = int(L/dx + 1)
    
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

    outstuff= euler(init_psi, init_zeta, N, dx, T, dt)
    outstuff2 = leapfrog(init_psi, init_zeta, N, dx, T, dt)
    
#    psiE_gauss = euler(init_psi_gauss, init_zeta_gauss, N, dx, T, dt)
#    psiLF_gauss = leapfrog(init_psi_gauss, init_zeta_gauss, N, dx, T, dt)
    
    x = np.linspace(0, 1, N-1)
    
    plt.figure()
    plt.plot(x, outstuff[:,0], 'r-')
    plt.plot(x, outstuff2[:,0], 'b-.')

    
#    plt.figure()
#    plt.plot(x, psiE_gauss[1:N-3], 'r-')
#    plt.plot(x, psiLF_gauss[1:N-3], 'b-.')       
        
        
        
        
        
        
        

    
    