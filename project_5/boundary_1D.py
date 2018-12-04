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
    
    
    psi_0 = init_psi
    zeta_0 = init_zeta
    
    return psi_0, zeta_0, bc_0, bc_N

def tridiag(b, y, N, soltn):
    
    b[0] = 2.0
    
    #forward substitution
    for i in range(1, N):
         b[i] = float(i + 2)/float(i + 1)
         y[i] = y[i] + (float(y[i-1])/float(b[i-1]))
         
    #backward substitution
    soltn[N-1] = float(y[N-1])/float(b[N-1])
    
    for i in range(N-2, 0, -1):
        soltn[i] = float(y[i] + soltn[i+1])/float(b[i])
    
    return soltn

def euler(N_x, dx, T, dt):
#    psi_0, zeta_0,  bc_0, bc_N = assertion(init_psi, init_zeta, N_x)
    psi_0, zeta_0 = initialize(N_x, dx)
    alpha = dt/(2.0*dx)
    dx2 = dx**2
    bc_0 = 0.0
    bc_N = 0.0
    
    psi_prev = np.zeros(N_x)
    psi_curr = np.zeros(N_x)
    zeta_prev = np.zeros(N_x)
    zeta_curr = np.zeros(N_x)
    
    
    # Arrays for tridiagonal solver
    diag = np.ones(N_x-2)*(-1)
    rhs_diag = np.zeros(N_x-2)


    # initial condition and boundary conditions
    psi_prev  = psi_0
    zeta_prev = zeta_0
    
    psi_curr[0] = bc_0; psi_curr[N_x-1] = bc_N 
    zeta_curr[0] = zeta_prev[0]; zeta_curr[N_x-1] = zeta_prev[N_x-1]

    outresults = np.zeros((N_x, int(float(T)/dt)))
    t = 0.0
    n = 0

    while t < T:
        #forward Euler:
        for i in range(1, N_x-1):
            zeta_curr[i] = zeta_prev[i] - alpha*(psi_prev[i+1] - psi_prev[i-1])
        for i in range(1, N_x-1):
            rhs_diag[i-1] = -dx2*zeta_curr[i]
        psi_curr = tridiag(diag, rhs_diag, N_x -2, psi_curr)
        for i in range(1, N_x-1):
            psi_prev[i] = psi_curr[i]
            zeta_prev[i] = zeta_curr[i]
        
        t += dt
        if (n % 20 == 0):
            outresults[:, n] = psi_curr[:]
        n += 1

    return outresults   


def leapfrog(N_x, dx, T, dt):
    psi_0, zeta_0 = initialize(N_x, dx)
    alpha = dt/(2*dx)
    dx2 = dx**2
    bc_0 = 0.0
    bc_N = 0.0
    gamma = dt/dx

    psi_prev = np.zeros(N_x)
    psi_curr = np.zeros(N_x)
    zeta_prev = np.zeros(N_x)
    zeta_pp = np.zeros(N_x)
    zeta_curr = np.zeros(N_x)
    
    
    # Arrays for tridiagonal solver
    diag = np.ones(N_x-2)*(-1)
    rhs_diag = np.zeros(N_x-2)
    
    # initial condition and boundary conditions
    psi_prev  = psi_0
    zeta_pp = zeta_0

    
    psi_curr[0] = bc_0; psi_curr[N_x-1] = bc_N 
    zeta_curr[0] = zeta_pp[0]; zeta_curr[N_x-1] = zeta_pp[N_x-1]


    #initial Euler:
    for i in range(1, N_x-1):
        zeta_prev[i] = zeta_0[i] - alpha*(psi_0[i+1] - psi_0[i-1])
    for i in range(1, N_x-1):
        rhs_diag[i-1] = -dx2*zeta_prev[i]

   
    psi_prev = tridiag(diag, rhs_diag, N_x-2, psi_prev)

    outstuff = np.zeros((N_x, int(float(T)/dt)))
    t = 0.0
    n = 0

    while t < T:
        for i in range(1, N_x-1):
            zeta_curr[i] = zeta_pp[i] - gamma*(psi_prev[i+1] - psi_prev[i-1])
        for i in range(1,N_x-1):
            rhs_diag[i-1] = -dx2*zeta_curr[i]
        psi_curr = tridiag(diag, rhs_diag, N_x -2, psi_curr)
    
        for i in range(1, N_x -1):
            psi_prev[i] = psi_curr[i]
            zeta_pp[i] = zeta_prev[i]
            zeta_prev[i] = zeta_curr[i]
        t += dt
        if (n % 20 == 0):
            outstuff[:, n] = psi_curr[:]
        n += 1

    return outstuff   

def tridiagonal_solidb(psi, rhs, dx, N):
    e = np.ones(N)
    d = np.ones(N)*-2
    
    psi[N-1] = 0.0
    psi[0] = 0.0
    
    for i in range(2, N):
        d[i] -=(e[i]*e[i-1]/d[i-1])
        rhs[i] -= ((e[i])*(rhs[i-1]))/d[i-1]
        
    for i in range(N-1, 1, -1):
        psi[i] = (rhs[i] - (e[i]*psi[i+1]))/d[i]

def center_diff(N_x, dx, T, dt):
    init_psi, init_zeta = initialize(N, dx)
    alpha = dt/(2*dx) #Ola
    gamma = dt/dx #ingvild
    dx2 = dx**2 #jostein
    
    zeta_0 = init_zeta
    zeta_curr = init_zeta
    for i in range(1, N_x -1):
        temp_zeta = init_zeta[i]
        temp_psi1 = psi[i+1]
        temp_psi2 = psi[i-1]
        zeta_curr[i] = temp_zeta +(temp_psi2 - temp_psi1)*alpha
        
    while t < T:
        pass
    #bruker center_diff_solidb og tridiagonal_solidb
    
    



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

if __name__ == "__main__":

    T = 150
    dt = 1

    dx = 1.0/40
    L = 1.0
    N = int(L/dx + 1)
 

    outstuff2 = leapfrog(N, dx, T, dt)

    outstuff= euler(N, dx, T, dt)

    
    
#    psiE_gauss = euler(init_psi_gauss, init_zeta_gauss, N, dx, T, dt)
#    psiLF_gauss = leapfrog(init_psi_gauss, init_zeta_gauss, N, dx, T, dt)
#    
    x = np.linspace(0, 1, N-1)
    
    plt.figure()
    plt.plot(x, outstuff[:-1,0], 'r-', label = 'Euler')
    plt.plot(x, outstuff2[:-1,0], 'b-.', label = 'Leapfrog')
    plt.legend()
    plt.title(r'Streamfunction $\psi(x, t)$ at $t = {}$ with $\Delta t = {:.1f}$'\
              .format(T, dt), fontsize = 15)
    plt.xlabel('x', fontsize = 12)
    plt.ylabel(r'$\psi(x,t)$', fontsize = 12)
    plt.grid()
    plt.savefig('figs/boundary_T{}_dt{}.pdf'.format(T, str(dt)), bbox_inches = 'tight')
    plt.show()
#    plt.figure()
#    plt.plot(x, psiE_gauss[1:N-3], 'r-')
#    plt.plot(x, psiLF_gauss[1:N-3], 'b-.')




