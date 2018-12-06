#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec  3 10:33:10 2018

@author: maritkollstuen
"""

import numpy as np
import matplotlib.pyplot as plt


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

def initialize(N, dx):
    init_psi = np.zeros(N)
    init_zeta = np.zeros(N)

    for i in range(0, N-1):
        x = i*dx
        init_psi[i] = np.sin(4.0*np.pi*x) 
        init_zeta[i] = -16.0*np.pi**2*np.sin(4.0*np.pi*x)

    return init_psi, init_zeta



def euler_fwd(N_x, dx, T, dt):

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

    out_data = np.zeros((N_x, int(float(T)/dt)))
    t = 0.0
    n = 0

    while t < T:
        # forward Euler:
        for i in range(1, N_x-2):
            zeta_curr[i] = zeta_prev[i] - alpha*(psi_prev[i+1] - psi_prev[i-1])
        for i in range(1, N_x-2):
            rhs_diag[i-1] = -dx2*zeta_curr[i]
        psi_curr = tridiag(diag, rhs_diag, N_x -2, psi_curr)
        for i in range(1, N_x-2):
            psi_prev[i] = psi_curr[i]
            zeta_prev[i] = zeta_curr[i]
        
        t += dt
        if (n % 20 == 0):
            out_data[:, n] = psi_curr[:]
        n += 1
 
    return out_data
#    return psi_curr

def center(N_x, dx, T, dt):
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
    for i in range(1, N_x-2):
        zeta_prev[i] = zeta_0[i] - alpha*(psi_0[i+1] - psi_0[i-1])
#    for i in range(1, N_x-2):
#        rhs_diag[i-1] = -dx2*zeta_prev[i]
#
#   
#    psi_prev = tridiag(diag, rhs_diag, N_x-2, psi_prev)
#    print(psi_prev[1:10])
    out_data = np.zeros((N_x, int(float(T)/dt)))
    t = 0.0
    n = 0

    while t < T:
        for i in range(1, N_x-2):
            zeta_curr[i] = zeta_pp[i] - gamma*(psi_prev[i+1] - psi_prev[i-1])
        for i in range(1,N_x-2):
            rhs_diag[i-1] = -dx2*zeta_curr[i]
        psi_curr = tridiag(diag, rhs_diag, N_x -2, psi_curr)
    
        for i in range(1, N_x -2):
            psi_prev[i] = psi_curr[i]
            zeta_pp[i] = zeta_prev[i]
            zeta_prev[i] = zeta_curr[i]
        t += dt
        if (n % 20 == 0):
            out_data[:, n] = psi_curr[:]
        n += 1

    return out_data   



if __name__ == "__main__":

    T = 150
    dt = 0.2

    dx = 1.0/40
    L = 1.0
    N = int(L/dx + 1)
 

    psi_center = center(N, dx, T, dt)

    psi_euler= euler_fwd(N, dx, T, dt)

    dt2 = 0.01
    
    psi_center2 = center(N, dx, T, dt2)
    
    psi_euler2 = euler_fwd(N, dx, T, dt2)
    
    dt3 = 0.2
    
    psi_center3 = center(N, dx, T, dt3)
    
    psi_euler3 = euler_fwd(N, dx, T, dt3)
    
    dt4 = 1.0
    
    psi_center4 = center(N, dx, T, dt4)
    
    psi_euler4 = euler_fwd(N, dx, T, dt4)
    
    x = np.linspace(0, 1, N-1)
#
#    plt.figure(1)
#    plt.plot(x, psi_euler[:-1,0], 'r-', label = 'Euler')
#    plt.plot(x, psi_center[:-1,0], 'b-.', label = 'Centered')
#    plt.legend()
#    plt.title(r'Streamfunction $\psi(x, t)$ at $t = {}$ with $\Delta t = {:.3f}$'\
#              .format(T, dt), fontsize = 15)
#    plt.xlabel('x', fontsize = 12)
#    plt.ylabel(r'$\psi(x,t)$', fontsize = 12)
#    plt.grid()
##    plt.savefig('figs/boundary_T{}_dt{}.pdf'.format(T, str(dt)), bbox_inches = 'tight')
#    plt.show()

#
    plt.figure(2, figsize = (10, 8))
    plt.subplot(221)
    plt.plot(x, psi_euler[:-1,0], 'r-', label = 'Euler')
    plt.plot(x, psi_center[:-1,0], 'b-.', label = 'Centered')
    plt.legend()
    plt.title(r'$\Delta t = {:.3f}$'\
              .format(dt), fontsize = 15)
#    plt.xlabel('x', fontsize = 12)
    plt.ylabel(r'$\psi(x,t)$', fontsize = 12)
    plt.grid()
    plt.subplot(222)
    plt.plot(x, psi_euler2[:-1,0], 'r-', label = 'Euler')
    plt.plot(x, psi_center2[:-1,0], 'b-.', label = 'Centered')
    plt.legend()
    plt.title(r'$\Delta t = {:.3f}$'\
              .format(dt2), fontsize = 15)
#    plt.xlabel('x', fontsize = 12)
#    plt.ylabel(r'$\psi(x,t)$', fontsize = 12)
    plt.grid()
    plt.subplot(223)
    plt.plot(x, psi_euler3[:-1,0], 'r-', label = 'Euler')
    plt.plot(x, psi_center3[:-1,0], 'b-.', label = 'Centered')
    plt.legend()
    plt.title(r'$\Delta t = {:.2f}$'\
              .format(dt3), fontsize = 15)
    plt.xlabel('x', fontsize = 12)
    plt.ylabel(r'$\psi(x,t)$', fontsize = 12)
    plt.grid()
    plt.subplot(224)
    plt.plot(x, psi_euler4[:-1,0], 'r-', label = 'Euler')
    plt.plot(x, psi_center4[:-1,0], 'b-.', label = 'Centered')
    plt.legend()
    plt.title(r'$\Delta t = {:.1f}$'\
              .format(dt4), fontsize = 15)
    plt.xlabel('x', fontsize = 12)
#    plt.ylabel(r'$\psi(x,t)$', fontsize = 12)
    plt.grid()
#    plt.savefig('figs/subplots_boundary1D.pdf')











