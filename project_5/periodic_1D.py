#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec  3 15:22:26 2018

@author: maritkollstuen
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.sparse import linalg
import sys, os


def periodic_matrix(n_rows, n_cols):
    A = np.zeros((n_rows, n_cols))
    for i in range(n_rows):
        for j in range(n_cols):
            if i == j:
                A[i, j] = -2.0
            elif abs(i-j) ==1:
                A[i,j] = 1.0
    A[0, n_cols-1] = 1.0
    A[n_rows -1, 0] = 1.0
    
    return A


def initialize(N, dx, case):
    init_psi = np.zeros(N)
    init_zeta = np.zeros(N)
    if case == 'sine':
        for i in range(0, N):
            x = i*dx
            init_psi[i] = np.sin(4.0*np.pi*x)
            init_zeta[i] = -16.0*np.pi**2*np.sin(4.0*np.pi*x)
    if case == 'gauss':
        for i in range(0, N):
            x = i*dx
            sigma = 0.1
            init_psi[i] = np.exp(-((x-0.5)/sigma)**2)
            init_zeta[i] = 4.0*((x-0.5)/sigma**2)**2 - (2/sigma**2)*np.exp(-((x-0.5)/sigma)**2)

    return init_psi, init_zeta

def euler_fwd(N_x, dx, T, dt, case):
    psi_0, zeta_0 = initialize(N_x, dx, case)
    alpha = dt/(dx*2)
    dx2 = dx**2
    
    psi_prev = np.zeros(N_x)
    psi_curr = np.zeros(N_x)
    zeta_prev = np.zeros(N_x)
    zeta_curr = np.zeros(N_x)

    
    rhs_poisson = np.zeros(N_x-2)
    
    
    
    psi_prev  = psi_0
    zeta_prev = zeta_0
    
    out_data = np.zeros((N_x+1, int(float(T)/dt)+2))
    t = 0.0
    out_data[0,0] = t
    out_data[1:, 0] = psi_0[:]
    n = 0
    n2 = 1
    while t < T:
#        #forward Euler:
        for i in range(1, N_x-1):
            zeta_curr[i] = zeta_prev[i] - alpha*(psi_prev[i+1] - psi_prev[i-1])
            
        zeta_curr[0] = zeta_prev[0] - alpha*(psi_prev[1] - psi_prev[-2])
        zeta_curr[-1] = zeta_curr[0]
    
        rhs_poisson = dx2*zeta_curr[1:-1]
            
        A = periodic_matrix(int(N_x-2),int(N_x-2))
        
#        psi_curr[1:-1] = np.linalg.solve(A, rhs_poisson)
        sys.stdout = open(os.devnull, "w")
        psi_curr[1:-1],istop, itn, normr, normar, norma, conda, normx, bla, bkabla = linalg.lsqr(A, rhs_poisson)
        sys.stdout = sys.__stdout__
            
        
        psi_curr[-1] = psi_curr[0]
        
        for i in range(0, N_x-1):
            psi_prev[i] = psi_curr[i]
            zeta_prev[i] = zeta_curr[i]
        
        t += dt
        if (n % 20 == 0):
            out_data[0, n2] = t
            out_data[1:, n2] = psi_curr[:]
            n2 += 1
 
        n += 1

    return out_data   
        
def center(N_x, dx, T, dt, case):
    psi_0, zeta_0 = initialize(N_x, dx, case)
    alpha = dt/(2*dx)
    gamma =  dt/dx
    dx2 = dx**2
    
    psi_prev = np.zeros(N_x)
    psi_curr = np.zeros(N_x)
    zeta_prev = np.zeros(N_x)

    zeta_curr = np.zeros(N_x)

    # To solve Ax = b
    rhs_poisson = np.zeros(N_x-2)
    
    
    
    # Initial conditions
    psi_prev  = psi_0
    zeta_prev = zeta_0

    # initial Euler:
    for i in range(1, N_x-1):
        zeta_prev[i] = zeta_0[i] - alpha*(psi_0[i+1] - psi_0[i-1])

    zeta_prev[0] = zeta_0[0] - alpha*(psi_0[1] - psi_0[-2])
    zeta_prev[-1] = zeta_0[0]    

    
    # store data
    out_data = np.zeros((N_x+1, int(float(T)/dt)+1))
    t = 0.0
    out_data[0,0] = t
    out_data[1:, 0] = psi_0[:]
    n = 0
    n2 = 1
    
    # Loop over time using centered difference
    while t < T:
        # forward Euler:
        for i in range(1, N_x-1):
            zeta_curr[i] = zeta_prev[i] - gamma*(psi_prev[i+1] - psi_prev[i-1])
        
        zeta_curr[0] = zeta_prev[0] - gamma*(psi_prev[1] - psi_prev[-2])
        zeta_curr[-1] = zeta_curr[0]
        

        rhs_poisson = dx2*zeta_curr[1:-1]
        
        A = periodic_matrix(int(N_x-2),int(N_x-2))
        sys.stdout = open(os.devnull, "w")
        psi_curr[1:-1],istop, itn, normr, normar, norma, conda, normx, bla, bkabla = linalg.lsqr(A, rhs_poisson)
        sys.stdout = sys.__stdout__
            
        psi_curr[-1] = psi_curr[0]
        
        
        for i in range(0, N_x-1):
            psi_prev[i] = psi_curr[i]
            zeta_prev[i] = zeta_curr[i]
        
        t += dt
        #saving every 20th point
        if (n % 20 == 0):
            out_data[0, n2] = t
            out_data[1:, n2] = psi_curr[:]
            
            n2 += 1

        n += 1

    return out_data   


        
if __name__ == "__main__":

    T = 200
    dt = 0.2
    
    dx = 1.0/40
    L = 1.0
    N = int(L/dx + 1)
    


    psi_center_sine = center(N, dx, T, dt, case = 'sine')
#    psi_center_gauss = euler_fwd(N, dx, T, dt, case = 'gauss')
    
##    psi_center_sine = psi_center_sine[:,0::20]
##    psi_center_gauss = psi_center_gauss[:,0::20]
##    
    x = np.linspace(0,L,N)
    t = np.linspace(0,T,N)



#    plt.figure(1)
#    plt.contourf(t, x, psi_center_sine[1:, :40])
#    plt.colorbar()
#    plt.show()
#    
#    plt.figure(2)
#    plt.contourf(t, x,psi_center_gauss[1:, :40])
#    plt.colorbar()
#    plt.show()  
#    
#    plt.figure(3)
#    plt.style.use("ggplot")
#    fig = plt.figure(figsize = (9,7))
#    CS = plt.contourf(x, t, psi_center_gauss[1:, :41].transpose(), 20, cmap = plt.cm.RdBu_r)
#    plt.colorbar(CS, orientation = "vertical")
#    plt.xlabel('x', fontsize = 13)
#    plt.ylabel('time, t', fontsize = 13)
#    plt.title(r'Hovmüller diagram of $\psi(x, t)$')
#    
    plt.figure(4)
    plt.style.use("ggplot")
    fig = plt.figure(figsize = (9,7))
    CS = plt.contourf(x, t, psi_center_sine[1:, :41].transpose(), 20, cmap = plt.cm.RdBu_r)
    plt.colorbar(CS, orientation = "vertical")
    plt.xlabel('x', fontsize = 13)
    plt.ylabel('time, t', fontsize = 13)
    plt.title(r'Hovmüller diagram of $\psi(x, t)$')
    
#    psi_center = center(N, dx, T, dt, case = 'sine')
#
#    psi_euler= euler_fwd(N, dx, T, dt, case = 'sine')
#
#    dt2 = 0.01
#    
#    psi_center2 = center(N, dx, T, dt2, case = 'sine')
#    
#    psi_euler2 = euler_fwd(N, dx, T, dt2, case = 'sine')
#    
#    dt3 = 0.2
#    
#    psi_center3 = center(N, dx, T, dt3, case = 'sine')
#    
#    psi_euler3 = euler_fwd(N, dx, T, dt3, case = 'sine')
#    
#    dt4 = 1.0
#    
#    psi_center4 = center(N, dx, T, dt4, case = 'sine')
#    
#    psi_euler4 = euler_fwd(N, dx, T, dt4, case = 'sine')
#    
#    x = np.linspace(0, 1, N)
#
#
#    
# 
#
#    plt.figure(2, figsize = (10, 8))
#    plt.subplot(221)
#    plt.plot(x, psi_euler[1:,1], 'r-', label = 'Euler')
#    plt.plot(x, psi_center[1:,1], 'b-.', label = 'Centered')
#    plt.legend()
#    plt.title(r'$\Delta t = {:.3f}$'\
#              .format(dt), fontsize = 15)
##    plt.xlabel('x', fontsize = 12)
#    plt.ylabel(r'$\psi(x,t)$', fontsize = 12)
#    plt.grid()
#    plt.subplot(222)
#    plt.plot(x, psi_euler2[1:,1], 'r-', label = 'Euler')
#    plt.plot(x, psi_center2[1:,1], 'b-.', label = 'Centered')
#    plt.legend()
#    plt.title(r'$\Delta t = {:.3f}$'\
#              .format(dt2), fontsize = 15)
##    plt.xlabel('x', fontsize = 12)
##    plt.ylabel(r'$\psi(x,t)$', fontsize = 12)
#    plt.grid()
#    plt.subplot(223)
#    plt.plot(x, psi_euler3[1:,1], 'r-', label = 'Euler')
#    plt.plot(x, psi_center3[1:,1], 'b-.', label = 'Centered')
#    plt.legend()
#    plt.title(r'$\Delta t = {:.2f}$'\
#              .format(dt3), fontsize = 15)
#    plt.xlabel('x', fontsize = 12)
#    plt.ylabel(r'$\psi(x,t)$', fontsize = 12)
#    plt.grid()
#    plt.subplot(224)
#    plt.plot(x, psi_euler4[1:,1], 'r-', label = 'Euler')
#    plt.plot(x, psi_center4[1:,1], 'b-.', label = 'Centered')
#    plt.legend()
#    plt.title(r'$\Delta t = {:.1f}$'\
#              .format(dt4), fontsize = 15)
#    plt.xlabel('x', fontsize = 12)
##    plt.ylabel(r'$\psi(x,t)$', fontsize = 12)
#    plt.grid()
#    plt.savefig('figs/subplots_periodic1D.pdf')





#   
#    outstuff = euler_fwd(N, dx, T, dt, case = 'sine')
#    outstuff2 = center(N, dx, T, dt, case = 'sine')
#    plt.figure()
#    plt.plot(outstuff, 'r-', label = 'Euler')
#    plt.plot(outstuff2[1:,1], 'b-.', label = 'Centered')
#    plt.legend()
#    plt.title(r'Streamfunction $\psi(x, t)$ at $t = {}$ with $\Delta t = {:.3f}$'\
#          .format(T, dt), fontsize = 15)
#    plt.xlabel('x', fontsize = 12)
#    plt.ylabel(r'$\psi(x,t)$', fontsize = 12)
#    plt.grid()

        
        
        
        
        
        
        

    
    