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
        """
        Function to make periodic matrix to solve the 1D PDE in a periodic domain.
        The function returns a n_rows times n_cols matrix which is used in a matrix solver.

        Input:
            n_rows      <int>       number of rows in matrix
            n_cols      <int>       number of columns in matrix

        Output:
            A           <float64>   matrix with size a n_rows times n_cols
        """
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


def initialize(N, dx, case, sigma = None):
    """
    Function to set the inital values of the sine and the gaussian wave.
    This will be used in a function below to calculate the full solution.

    Input:
        N          <int>        number of points for loop
        dx         <float>      stepsize in x-direction
        case       <str>        select which case to initialize, sine or gaussian
        sigma      <float>      width of the gaussian wave. Only used for case == gauss

    Output
        init_psi    <float64>   array with initial values for psi
        init_zeta   <float64>   array with inital values for zeta

    """
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
            sigma = sigma
            init_psi[i] = np.exp(-((x-0.5)/sigma)**2)
            init_zeta[i] = 4.0*((x-0.5)/sigma**2)**2 - (2/sigma**2)*np.exp(-((x-0.5)/sigma)**2)

    return init_psi, init_zeta

def euler_fwd(N_x, dx, T, dt, case, sigma = None):
    """
    Function to use in forward euler (forward difference) method.
    Uses functions initialize and tridiag.
    linalg.lsqr is used in order to solve our sparse matrix and find psi

    Input:
        N_x     <int>       number of points in x-direction
        dx      <float>     stepsize in x-direction
        T       <int>       time limit
        dt      <float>     time stepsize
        case    <str>       select which case to initialize, sine or gaussian
        sigma   <float>     width of the gaussian wave. Only used for case == gauss

    Output:
        out_data <float64>  array with time in first row and psi (solution)
    """

    psi_0, zeta_0 = initialize(N_x, dx, case, sigma)
    alpha = dt/(dx*2)
    dx2 = dx**2

    psi_prev = np.zeros(N_x)
    psi_curr = np.zeros(N_x)
    zeta_prev = np.zeros(N_x)
    zeta_curr = np.zeros(N_x)


    rhs = np.zeros(N_x-2)



    psi_prev  = psi_0
    zeta_prev = zeta_0

    out_data = np.zeros((N_x+1, int(float(T)/dt)+2))
    t = 0.0
    out_data[0,0] = t
    out_data[1:, 0] = psi_0[:]
    n = 0
    n2 = 1
    while t < T:
        for i in range(1, N_x-1):
            zeta_curr[i] = zeta_prev[i] - alpha*(psi_prev[i+1] - psi_prev[i-1])

        zeta_curr[0] = zeta_prev[0] - alpha*(psi_prev[1] - psi_prev[-2])
        zeta_curr[-1] = zeta_curr[0]

        rhs = dx2*zeta_curr[1:-1]

        A = periodic_matrix(int(N_x-2),int(N_x-2))


        sys.stdout = open(os.devnull, "w")
        psi_curr[1:-1],istop, itn, normr, normar, norma, conda, normx, bla, bkabla = linalg.lsqr(A, rhs)
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

def leapfrog(N_x, dx, T, dt, case, sigma = None):
        """
        Function to use for time stepping using leapfrog (centered difference) to solve the PDE
        in question.
        Uses functions initialize and tridiag.
        linalg.lsqr is used in order to solve our sparse matrix and find psi



        Input:
            N_x     <int>       number of points in x-direction
            dx      <float>     stepsize in x-direction
            T       <int>       time limit
            dt      <float>     time stepsize
            case    <str>       select which case to initialize, sine or gaussian
            sigma   <float>     width of the gaussian wave. Only used for case == gauss

        Output:
            out_data <float64>  array with time in first row and psi (solution)

        """
    psi_0, zeta_0 = initialize(N_x, dx, case, sigma)
    alpha = dt/(2*dx)
    gamma =  dt/dx
    dx2 = dx**2

    psi_prev = np.zeros(N_x)
    psi_curr = np.zeros(N_x)
    zeta_prev = np.zeros(N_x)

    zeta_curr = np.zeros(N_x)

    # To solve Ax = b
    rhs = np.zeros(N_x-2)


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


        rhs = dx2*zeta_curr[1:-1]

        A = periodic_matrix(int(N_x-2),int(N_x-2))
        sys.stdout = open(os.devnull, "w")
        psi_curr[1:-1],istop, itn, normr, normar, norma, conda, normx, bla, bkabla = linalg.lsqr(A, rhs)
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
    """
    Calling function with variables as set below.
    Plotting Hovmuller diagrams in 1D and
    """

    T = 200
    dt = 0.001

    dx = 1.0/40
    L = 1.0
    N = int(L/dx + 1)


    psi_euler_sine = euler_fwd(N, dx, T, dt,case = 'sine')

    psi_center_sine = leapfrog(N, dx, T, dt, case = 'sine')
    psi_center_gauss_10 = leapfrog(N, dx, T, dt, case = 'gauss', sigma = 0.1)
    psi_center_gauss_25 = leapfrog(N, dx, T, dt, case = 'gauss', sigma = 0.25)
    psi_center_gauss_09 = leapfrog(N, dx, T, dt, case = 'gauss', sigma = 0.09)

    psi_center_sine = psi_center_sine[:,0::20]
    psi_center_gauss = psi_center_gauss[:,0::20]
   
    x = np.linspace(0,L,N)
    t = psi_center_sine[0,:51]

    plt.figure(1, figsize = (8, 10))
    plt.subplot(221)
    plt.style.use("ggplot")
    CS = plt.contourf(x, t, psi_center_sine_10[1:, :51].transpose(), 20, cmap = plt.cm.coolwarm)
    plt.colorbar(CS, orientation = "horizontal")
    plt.xlabel('x', fontsize = 13)
    plt.ylabel('time, t', fontsize = 13)
    plt.title(r' $\psi(x, t)$ sine wave', fontsize=12)
#
    plt.subplot(222)
    plt.style.use("ggplot")
    CS = plt.contourf(x, t, psi_center_gauss_09[1:, :51].transpose(), 20, cmap = plt.cm.coolwarm)
    plt.colorbar(CS, orientation = "horizontal")
    plt.xlabel('x', fontsize = 13)
    plt.title(r'$\psi(x, t)$ gaussian wave $\sigma = 0.09$', fontsize=12)
##
    plt.subplot(223)
    plt.style.use("ggplot")
    CS = plt.contourf(x, t, psi_center_gauss[1:, :51].transpose(), 20, cmap = plt.cm.coolwarm)
    plt.colorbar(CS, orientation = "horizontal")
    plt.xlabel('x', fontsize = 13)
    plt.ylabel('time, t', fontsize = 13)
    plt.title(r'$\psi(x, t)$ gaussian wave $\sigma = 0.1$', fontsize=12)

    plt.subplot(224)
    plt.style.use("ggplot")
    CS = plt.contourf(x, t, psi_center_gauss_25[1:, :51].transpose(), 20, cmap = plt.cm.coolwarm)
    plt.colorbar(CS, orientation = "horizontal")
    plt.xlabel('x', fontsize = 13)
    plt.title(r'$\psi(x, t)$ gaussian wave $\sigma = 0.25$', fontsize=12)

    plt.savefig('figs/periodic_1D_hovmuller.pdf', bbox_inches = 'tight')

    psi_center = leapfrog(N, dx, T, dt, case = 'sine')

    psi_euler= euler_fwd(N, dx, T, dt, case = 'sine')

    dt2 = 0.01

    psi_center2 = leapfrog(N, dx, T, dt2, case = 'sine')

    psi_euler2 = euler_fwd(N, dx, T, dt2, case = 'sine')

    dt3 = 0.2

    psi_center3 = leapfrog(N, dx, T, dt3, case = 'sine')

    psi_euler3 = euler_fwd(N, dx, T, dt3, case = 'sine')

    dt4 = 1.0

    psi_center4 = leapfrog(N, dx, T, dt4, case = 'sine')

    psi_euler4 = euler_fwd(N, dx, T, dt4, case = 'sine')

    x = np.linspace(0, 1, N)





    plt.figure(2, figsize = (10, 8))
    plt.subplot(221)
    plt.plot(x, psi_euler[1:,1], 'r-', label = 'Euler')
    plt.plot(x, psi_center[1:,1], 'b-.', label = 'Leapfrog')
    plt.legend()
    plt.title(r'$\Delta t = {:.3f}$'\
              .format(dt), fontsize = 15)
    plt.ylabel(r'$\psi(x,t)$', fontsize = 12)
    plt.grid()
    plt.subplot(222)
    plt.plot(x, psi_euler2[1:,1], 'r-', label = 'Euler')
    plt.plot(x, psi_center2[1:,1], 'b-.', label = 'Leapfrog')
    plt.legend()
    plt.title(r'$\Delta t = {:.3f}$'\
              .format(dt2), fontsize = 15)

    plt.grid()
    plt.subplot(223)
    plt.plot(x, psi_euler3[1:,1], 'r-', label = 'Euler')
    plt.plot(x, psi_center3[1:,1], 'b-.', label = 'Leapfrog')
    plt.legend()
    plt.title(r'$\Delta t = {:.2f}$'\
              .format(dt3), fontsize = 15)
    plt.xlabel('x', fontsize = 12)
    plt.ylabel(r'$\psi(x,t)$', fontsize = 12)
    plt.grid()
    plt.subplot(224)
    plt.plot(x, psi_euler4[1:,1], 'r-', label = 'Euler')
    plt.plot(x, psi_center4[1:,1], 'b-.', label = 'Leapfrog')
    plt.legend()
    plt.title(r'$\Delta t = {:.1f}$'\
              .format(dt4), fontsize = 15)
    plt.xlabel('x', fontsize = 12)
    plt.grid()

    plt.savefig('figs/subplots_periodic1D.pdf')
