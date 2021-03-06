#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec  3 10:33:10 2018

@author: maritkollstuen
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation


def tridiag(psi, d, N):

    """
    Function to make tridiagonal matrix to solve our 1D PDE with boundary conditions.
    Returns psi which is our solution

    Input:
        psi     <float64> array with inital values of psi
        d       <float64> array, right handside of Ax = y, in our case zeta
        N       <int>  number of points to loop through

    Output:
        psi     <float64> array with solution of our PDE
     """

    b = np.ones(N-1)*(1.0)
    a = np.ones(N)*(-2.0)

    for i in range(1, N):
        a[i] -= b[i-1]**2/a[i-1]
        d[i] -= d[i-1]*b[i-1]/a[i-1]

    psi[N-1] = d[N-1]/a[N-1]

    for i in range(N-2, -1, -1):
        psi[i] = (d[i]-b[i-1]*psi[i+1])/a[i]

    return psi

def periodic_matrix(n_rows, n_cols):
    """
    Function to make periodic matrix to solve the 1D PDE in a periodic domain.
    The function returns a n_rows times n_cols matrix which is used in a matrix solver.

    Input:
        n_rows      <int> number of rows in matrix
        n_cols      <int> number of columns in matrix

    Output:
        A           <float64> matrix with size a n_rows times n_cols
    """
    A = np.zeros((n_rows, n_cols))
    for i in range(n_rows):
        for j in range(n_cols):
            if i == j:
                A[i, j] = -2.0
            elif abs(i-j) ==1:
                A[i,j] = 1.0

    return A


def initialize(N, dx, case, sigma = None):
    """
    Function to set the inital values of the sine and the gaussian wave.
    This will be used in a function below to calculate the full solution.

    Input:
        N          <int> number of points for loop
        dx         <float> stepsize in x-direction
        case       <str> select which case to initialize, sine or gaussian
        sigma      <float> width of the gaussian wave. Only used for case == gauss

    Output
        init_psi    <float64> array with initial values for psi
        init_zeta   <float64> array with inital values for zeta

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
            init_psi[i] = np.exp(-((x-0.5)/sigma)**2)
            init_zeta[i] = 4.0*((x-0.5)/sigma**2)**2 - (2/sigma**2)*np.exp(-((x-0.5)/sigma)**2)


    return init_psi, init_zeta




def euler_fwd(N_x, dx, T, dt, case, sigma = None):
    """
    Function to use in forward euler (forward difference) method.
    Uses functions initialize and tridiag in order to solve the PDE in question.
    Boundary conditions are set in order to have no flow into walls.


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
    alpha = dt/(2.0*dx)
    dx2 = dx**2
    bc_0 = 0.0
    bc_N = 0.0

    psi_prev = np.zeros(N_x)
    psi_curr = np.zeros(N_x)
    zeta_prev = np.zeros(N_x)
    zeta_curr = np.zeros(N_x)


    psi_prev  = psi_0
    zeta_prev = zeta_0

    psi_curr[0] = bc_0; psi_curr[-1] = bc_N
    zeta_curr[0] = zeta_prev[0]; zeta_curr[-1] = zeta_prev[-1]

    out_data = np.zeros((N_x+1, int(float(T)/dt)+2))
    t = 0.0
    out_data[0,0] = t
    out_data[1:, 0] = psi_0[:]
    n = 0
    n2 = 1
    while t < T:
        # forward Euler:
        for i in range(1, N_x-2):
            zeta_curr[i] = zeta_prev[i] - alpha*(psi_prev[i+1] - psi_prev[i-1])


        rhs = dx2*zeta_curr[1:-1]


        psi_curr[1:-1] = tridiag(psi_curr[1:-1], rhs, N_x-2)
        for i in range(1, N_x-2):
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
    in question. Boundary conditions are set in order to have to flow into the walls.
    Uses functions initialize and tridiag in order to solve the PDE in question.


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
    dx2 = dx**2
    bc_0 = 0.0
    bc_N = 0.0
    gamma = dt/dx

    psi_prev = np.zeros(N_x)
    psi_curr = np.zeros(N_x)
    zeta_prev = np.zeros(N_x)
    zeta_pp = np.zeros(N_x)
    zeta_curr = np.zeros(N_x)



    psi_prev  = psi_0
    zeta_pp = zeta_0
    #initial conditions are set with no flow into walls
    psi_curr[0] = bc_0; psi_curr[-1] = bc_N
    zeta_curr[0] = zeta_pp[0]; zeta_curr[-1] = zeta_pp[-1]


    #initial Euler:
    for i in range(1, N_x-1):
        zeta_prev[i] = zeta_0[i] - alpha*(psi_0[i+1] - psi_0[i-1])

    out_data = np.zeros((N_x+1, int(float(T)/dt)))
    t = 0.0
    out_data[0,0] = t
    out_data[1:, 0] = psi_0[:]
    n = 0
    n2 = 1

    while t < T:
        for i in range(1, N_x-1):
            zeta_curr[i] = zeta_pp[i] - gamma*(psi_prev[i+1] - psi_prev[i-1])


        rhs_diag = dx2*zeta_curr[1:-1]
        psi_curr[1:-1] = tridiag(psi_curr[1:-1], rhs_diag, N_x-2)


        psi_prev = psi_curr
        zeta_pp = zeta_prev
        zeta_prev = zeta_curr

        t += dt
        if (n % 20 == 0):
            out_data[0, n2] = t
            out_data[1:, n2] = psi_curr[:]

            n2 += 1

        n += 1

    return out_data

def animate_wave(x, t, psi):
    """
    Function to animate wave solution. Not included in report but used during the work
    to ensure our wave propagate west
    """
    plt.style.use("ggplot")
    fig = plt.figure(figsize = (8, 8))
    ax = plt.axes(xlim = (x[0], x[-1]), ylim = (psi.min(), psi.max()))
    ax.set_title("Streamfunction $\psi$ after {:.2f} time".format(t[0]),
        fontname = "serif", fontsize = 18)
    ax.set_xlabel("x [dim-less]", fontname = "serif", fontsize = 13)
    ax.set_ylabel("$\psi(x,t)$ [dim-less]", fontname = "serif", fontsize = 13)

    wave = ax.plot(x, psi[0, :], linewidth = 2)[0]


    def animate(frame):
        wave.set_data(x, psi[frame, :])
        ax.set_title("Streamfunction $\psi(x,t)$ after {:.2f} time".format(t[frame]),
            fontname = "serif", fontsize = 18)
        return []

    anim = animation.FuncAnimation(fig, animate, frames = psi.shape[0], interval = 10, blit = True)

    return anim



if __name__ == "__main__":
    # Calling function and make output

    T = 200
    dt = 0.001

    dx = 1.0/40
    L = 1.0
    N = int(L/dx + 1)


    psi_center_sine = leapfrog(N, dx, T, dt, 'sine')

    x = np.linspace(0,L,N)
    t = psi_center_sine[0,:400]

    # Make variables using leapfrog and sine/gaussian

    psi_center_sine = leapfrog(N, dx, T, dt, 'sine')
    psi_center_gauss_10 = leapfrog(N, dx, T, dt, 'gauss', sigma = 0.1)
    psi_center_gauss_25 = leapfrog(N, dx, T, dt, 'gauss', sigma = 0.25)
    psi_center_gauss_50 = leapfrog(N, dx, T, dt, 'gauss', sigma = 0.5)

    x = np.linspace(0,L,N)
    t = psi_center_sine[0,:400]

    # Plotting Hovmuller diagrams for sine and different gaussian

    plt.figure(1, figsize = (8, 10))
    plt.subplot(221)
    plt.style.use("ggplot")
    CS = plt.contourf(x, t, psi_center_sine[1:, :400].transpose(), 20, cmap = plt.cm.coolwarm)
    plt.colorbar(CS, orientation = "horizontal")
    plt.xlabel('x', fontsize = 13)
    plt.ylabel('time, t', fontsize = 13)
    plt.title(r'$\psi(x, t)$ sine wave', fontsize = 12)

    plt.subplot(222)
    plt.style.use("ggplot")
    CS = plt.contourf(x, t, psi_center_gauss_10[1:, :400].transpose(), 20, cmap = plt.cm.coolwarm)
    plt.colorbar(CS, orientation = "horizontal")
    plt.xlabel('x', fontsize = 13)
    plt.title(r'$\psi(x, t)$ gaussian wave $\sigma=0.1$', fontsize=12)

    plt.subplot(223)
    plt.style.use("ggplot")
    CS = plt.contourf(x, t, psi_center_gauss_25[1:, :400].transpose(), 20, cmap = plt.cm.coolwarm)
    plt.colorbar(CS, orientation = "horizontal")
    plt.xlabel('x', fontsize = 13)
    plt.ylabel('time, t', fontsize = 13)
    plt.title(r'$\psi(x, t)$ gaussian wave $\sigma=0.25$', fontsize=12)

    plt.subplot(224)
    plt.style.use("ggplot")
    CS = plt.contourf(x, t, psi_center_gauss_50[1:, :400].transpose(), 20, cmap = plt.cm.coolwarm)
    plt.colorbar(CS, orientation = "horizontal")
    plt.xlabel('x', fontsize = 13)
    plt.title(r'$\psi(x, t)$ gaussian wave $\sigma=0.5$', fontsize=12)

    plt.savefig('figs/boundary_1D_hovmuller.pdf', bbox_inches = 'tight')



    psi_center = leapfrog(N, dx, T, dt, 'sine')
    psi_euler= euler_fwd(N, dx, T, dt, 'sine')

    dt2 = 0.01

    psi_center2 = leapfrog(N, dx, T, dt2, 'sine')

    psi_euler2 = euler_fwd(N, dx, T, dt2, 'sine')

    dt3 = 0.2

    psi_center3 = leapfrog(N, dx, T, dt3, 'sine')

    psi_euler3 = euler_fwd(N, dx, T, dt3, 'sine')

    dt4 = 1.0

    psi_center4 = leapfrog(N, dx, T, dt4, 'sine')

    psi_euler4 = euler_fwd(N, dx, T, dt4, 'sine')

    x = np.linspace(0, 1, N-2)

    # Plotting waves with euler and leapfrog to see stability


    plt.figure(2, figsize = (10, 8))
    plt.subplot(221)
    plt.plot(x, psi_euler[2:-1,1], 'r-', label = 'Euler')
    plt.plot(x, psi_center[2:-1:,1], 'b-.', label = 'Leapfrog')
    plt.legend()
    plt.title(r'$\Delta t = {:.3f}$'\
              .format(dt), fontsize = 15)
    plt.ylabel(r'$\psi(x,t)$', fontsize = 12)
    plt.grid()
    plt.subplot(222)
    plt.plot(x, psi_euler2[2:-1,1], 'r-', label = 'Euler')
    plt.plot(x, psi_center2[2:-1,1], 'b-.', label = 'Leapfrog')
    plt.legend()
    plt.title(r'$\Delta t = {:.3f}$'\
              .format(dt2), fontsize = 15)
    plt.grid()

    plt.subplot(223)
    plt.plot(x, psi_euler3[2:-1,1], 'r-', label = 'Euler')
    plt.plot(x, psi_center3[2:-1,1], 'b-.', label = 'Leapfrog')
    plt.legend()
    plt.title(r'$\Delta t = {:.2f}$'\
              .format(dt3), fontsize = 15)
    plt.xlabel('x', fontsize = 12)
    plt.ylabel(r'$\psi(x,t)$', fontsize = 12)
    plt.grid()

    plt.subplot(224)
    plt.plot(x, psi_euler4[2:-1,1], 'r-', label = 'Euler')
    plt.plot(x, psi_center4[2:-1,1], 'b-.', label = 'Leapfrog')
    plt.legend()
    plt.title(r'$\Delta t = {:.1f}$'\
              .format(dt4), fontsize = 15)
    plt.xlabel('x', fontsize = 12)
    plt.grid()
    plt.savefig('figs/subplots_boundary1D.pdf')

#
