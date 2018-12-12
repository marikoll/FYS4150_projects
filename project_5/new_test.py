#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 11 09:05:04 2018

@author: maritkollstuen
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.sparse import linalg
import sys, os


def initialize(N, dx, dt):
    psi_0 = np.zeros(N)
    zeta_0 = np.zeros(N)
    for i in range(0, N):
        x = i*dx
        t = i*dt
        psi_0[i] = np.sin(np.pi*x)*np.cos(np.pi*x + t/np.pi)
        zeta_0[i] = -2*np.pi**2*(np.sin(np.pi*x)*np.cos(np.pi*x + t/np.pi) + np.cos(np.pi*x)*np.sin(np.pi*x + t/np.pi))
    return psi_0, zeta_0

def tridiag(psi, d, N): 
    
    b = np.ones(N-1)*(-1.0)
    a = np.ones(N)*2.0
    
    for i in range(1, N):
        a[i] -= b[i-1]**2/a[i-1]
        d[i] -= d[i-1]*b[i-1]/a[i-1]
    
    psi[N-1] = d[N-1]/a[N-1]

    for i in range(N-2, 0, -1):
        psi[i] = (d[i]-b[i-1]*psi[i+1])/a[i]
    psi[-1] = 0.0
    psi[0] = 0.0   
    return psi

def periodic_matrix(n_rows, n_cols):
    A = np.zeros((n_rows, n_cols))
    for i in range(n_rows):
        for j in range(n_cols):
            if i == j:
                A[i, j] = 2.0
            elif abs(i-j) ==1:
                A[i,j] = -1.0
    A[0, n_cols-1] = -1.0
    A[n_rows -1, 0] = -1.0
    
    return A

def forward_euler_basin(N, dx, T, dt):
    alpha = dt/(2*dx)
    dx2 = dx**2
    psi_0, zeta_0 = initialize(N, dx, dt)
    
    psi = psi_0
    zeta = zeta_0
    
    psi[0] = psi[-1] = 0.0
    
    time = np.linspace(dt, T, T)
    
    out = np.zeros((N+1, len(time)+1))
    out[0, 1:] = time
    out[1:, 0] = psi_0
    
    for k, t in enumerate(time):
        for i in range(1, N-1):
            zeta[i+1] = zeta[i] -alpha*(psi[i-1] - psi[i+1])
            psi = tridiag(psi, -zeta*dx2, N)
            out[1:, k+1] = psi
    return psi, zeta, out  

def forward_euler_periodic(N, dx, T, dt):
    alpha = dt/(2*dx)
    dx2 = dx**2
    psi_0, zeta_0 = initialize(N, dx, dt)
    
    psi = psi_0
    zeta = zeta_0
    
    psi[0] = psi[-1] = 0.0
    
    time = np.linspace(dt, T, T)
    
    out = np.zeros((N+1, len(time)+1))
    out[0, 1:] = time
    out[1:, 0] = psi_0
    
    for k, t in enumerate(time):
        for i in range(1, N-1):
            zeta[i+1] = zeta[i] -alpha*(psi[i-1] - psi[i+1])
            A = periodic_matrix(int(N-2),int(N-2))
            #print(len(psi[1:-1]), len(-zeta[1:-1]*dx2))
            sys.stdout = open(os.devnull, "w")
            psi[1:-1],istop, itn, normr, normar, norma, conda, normx, bla, bkabla = linalg.lsqr(A, -zeta[1:-1]*dx2)
            sys.stdout = sys.__stdout__
            
            out[1:, k+1] = psi
        zeta[0] = zeta_0[0] + alpha*(psi[N-2] - psi[1])
        zeta[-1] = zeta[0]
    return psi, zeta, out  


def centered_diff_basin(N, dx, T, dt):
    alpha = dt/(2*dx)
    gamma = dt/dx
    dx2 = dx**2
    psi_0, zeta_0 = initialize(N, dx, dt)
    
    psi = psi_0
    zeta = zeta_0
    
    psi[0] = psi[-1] = 0.0
    
    for i in range(1, N-1):
        zeta[i] = zeta[i] - alpha*(psi[i+1] - psi[i-1])
    
#    psi[1:-1] = tridiag(psi[1:-1], -zeta[1:-1]*dx2, N-2)
#    
    time = np.linspace(dt, T, T)

    out = np.zeros((N+1, len(time)+1))
    out[0, 1:] = time
    out[1:, 0] = psi_0
    
    
    for k, t in enumerate(time):
        for i in range(1, N-1):
            zeta[i+1] = zeta[i] -gamma*(psi[i-1] - psi[i+1])
            psi[1:-1] = tridiag(psi[1:-1], -zeta[1:-1]*dx2, N-2)
            out[1:, k+1] = psi
    return psi, zeta, out  



def centered_diff_periodic(N, dx, T, dt):
    alpha = dt/(2*dx)
    gamma = dt/dx
    dx2 = dx**2
    psi_0, zeta_0 = initialize(N, dx, dt)
    
    psi = psi_0
    zeta = zeta_0
    
    psi[0] = psi[-1] = 0.0
    
    for i in range(1, N-1):
        zeta[i] = zeta[i] - alpha*(psi[i+1] - psi[i-1])
    
#    A = periodic_matrix(int(N-2),int(N-2))
#    psi[1:-1] = linalg.solve(A, -zeta[1:-1]*dx2)
#    
    time = np.linspace(dt, T, T)

    out = np.zeros((N+1, len(time)+1))
    out[0, 1:] = time
    out[1:, 0] = psi_0
    
    
    for k, t in enumerate(time):
        for i in range(1, N-1):
            zeta[i+1] = zeta[i] -gamma*(psi[i-1] - psi[i+1])
            A = periodic_matrix(int(N-2),int(N-2))
            sys.stdout = open(os.devnull, "w")
            psi[1:-1],istop, itn, normr, normar, norma, conda, normx, bla, bkabla  = linalg.lsqr(A, -zeta[1:-1]*dx2)
            sys.stdout = sys.__stdout__
            out[1:, k+1] = psi
    return psi, zeta, out 


if __name__ == "__main__":
    T = 150
    dx = 1.0/40
    L = 1.0
    N = int(1.0/dx + 1)
    dt = 0.001
    
    psi1, zeta1, out_basin = forward_euler_basin(N, dx, T, dt)
    psi2, zeta2, out_periodic = forward_euler_periodic(N, dx, T, dt)    
    
    
    psi3, zeta3, out_c_periodic = centered_diff_periodic(N, dx, T, dt)
    
    psi4, zeta4, out_c_basin = centered_diff_basin(N, dx, T, dt)
    
    x = np.linspace(0,L,N)
    t = out_basin[0,:]
    
    
    plt.figure(1)
    plt.style.use("ggplot")
    fig = plt.figure(figsize = (9,7))
    CS = plt.contourf(x, t, out_basin[1:,:].transpose(), 20, cmap = plt.cm.RdBu_r)
    plt.colorbar(CS, orientation = "vertical")
    plt.xlabel('x', fontsize = 13)
    plt.ylabel('time, t', fontsize = 13)
    plt.title(r'Hovmüller diagram of $\psi(x, t)$, fwd Euler, basin')
 
    plt.figure(4)
    plt.style.use("ggplot")
    fig = plt.figure(figsize = (9,7))
    CS = plt.contourf(x, t, out_c_basin[1:,:].transpose(), 20, cmap = plt.cm.RdBu_r)
    plt.colorbar(CS, orientation = "vertical")
    plt.xlabel('x', fontsize = 13)
    plt.ylabel('time, t', fontsize = 13)
    plt.title(r'Hovmüller diagram of $\psi(x, t)$, centered, basin')
    
    plt.figure(2)
    plt.style.use("ggplot")
    fig = plt.figure(figsize = (9,7))
    CS = plt.contourf(x, t, out_periodic[1:,:].transpose(), 20, cmap = plt.cm.RdBu_r)
    plt.colorbar(CS, orientation = "vertical")
    plt.xlabel('x', fontsize = 13)
    plt.ylabel('time, t', fontsize = 13)
    plt.title(r'Hovmüller diagram of $\psi(x, t)$, fwd Euler, periodic')
    
    plt.figure(3)
    plt.style.use("ggplot")
    fig = plt.figure(figsize = (9,7))
    CS = plt.contourf(x, t, out_c_periodic[1:,:].transpose(), 20, cmap = plt.cm.RdBu_r)
    plt.colorbar(CS, orientation = "vertical")
    plt.xlabel('x', fontsize = 13)
    plt.ylabel('time, t', fontsize = 13)
    plt.title(r'Hovmüller diagram of $\psi(x, t)$, centered, periodic')
    
