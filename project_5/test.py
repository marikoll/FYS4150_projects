#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 30 10:31:15 2018

@author: maritkollstuen
"""

import numpy as np
import matplotlib.pyplot as plt

def forward_step(alpha,u,uPrev,N):
    """
    Steps forward-euler algo one step ahead.
    Implemented in a separate function for code-reuse from crank_nicolson()
    """
    
    for x in range(1,N+1): #loop from i=1 to i=N
        u[x] = alpha*uPrev[x-1] + (1.0-2*alpha)*uPrev[x] + alpha*uPrev[x+1]


def forward_euler(alpha,u,N,T):
    """
    Implements the forward Euler sheme, results saved to
    array u
    """
    u_fwdE = u
    #Skip boundary elements
    for t in range(1,T):
        forward_step(alpha,u_fwdE[t],u_fwdE[t-1],N)
    
    return u_fwdE

def tridiag(alpha,u,N):
    """
    Tridiagonal gaus-eliminator, specialized to diagonal = 1+2*alpha,
    super- and sub- diagonal = - alpha
    """
    d = np.zeros(N) + (1-2*alpha)
    b = np.zeros(N-1) + alpha
    u_new = u
    #Forward eliminate
    for i in range(1,N):
        #Normalize row i (i in u convention):
        b[i-1] /= d[i-1];
        u_new[i] /= d[i-1] #Note: row i in u = row i-1 in the matrix
        d[i-1] = 1.0
        #Eliminate
        u_new[i+1] += u[i]*alpha
        d[i] += b[i-1]*alpha
    #Normalize bottom row
    u_new[N] /= d[N-1]
    d[N-1] = 1.0

    #Backward substitute
    for i in range(N,1,-1): #loop from i=N to i=2
        u_new[i-1] -= u_new[i]*b[i-2]
    return u_new    

       
def crank_nicolson(alpha,u,N,T):
    """
    Implents crank-nicolson scheme, reusing code from forward- and backward euler
    """
    
    u_CN  = u
    for t in range(1,T):
        forward_step(alpha/2,u_CN[t],u_CN[t-1],N)
        tridiag(alpha/2,u_CN[t],N)
    return u_CN


def psi(x):
    """Initial condition u(x,0) = g(x), x \in [0,1]"""
    
    return np.sin(4*np.pi*x)


if __name__ == "__main__":
    # Number of integration points along x-axis
    N       =   100
    # Step length in time

    # Number of time steps till final time 
    T       =   150

    #dx = 1/float(N+1)
    u = np.zeros((T,N+2),np.double)
    x = np.linspace (0,1,N+2)
    dx = 1/40
    
    dt= 0.001
    t = np.linspace(0, T, T+2)
    alpha = dt/(dx*2)
    
    #Initial codition
    u[0,:] = psi(x)
    u[0,0] = u[0,N+1] = 0.0 #Implement boundaries rigidly
    
    zeta = tridiag(dx**2, u, N)
    
    uE = forward_euler(alpha, u, N, T)
    uCN = crank_nicolson(alpha,u,N,T)
    
    
    
    
    plt.figure()
    plt.plot(x, uE[1,:], 'r-')
    plt.plot(x, uCN[1,:], 'b-.')

    