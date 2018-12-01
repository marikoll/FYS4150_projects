#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 30 15:09:36 2018

@author: maritkollstuen
"""

#Program which solves the 2+1-dimensional wave equation by a finite difference scheme

import numpy as np
#Define the grid
N = 31
h = 1.0 / (N-1)
dt = .0005
t_steps = 10000

x1 = np.linspace(0,1,N)
y1 = np.linspace(0,1,N)
x,y = np.meshgrid(x1,y1,sparse=False)

alpha = dt**2 / h**2

#Initial conditions with du/dt = 0
u = np.sin(x*np.pi)*np.cos(y*np.pi-np.pi/2)
u_old = np.zeros(u.shape,type(u[0,0]))
for i in range(1,N-1):
    for j in range(1,N-1):
        u_old[i,j] = u[i,j] + (alpha/2)*(u[i+1,j] - 4*u[i,j] + u[i-1,j] + u[i,j+1] + u[i,j-1])
u_new = np.zeros(u.shape,type(u[0,0]))

#We don't necessarily want to plot every time step. We plot every n'th step where
n = 100
plotnr = 0

#Iteration over time steps
for k in range(t_steps):
    for i in range(1,N-1): #1 - N-2 because we don't want to change the boundaries
        for j in range(1,N-1):
            u_new[i,j] = 2*u[i,j] - u_old[i,j] + alpha*(u[i+1,j] - 4*u[i,j] + u[i-1,j] + u[i,j+1] + u[i,j-1])

    #Prepare for next time step by manipulating pointers
    temp = u_new
    u_new = u_old
    u_old = u
    u = temp