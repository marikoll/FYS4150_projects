#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 12 10:18:33 2018

@author: maritkollstuen
"""

import numpy as np
import matplotlib.pyplot as plt


class planets(object):
    
        
    def __init__(self, vel, pos, mass, planet, N):
        self.vel = [i*365 for i in vel] # AU/year
        self.pos = np.array(pos)        # position in astronomical units
        self.mass = mass/2.10E30
        self.planet = planet
        self.mass_of_sun = 1
        self.pos_sun = [0,0,0]
        self.vel_sun = [0,0,0]
        self.gravit_force = 4*np.pi**2
        self.N = N
        self.GMs = self.gravit_force*self.mass_of_sun
        


    def acceleration(self, pos):
        
        r = np.linalg.norm(pos)
        
        self.acc = (-pos*(self.GMs)/(r**3))
        
        return self.acc
    

    def euler_fwd(self):
        t_f = 1
        t = 0
        h = np.linspace(t, t_f, self.N*t_f)
        dt = h[1]-h[0]  
        
        pos = np.empty((len(h),3))
        pos[0,:] = self.pos
        
        vel = np.empty((len(h),3))
        vel[0,:] = self.vel
        
        acc = np.empty((len(h),3))
        acc[0,:] = self.acceleration(self.pos)
        
        
        for i in range(1, len(h)):
            pos[i] = pos[i-1] + dt*vel[i-1]
            acc[i] = self.acceleration(pos[i])
            vel[i] = vel[i-1] + dt*acc[i-1]
            
        return pos, acc, vel, dt
        
    def velocity_verlet(self):

        t_f = 1
        t = 0
        h = np.linspace(t, t_f, self.N*t_f)
        dt = h[1]-h[0]
        
        pos = np.empty((len(h),3))
        pos[0,:] = self.pos
        
        vel = np.empty((len(h),3))
        vel[0,:] = self.vel
        
        acc = np.empty((len(h),3))
        acc[0,:] = self.acceleration(self.pos)
        
        
        for i in range(1, len(h)):
            pos[i] = pos[i-1] + dt*vel[i-1] + (dt**2/2)*acc[i-1]
            acc[i] = self.acceleration(pos[i])
            vel[i] = vel[i-1] + (dt/2)*(acc[i] + acc[i-1])
            
        return pos, acc, vel, dt
    
#    def make_dict(self):
#        self.planetDict = {self.planet: }

def plot_figures_3c(pos_v, pos_e, dt_v, dt_e, N):
    
    
    plt.figure()
    plt.subplot(121)
    plt.axis('equal')
    plt.plot(pos_v[:,0], pos_v[:,1], label = 'dt = {:.4}\nN = {}'.format(dt_v, i))
    plt.legend()
#    axes = plt.gca()
#    axes.set_ylim([-1,1])
#    axes.set_xlim([-1.15,1.15])
    plt.title('Verlet method')#\ndt = {:.4f}\nN = {}'.format(N))#, fontsize = 15)
    
    plt.subplot(122)
    plt.axis('equal')
    plt.plot(pos_e[:,0], pos_e[:,1], label = 'dt = {:.4}\nN = {}'.format(dt_e, i))
    plt.legend()
#    axes = plt.gca()
#    axes.set_ylim([-1,1])
#    axes.set_xlim([-1.15,1.15])
    plt.title('Forward Euler')# method\ndt = {:.4f}\nN = {}'.format(dt_e, N))#, fontsize = 15)
    
#       plt.plot(pos_J[:,0], pos_J[:,1])
#       plt.plot(pos_m[:,0], pos_m[:,1])
    plt.savefig('euler_verlet_N{}.pdf'.format(i))

if __name__ == "__main__":
    
    planetE = 'Earth'
    velE = [-5.71E-3, 1.62E-2, -3.03E-7]
    posE = np.array([9.47E-1, 3.21E-1, -9.31E-5])
    massE = 6E24
#    N = 50
#    Earth = planets(velE, posE, massE, planetE, N)
#    F_g = Earth.acceleration(posE)
#    pos_v, acc_v, vec_v, dt_v = Earth.velocity_verlet()
#    pos_e, acc_e, vec_e, dt_e = Earth.euler_fwd()
    
    ### Plot difference in Verlet and Euler (task 3c)
    N = [10, 50, 100, 1000, 10000]
    for i in N:
        Earth = planets(velE, posE, massE, planetE, i)
        F_g = Earth.acceleration(posE)
        pos_v, acc_v, vec_v, dt_v = Earth.velocity_verlet()
        pos_e, acc_e, vec_e, dt_e = Earth.euler_fwd()
        plot_figures_3c(pos_v, pos_e, dt_v, dt_e, N)
#    
#    planetJ = 'Jupiter'
#    velJ = [6.45E-3, -3.4E-3, -1.3E-4]
#    posJ = np.array([-2.67, -4.65, 7.9E-2])
#    massJ = 1898.13E24
#    N_j = 100
#    Jupiter = planets(velJ, posJ, massJ, planetJ, N_j)
#    pos_J, acc_j, vec_j = Jupiter.velocity_verlet()
#    
#    planetM = 'Mercury'
#    velM = [1.98E-1, -1.04E-2, -2.67E-3]
#    posM = np.array([-1.94E-1, -4.13E-1, -1.67E-2])
#    massM = 3.302E23
#    N_m = 50
#    Mercury = planets(velM, posM, massM, planetM, N_m)
#    pos_m, acc_m, vel_m = Mercury.velocity_verlet()
    

#    earth = planets(1, 0, 0, 0, 1.2, 0, 6E24)
#    r = earth.length()
#    f = earth.force()
    
    
    
    
    
    