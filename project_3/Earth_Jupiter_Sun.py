#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 19 10:38:33 2018

@author: maritkollstuen
"""

import numpy as np
# import matplotlib.pyplot as plt
import time


class planets:   
    def __init__(self, vel0, pos0, mass):
        self.vel0 = vel0        # AU/year
        self.pos0 = pos0        # position in astronomical units
        self.mass = mass/2.10E30
        self.G = 4*np.pi**2
        
    def kinetic_energy(self, vel):
        KE =  0.5*self.mass*np.linalg.norm(vel,axis=1)**2
        return KE
    
    def potential_energy(self, r, mass_of_sun=1):
        PE = -self.G*mass_of_sun*self.mass/np.linalg.norm(r, axis=1)
        return PE
    
    def angular_momentum(self, vel, r):
        L = np.cross(r, vel*self.mass,axis=1)
        return L
    
    def make_pos_vec(self, N):
        self.pos_vec = np.zeros((N,3))
        self.pos_vec[0] = self.pos0

    def make_vel_vec(self,N):
        self.vel_vec = np.zeros((N,3))
        self.vel_vec[0] = self.vel0

    def make_acc_vec(self,N,acc0):
        self.acc_vec = np.zeros((N,3))
        self.acc_vec[0] = acc0
    
 
    
def system(*args):
    system = []
    for i in args:
        system.append(i)

    return system
    
    
class solver: 
    
    
    def __init__(self, system, N, dt):
        self.dt = dt
        self.N = N
        self.system = system
        self.G = 4*np.pi**2
        
        for i in self.system: 
            i.make_pos_vec(N)
            i.make_vel_vec(N)
        for i in self.system:
            i.make_acc_vec(N, self.acceleration(i, 0))
    
    def acceleration(self, current, i):
        acc = np.zeros(3)
        
        for planet in self.system:
            if planet != current:
                r = current.pos_vec[i] - planet.pos_vec[i]
                acc -= r*self.G*planet.mass/(np.linalg.norm(r)**3)
        
        return acc

        
    def euler_fwd(self):
        dt = self.dt

        start = time.time()

        for i in range(self.N-1):
            for planet in self.system:
                planet.pos_vec[i+1] = planet.pos_vec[i] + dt*planet.vel_vec[i-1]
            for planet in self.system:    
                planet.acc_vec[i+1] = self.acceleration(planet, i+1)
            for planet in self.system:
                planet.vel_vec[i+1] = planet.vel_vec[i] + dt*planet.acc_vec[i]
        
        stop = time.time() 
        print('Euler: {:.3f} sec'.format(stop - start))
    
    
    def velocity_verlet(self):

        dt = self.dt
        dt05 = 0.5*dt
        dt2 = dt**2
        dt205 = 0.5*dt2

        
        start = time.time()

        for i in range(self.N-1):
            for planet in self.system:
                planet.pos_vec[i+1] = planet.pos_vec[i] + dt*planet.vel_vec[i] + dt205*planet.acc_vec[i]
            for planet in self.system:
                planet.acc_vec[i+1] = self.acceleration(planet, i+1)
            for planet in self.system:
                planet.vel_vec[i+1] = planet.vel_vec[i] + dt05*(planet.acc_vec[i] + planet.acc_vec[i+1])
        stop = time.time()
        print("Verlet: {:.3f} sec".format(stop-start))
    
if __name__ == "__main__":
    pass
    
    
    
    
    