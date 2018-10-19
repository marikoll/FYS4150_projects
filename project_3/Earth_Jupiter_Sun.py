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
#        self.mass_of_sun = 1
#        self.pos_sun = [0,0,0]
#        self.vel_sun = [0,0,0]

#    def distance(self, pos):
#        self.r = np.linalg.norm(pos,axis=1)
#        return self.r
        
    def kinetic_energy(self, vel):
        KE =  0.5*self.mass*np.linalg.norm(vel,axis=1)**2
        return KE
    
    def potential_energy(self, r, mass_of_sun=1):
        PE = -self.G*mass_of_sun*self.mass/np.linalg.norm(r, axis=1)#self.distance(pos)
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

        
    def velocity_verlet(self):
#
#        t_f = 1
#        t = 0
#        h = np.linspace(t, t_f, self.N*t_f)
#        dt = h[1]-h[0]
        dt = self.dt
        dt2 = dt**2

        
        start = time.time()

        for i in range(self.N-1):
            for planet in self.system:
                planet.pos_vec[i+1] = planet.pos_vec[i] + dt*planet.vel_vec[i]\
                + 0.5*dt2*planet.acc_vec[i]
            for planet in self.system:
                planet.acc_vec[i+1] = self.acceleration(planet, i+1)
            for planet in self.system:
                planet.vel_vec[i+1] = planet.vel_vec[i] + \
                0.5*dt*(planet.acc_vec[i] + planet.acc_vec[i+1])
        stop = time.time()
        print("Verlet: {:.3f} sec".format(stop-start))
    
if __name__ == "__main__":
    
#    velE = np.array([-5.71E-3, 1.62E-2, -3.03E-7])
#    posE = np.array([9.47E-1, 3.21E-1, -9.31E-5])
#    massE = 6E24
#    N = 220
#    Earth = planets(velE, posE, massE)
#
#
#    velJ = np.array([6.45E-3, -3.4E-3, -1.3E-4])
#    posJ = np.array([-2.67, -4.65, 7.9E-2])
#    massJ = 1898.13E24
#    Jupiter = planets(velJ, posJ, massJ)
#
#
#    velM = np.array([1.98E-1, -1.04E-2, -2.67E-3])
#    posM = np.array([-1.94E-1, -4.13E-1, -1.67E-2])
#    massM = 3.302E23
#    Mercury = planets(velM, posM, massM)
#    massS = 2E30
#    posS = np.array([0,0,0])
#    velS = -(massE*velE + massJ*velJ + massM*velM)/massS
#    Sun = planets(velS, posS, massS)
    
    pass
    
    
    
    
    