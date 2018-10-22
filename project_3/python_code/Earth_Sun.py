#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
To run, use python master_document.py

Project 3, task 3b,

Module that calculates the trajectory of all planets defined in the system
by using the velocity Verlet method and forward Euler method and using the 
Sun as center of mass
"""

import numpy as np
import matplotlib.pyplot as plt
import time


class planets(object):
    """
    Class that creates the planet-object
    
    Input: 
        vel0    - <numpy array> initial velocity of the planet in AU/year
        pos0    - <numpy array> initial position of the planet in AU
        mass    - <float> mass of planet
    """       
        
    def __init__(self, vel0, pos0, mass):
        self.vel0 = vel0.copy()        # AU/year
        self.pos0 = pos0.copy()        # position in astronomical units
        self.mass = mass/2.10E30
        self.G = 4*np.pi**2
    
    def kinetic_energy(self, vel):
        """
        Calculates the kinetic energy of the planet
        
        Input: 
            vel     - <numpy array> velocity of the planet
            
        Output: 
            KE      - <numpy array> kinetic energy of the planet
        """
        KE =  0.5*self.mass*np.linalg.norm(vel,axis=1)**2
        return KE
    
    def potential_energy(self, pos, mass_of_sun=1):
        """
        Calculates the potential energy of the planet
        
        Input: 
            pos         - <numpy array> position of the planet
            mass_of_sun - <float> mass of the Sun, default = 1 (scaled)
            
        Output: 
            PE      - <numpy array> potential energy of the planet
        """        
        PE = -self.G*mass_of_sun*self.mass/np.linalg.norm(pos, axis=1)
        return PE
    
    def angular_momentum(self, vel, pos):
        """
        Calculates the angular momentum of the planet
        
        Input: 
            vel     - <numpy array> velocity of the planet
            pos     - <numpy array> position of the planet
            
        Output: 
            L       - <numpy array> angular momentum of the planet
        """
        L = np.cross(pos, vel*self.mass,axis=1)
        return L
    
    def make_pos_vec(self, N):
        """
        Creates a position vector of size (N,3) with pos0 as first entry
        
        Input: 
            N       - <int> length of vector
        """        
        self.pos_vec = np.zeros((N,3))
        self.pos_vec[0] = self.pos0

    def make_vel_vec(self,N):
        """
        Creates a velocity vector of size (N,3) with vel0 as first entry
        
        Input: 
            N       - <int> length of vector
        """        
        self.vel_vec = np.zeros((N,3))
        self.vel_vec[0] = self.vel0

    def make_acc_vec(self,N,acc0):
        """
        Creates an acceleration vector of size (N,3) with acc0 as first entry
        
        Input: 
            N       - <int> length of vector
            acc0    - <numpy array> initial acceleration of the planet
        """        
        self.acc_vec = np.zeros((N,3))
        self.acc_vec[0] = acc0

    
def system(*args):
    """
    Function that generates a system of all input planets
    
    Input: 
        *args   - <class object of planets> planets to be used in further 
                  calculations

    Output: 
        system  - <list> all planets in the system
    """  
    system = []
    for i in args:
        system.append(i)

    return system
    

class solver: 
    """
    Class that calculates the trajectory of a planet by applying the velocity
    Verlet method
    
    Input: 
        system  - <list> all planets in the system
        N       - <int> length of list
        dt      - <float> step size
    """       
    
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
        """
        Calculates the acceleration of the planet by adjusting the 
        gravitational force due to the presence of other planets
        
        Input: 
            current - <class object of planets> current planet
            i       - <int> index
        Output: 
            acc     - <numpy array> acceleration of planet
        """               
        
        acc = np.zeros(3)
        
        for planet in self.system:
            if planet != current:
                r = current.pos_vec[i] - planet.pos_vec[i]
                acc -= r*self.G*planet.mass/(np.linalg.norm(r)**3)
        
        return acc
    
    
    
    def euler_fwd(self):
        """
        Calculates the trajectory of the planet by using the forward Euler
        method
        """
        dt = self.dt

        start = time.time()

        for i in range(self.N-1):
            for planet in self.system:
                planet.pos_vec[i+1] = planet.pos_vec[i] + dt*planet.vel_vec[i]
            for planet in self.system:    
                planet.acc_vec[i+1] = self.acceleration(planet, i+1)
            for planet in self.system:
                planet.vel_vec[i+1] = planet.vel_vec[i] + dt*planet.acc_vec[i]
        
        stop = time.time() 
        print('Euler: {:.3f} sec'.format(stop - start))

    def velocity_verlet(self):
        """
        Calculates the trajectory of the planet by using the velocity Verlet
        method
        """    
        dt = self.dt 
        dt05 = 0.5*dt
        dt205 = 0.5*dt**2

        start = time.time()
        for i in range(self.N-1):
            for planet in self.system:
                planet.pos_vec[i+1] = planet.pos_vec[i] + dt*planet.vel_vec[i]\
                + dt205*planet.acc_vec[i]
            for planet in self.system:
                planet.acc_vec[i+1] = self.acceleration(planet, i+1)
            for planet in self.system:
                planet.vel_vec[i+1] = planet.vel_vec[i] + \
                dt05*(planet.acc_vec[i] + planet.acc_vec[i+1])
        stop = time.time() 
        print('Verlet: {:.3f} sec'.format(stop - start))    

if __name__ == "__main__":
    pass