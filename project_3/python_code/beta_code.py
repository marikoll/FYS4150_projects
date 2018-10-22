#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
To run, use python master_document.py

Project 3, task 3d,
changing the initial algorithm to include changes of the gravitational force.
The solver class now needs beta as input, and the acceleration can be changed easily.

"""


import numpy as np
import time


class planets:
    """
    Class that creates the planet-object
    
    Input: 
        vel0    - <numpy array> initial velocity of the planet in AU/year
        pos0    - <numpy array> initial position of the planet in AU
        mass    - <float> mass of planet
    """
    def __init__(self, vel0, pos0, mass):
        self.vel0 = vel0        
        self.pos0 = pos0        
        self.mass = mass/2.10E30
        self.G = 4*np.pi**2

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
        b       - <float> adjustment of gravitational force
    """
    def __init__(self, system, N, dt,b):
        self.dt = dt
        self.N = N
        self.system = system
        self.G = 4*np.pi**2
        self.b = b

        for i in self.system:
            i.make_pos_vec(N)
            i.make_vel_vec(N)
        for i in self.system:
            i.make_acc_vec(N, self.acceleration(i, 0, b))

    def acceleration(self, current, i, b):
        """
        Calculates the acceleration of the planet by adjusting the 
        gravitational force due to the presence of other planets
        
        Input: 
            current - <class object of planets> current planet
            i       - <int> index
            b       - <float> adjustment of gravitational force
        Output: 
            acc     - <numpy array> acceleration of planet    
        """
        acc = np.zeros(3)

        for planet in self.system:
            if planet != current:
                r = current.pos_vec[i] - planet.pos_vec[i]
                acc -= r*self.G*planet.mass/(np.linalg.norm(r)**b)

        return acc


    def velocity_verlet(self):
        """
        Calculates the trajectory of the planet by using the velocity Verlet
        method
        """

        b = self.b
        dt = self.dt
        dt2 = dt**2


        start = time.time()

        for i in range(self.N-1):
            for planet in self.system:
                planet.pos_vec[i+1] = planet.pos_vec[i] + dt*planet.vel_vec[i]\
                + 0.5*dt2*planet.acc_vec[i]
            for planet in self.system:
                planet.acc_vec[i+1] = self.acceleration(planet, i+1, b)
            for planet in self.system:
                planet.vel_vec[i+1] = planet.vel_vec[i] + \
                0.5*dt*(planet.acc_vec[i] + planet.acc_vec[i+1])
        stop = time.time()
        print("Verlet: {:.3f} sec".format(stop-start))


if __name__ == "__main__":
    pass
