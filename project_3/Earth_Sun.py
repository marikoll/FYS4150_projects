#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 12 10:18:33 2018
@author: maritkollstuen
"""

import numpy as np
import matplotlib.pyplot as plt
import time


class planets(object):
    
        
    def __init__(self, vel0, pos0, mass):
        self.vel0 = vel0.copy()        # AU/year
        self.pos0 = pos0.copy()        # position in astronomical units
        self.mass = mass/2.10E30
        self.G = 4*np.pi**2
#        self.N = N
#        self.dt = dt

        

    def kinetic_energy(self, vel):
        KE =  0.5*self.mass*np.linalg.norm(vel,axis=1)**2
        return KE
    
    def potential_energy(self, r, mass_of_sun=1):
        PE = -self.G*mass_of_sun*self.mass/np.linalg.norm(r, axis=1)#self.distance(pos)
        return PE
    
    def angular_momentum(self, vel, r):
        L = np.cross(r, vel*self.mass)#,axis=1)
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
#    def distance(self, pos):
#        self.r = np.linalg.norm(pos)
#        return self.r
#        
#        
#    def acceleration(self, pos):
#        r = self.distance(self.pos)
#        self.acc = (-pos*(self.GMs)/(r**3))
#        
#        return self.acc
    
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
#        t_f = 1
#        t = 0
#        h = np.linspace(t, t_f, self.N*t_f)
#        dt = h[1]-h[0]  
        dt = self.dt


        
#        pos_e = np.empty((N,3))
#        pos_e[0,:] = self.pos
#        
#        vel_e = np.empty((N,3))
#        vel_e[0,:] = self.vel
#        
#        acc_e = np.empty((N,3))
#        acc_e[0,:] = self.acceleration(self.pos)
#        
        start = time.time()
#        for i in range(1, N):
#            pos_e[i] = pos_e[i-1] + dt*vel_e[i-1]
#            acc_e[i] = self.acceleration(pos_e[i])
#            vel_e[i] = vel_e[i-1] + dt*acc_e[i-1]
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

#        t_f = 1
#        t = 0
#        h = np.linspace(t, t_f, self.N*t_f)
#        dt = h[1]-h[0]
        dt = self.dt 
        dt05 = 0.5*dt
        dt205 = 0.5*dt**2
        # N = self.N
        
#        pos_v = np.empty((N,3))
#        pos_v[0,:] = self.pos
#        
#        vel_v = np.empty((N,3))
#        vel_v[0,:] = self.vel
#        
#        acc_v = np.empty((N,3))
#        acc_v[0,:] = self.acceleration(self.pos)
#        
#        start = time.time()
#        for i in range(1, N):
#            pos_v[i] = pos_v[i-1] + dt*vel_v[i-1] + dt205*acc_v[i-1]
#            acc_v[i] = self.acceleration(pos_v[i])
#            vel_v[i] = vel_v[i-1] + dt05*(acc_v[i] + acc_v[i-1])
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

    
#    def potential_energy(self, method):
#        if method == 'v':
#            pos, acc, vel, dt = self.velocity_verlet()
#        elif method == 'e':
#            pos, acc, vel, dt = self.euler_fwd()
#        r_1 = self.distance(pos[0,:])
#        r_2 = self.distance(pos[-1,:])
#        PE_1 = -self.gravit_force*(self.mass_of_sun*self.mass)/r_1
#        PE_2 = -self.gravit_force*(self.mass_of_sun*self.mass)/r_2
#        return PE_1, PE_2 
#    
#    def kinetic_energy(self, method):
#        if method == 'v':
#            pos, acc, vel, dt = self.velocity_verlet()
#        elif method == 'e':
#            pos, acc, vel, dt = self.euler_fwd()
#        
#        KE_1 = np.linalg.norm(0.5*self.mass*vel[0,:]**2)
#        KE_2 = np.linalg.norm(0.5*self.mass*vel[-1,:]**2)
#        return KE_1, KE_2
#    
#    
#    def angular_momentum(self, method):
#        if method == 'v':
#            pos, acc, vel, dt = self.velocity_verlet()
#        elif method == 'e':
#            pos, acc, vel, dt = self.euler_fwd()
#        r_1 = self.distance(pos[0,:])
#        r_2 = self.distance(pos[-1,:])
#        L1 = np.linalg.norm(self.mass*vel[0,:]*r_1)
#        L2 = np.linalg.norm(self.mass*vel[-1,:]*r_2)
#        return L1, L2
#    
#    def test_conservation(self, method):
#        err = 1E-5   
#        PE_1, PE_2 = self.potential_energy(method)
#        KE_1, KE_2 = self.kinetic_energy(method)
#        L_1, L_2 = self.angular_momentum(method)
#
#        assert abs(PE_1 - PE_2) < err
#        assert abs(KE_1 - KE_2) < err
#        assert abs(L_1 - L_2) < err
#        print('Potential and kintetic energy and angular momentum is conserved')
        
#    def make_dict(self):
#        self.planetDict = {self.planet: }

def plot_figures_3c(pos_v, pos_e, dt_v, dt_e, N):
    
    
    plt.figure()
    plt.subplot(121)
    plt.axis('equal')
    plt.plot(pos_v[:,0], pos_v[:,1], label = 'dt = {:.4}\nN = {}'.format(dt_v, N))
    plt.legend()
    plt.title('Verlet method')
    
    plt.subplot(122)
    plt.axis('equal')
    plt.plot(pos_e[:,0], pos_e[:,1], label = 'dt = {:.4}\nN = {}'.format(dt_e, N))
    plt.legend()
    plt.title('Forward Euler')
    plt.savefig('euler_verlet_N{}.pdf'.format(i))


    

if __name__ == "__main__":
    
    planetE = 'Earth'
    velE = [-5.71E-3, 1.62E-2, -3.03E-7]
    posE = np.array([9.47E-1, 3.21E-1, -9.31E-5])
    massE = 6E24
    N = 7
    dt = 0.1677
    Earth = planets(velE, posE, massE, N, dt)

    pos_v = Earth.velocity_verlet()
    pos_e = Earth.euler_fwd()
    
    plot_figures_3c(pos_v, pos_e, dt, dt, N)
    
#    kin_energy1, kin_energy2 = Earth.kinetic_energy('e')
#    pot_energy1, pot_energy2 = Earth.potential_energy('e')
#    L1, L2 = Earth.angular_momentum('e')
#    Earth.test_conservation('v')
    
    
    ## Plot difference in Verlet and Euler (task 3c)
#    N = [7, 10, 50, 100, 220, 1000, 10000]
#    for i in N:
#        Earth = planets(velE, posE, massE, planetE, i)
#        F_g = Earth.acceleration(posE)
#        pos_v, acc_v, vec_v, dt_v = Earth.velocity_verlet()
#        pos_e, acc_e, vec_e, dt_e = Earth.euler_fwd()
#        plot_figures_3c(pos_v, pos_e, dt_v, dt_e, N)

        
    ##### OTHER PLANETS
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