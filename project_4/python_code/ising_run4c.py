#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Nov 11 16:55:29 2018

@author: maritkollstuen
"""

import numpy as np
import time
import matplotlib.pyplot as plt

from Ising_model import MC


spins       = 20
trials      = np.arange(100, 1e6, 5000, np.dtype(np.int64))

temp = [1.0, 2.4]

sampled_energies = np.zeros((len(trials), len(temp)))
sampled_magnets = np.zeros((len(trials), len(temp)))
sampled_cv = np.zeros((len(trials), len(temp)))
sampled_suscept = np.zeros((len(trials), len(temp)))
sampled_absmagn = np.zeros((len(trials), len(temp)))

start = time.time()


for k in range(len(temp)):
    for i in range(len(trials)):
        grid = np.random.choice([-1,1],size=(spins, spins))
        energy_avg, magnet_avg, C_v, susceptibility, abs_magnet, c = MC(grid, trials[i], temp[k])
        sampled_energies[i, k] = energy_avg
        sampled_magnets[i, k] = magnet_avg
        sampled_cv[i, k] = C_v
        sampled_suscept[i, k] = susceptibility
        sampled_absmagn[i, k] = abs_magnet
        
o_sampled_energies = np.zeros((len(trials), len(temp)))
o_sampled_magnets = np.zeros((len(trials), len(temp)))
o_sampled_cv = np.zeros((len(trials), len(temp)))
o_sampled_suscept = np.zeros((len(trials), len(temp)))
o_sampled_absmagn = np.zeros((len(trials), len(temp)))



for k in range(len(temp)):
    for i in range(len(trials)):
        grid = np.ones((spins, spins))
        energy_avg, magnet_avg, C_v, susceptibility, abs_magnet, c = MC(grid, trials[i], temp[k], ordered = True)
        o_sampled_energies[i, k] = energy_avg
        o_sampled_magnets[i, k] = np.abs(magnet_avg)
        o_sampled_cv[i, k] = C_v
        o_sampled_suscept[i, k] = susceptibility
        o_sampled_absmagn[i, k] = abs_magnet

stop = time.time()
print('CPU: ', stop-start)  

plt.figure()
ax1 = plt.subplot(211)
ax1.set_title('Expectation values for energy and magnetization as functions of MC-cycles, random configuration', fontsize = 12)
ax1.plot(trials, sampled_energies[:,0], 'r', label= 'T = 1.0')
ax1.plot(trials, sampled_energies[:,1], 'b', label = 'T = 2.4')
ax1.legend()
ax1.set_ylabel(r'$\langle E \rangle$', fontsize = 10)
ax1.grid()
ax2 = plt.subplot(212)
ax2.plot(trials, sampled_absmagn[:,0], 'r', label = 'T = 1.0')
ax2.plot(trials, sampled_absmagn[:,1], 'b', label = 'T = 2.4')
ax2.legend()
ax2.set_ylabel(r'$\langle |M|\rangle$', fontsize = 10)
ax2.set_xlabel('MC-cycles', fontsize = 10)
ax2.grid()
plt.savefig('figs/expectations_2020lattice_random.pdf', bbox_inches = 'tight')
plt.show()



plt.figure()
ax1 = plt.subplot(211)
ax1.set_title('Expectation values for energy and magnetization as functions of MC-cycles, ordered configuration', fontsize = 12)
ax1.plot(trials, o_sampled_energies[:,0], 'r', label= 'T = 1.0')
ax1.plot(trials, o_sampled_energies[:,1], 'b', label = 'T = 2.4')
ax1.legend()
ax1.set_ylabel(r'$\langle E \rangle$', fontsize = 10)
ax1.grid()
ax2 = plt.subplot(212)
ax2.plot(trials, o_sampled_absmagn[:,0], 'r', label = 'T = 1.0')
ax2.plot(trials, o_sampled_absmagn[:,1], 'b', label = 'T = 2.4')
ax2.legend()
ax2.set_ylabel(r'$\langle |M|\rangle$', fontsize = 10)
ax2.set_xlabel('MC-cycles', fontsize = 10)
ax2.grid()
plt.savefig('figs/expectations_2020lattice_ordered.pdf', bbox_inches = 'tight')
plt.show()

