#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Nov 11 14:15:51 2018

@author: maritkollstuen
"""

import numpy as np
import time
import matplotlib.pyplot as plt

from Ising_model import MC


spins       = 2
trials      = [int(1e2), int(1e3), int(1e4), int(1e5), int(1e6), int(1e7)]
temp = 1.0

sampled_energies = np.zeros(len(trials))
sampled_magnets = np.zeros(len(trials))
sampled_cv = np.zeros(len(trials))
sampled_suscept = np.zeros(len(trials))
sampled_absmagn = np.zeros(len(trials))


start = time.time()


for i in range(len(trials)):
    grid = np.ones((spins, spins))
    energy_avg, magnet_avg, C_v, susceptibility, abs_magnet, c= MC(grid, trials[i], temp)#, w)
    sampled_energies[i] = energy_avg
    sampled_magnets[i] = magnet_avg
    sampled_cv[i] = C_v
    sampled_suscept[i] = susceptibility
    sampled_absmagn[i] = abs_magnet


        
stop = time.time()
print('CPU: ', stop-start)  

plt.figure()
ax1 = plt.subplot(211)
ax1.set_title('Expectation values for energy and magnetization as functions of MC-cycles', fontsize = 12)
ax1.plot(trials, sampled_energies, 'r')
ax1.set_ylabel(r'$\langle E \rangle$', fontsize = 10)
ax1.grid()
ax2 = plt.subplot(212)
ax2.plot(trials, sampled_absmagn)
ax2.set_ylabel(r'$\langle |M|\rangle$', fontsize = 10)
ax2.set_xlabel('MC-cycles', fontsize = 10)
ax2.grid()
plt.savefig('figs/expectations_22lattice.pdf', bbox_inches = 'tight')
plt.show()



