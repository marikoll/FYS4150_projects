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
temp = np.linspace(0.1, 3.0, 101)
temp = 1.0

sampled_energies = np.zeros(len(trials))
sampled_magnets = np.zeros(len(trials))
sampled_cv = np.zeros(len(trials))
sampled_suscept = np.zeros(len(trials))
sampled_absmagn = np.zeros(len(trials))


num_trials = 5
start = time.time()


avg_energy = 0
avg_magnet = 0
avg_cv = 0
avg_suscept = 0
avg_absmagn = 0
j = 0

while j < num_trials:
    for i in range(len(trials)):
        grid = np.ones((spins, spins))
        energy_avg, magnet_avg, C_v, susceptibility, abs_magnet = MC(grid, trials[i], temp)#, w)
        sampled_energies[i] = energy_avg
        sampled_magnets[i] = magnet_avg
        sampled_cv[i] = C_v
        sampled_suscept[i] = susceptibility
        sampled_absmagn[i] = abs_magnet
    avg_energy += sampled_energies
    avg_magnet += sampled_magnets
    avg_cv += sampled_cv
    avg_suscept += sampled_suscept
    avg_absmagn += sampled_absmagn
    j += 1

avg_energy = avg_energy/num_trials
avg_magnet = np.abs(avg_magnet)/num_trials
avg_cv = avg_cv/num_trials
avg_suscept = avg_suscept/num_trials
avg_absmagn = avg_absmagn/num_trials
        
stop = time.time()
print('CPU: ', stop-start)  

plt.figure()
ax1 = plt.subplot(211)
ax1.set_title('Expectation values for energy and magnetization as functions of MC-cycles', fontsize = 12)
ax1.plot(trials, avg_energy, 'r')
ax1.set_ylabel(r'$\langle E \rangle$', fontsize = 10)
ax1.grid()
ax2 = plt.subplot(212)
ax2.plot(trials, avg_magnet)
ax2.set_ylabel(r'$\langle M\rangle$', fontsize = 10)
ax2.set_xlabel('MC-cycles', fontsize = 10)
ax2.grid()
plt.savefig('expectations_22lattice.pdf', bbox_inches = 'tight')
plt.show()
