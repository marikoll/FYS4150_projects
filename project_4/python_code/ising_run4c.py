#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Expectation values as funct. of MC-cycles
"""

import numpy as np
import time
import matplotlib.pyplot as plt

from Ising_model import MC


spins       = 20
trials      = int(1e6)# np.logspace(2, 7, 100, base = 10.0, dtype = np.dtype(np.int64))

temp = [1.0, 2.4]



sampled_energies = np.zeros((trials, len(temp)))
sampled_absmagn = np.zeros((trials, len(temp)))


start = time.time()


for k in range(len(temp)):
#    for i in range(len(trials)):
    grid = np.random.choice([-1,1],size=(spins, spins))
    energy_avg, magnet_avg, C_v, susceptibility, abs_magnet, c = MC(grid, trials, temp[k])#, cumsum = True)
    sampled_energies[:, k] = energy_avg
    sampled_absmagn[:, k] = abs_magnet

        

o_sampled_energies = np.zeros((trials, len(temp)))
o_sampled_absmagn = np.zeros((trials, len(temp)))


for k in range(len(temp)):
#    for i in range(len(trials)):
    grid = np.ones((spins, spins))
    energy_avg, magnet_avg, C_v, susceptibility, abs_magnet, c = MC(grid, trials, temp[k])#, cumsum = True)
    o_sampled_energies[:, k] = energy_avg
    o_sampled_absmagn[:, k] = abs_magnet

stop = time.time()
print('CPU: ', stop-start)  


T = np.linspace(0, trials, trials)


plt.figure()
ax1 = plt.subplot(211)
ax1.set_title('Expectation values for energy and magnetization as functions of MC-cycles, T = 1.0', fontsize = 12)
ax1.plot(T, o_sampled_energies[:,0], 'r', label= 'ordered configuration')
ax1.plot(T, sampled_energies[:,0], 'b', label = 'random configuration')
ax1.legend()
ax1.set_ylabel(r'$\langle E \rangle$', fontsize = 10)
ax1.grid()
ax2 = plt.subplot(212)
ax2.plot(T, o_sampled_absmagn[:,0], 'r', label = 'ordered configuration')
ax2.plot(T, sampled_absmagn[:,0], 'b', label = 'random configuration')
ax2.legend()
ax2.set_ylabel(r'$\langle |M|\rangle$', fontsize = 10)
ax2.set_xlabel('MC-cycles', fontsize = 10)
ax2.grid()
plt.savefig('figs/expectations_2020lattice_temp1.pdf', bbox_inches = 'tight')
plt.show()



plt.figure()
ax1 = plt.subplot(211)
ax1.set_title('Expectation values for energy and magnetization as functions of MC-cycles, T = 2.4', fontsize = 12)
ax1.plot(T, o_sampled_energies[:,1], 'r', label= 'ordered configuration')
ax1.plot(T, sampled_energies[:,1], 'b', label = 'random configuration')
ax1.legend()
ax1.set_ylabel(r'$\langle E \rangle$', fontsize = 10)
ax1.grid()
ax2 = plt.subplot(212)
ax2.plot(T, o_sampled_absmagn[:,1], 'r', label = 'ordered configuration')
ax2.plot(T, sampled_absmagn[:,1], 'b', label = 'random configuration')
ax2.legend()
ax2.set_ylabel(r'$\langle |M|\rangle$', fontsize = 10)
ax2.set_xlabel('MC-cycles', fontsize = 10)
ax2.grid()
plt.savefig('figs/expectations_2020lattice_temp24.pdf', bbox_inches = 'tight')
plt.show()