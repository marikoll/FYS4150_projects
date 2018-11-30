#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Nov 11 18:47:20 2018

@author: maritkollstuen
"""

import numpy as np
import time
import matplotlib.pyplot as plt
from numba import prange

from Ising_model import MC


spins       = [40, 60, 80, 100]
trials      = int(1e7)

temp = np.arange(2.1, 2.5, 0.05)

temp = np.linspace(2.1,2.5,8)

sampled_energies = np.zeros((len(spins), len(temp)))
sampled_cv = np.zeros((len(spins), len(temp)))
sampled_suscept = np.zeros((len(spins), len(temp)))
sampled_absmagn = np.zeros((len(spins), len(temp)))

start = time.time()
p = 0.1 # take away 10 percent

for k in prange(len(temp)):
    for i in range(len(spins)):
        grid = np.ones((spins[i], spins[i]))
        energy_avg, magnet_avg, C_v, susceptibility, abs_magnet, c = MC(grid, trials, temp[k])
        sampled_energies[i, k] = energy_avg[-1]
        sampled_cv[i, k] = C_v
        sampled_suscept[i, k] = susceptibility
        sampled_absmagn[i, k] = abs_magnet[-1]
        
        

        
stop = time.time()
print('CPU: ', stop-start)  



plt.figure(1)
for i in range(len(spins)):
    plt.plot(temp, sampled_energies[i, :], label = 'L{}'.format(spins[i]))
#plt.axvline(2.269, color='k', linestyle='--')
plt.legend()
plt.ylabel(r'$\langle E\rangle$')
plt.xlabel('Temperature')
plt.title('Expectation values of energy versus temperature')
plt.savefig('figs/Evs_energy_Test1.pdf')

plt.figure(2)
for i in range(len(spins)):
    plt.plot(temp, sampled_absmagn[i, :], label = 'L{}'.format(spins[i]))
#plt.axvline(2.269, color='k', linestyle='--')
plt.legend()
plt.ylabel(r'$\langle |M|\rangle$')
plt.xlabel('Temperature')
plt.title('Expectation values of magnetization versus temperature')
plt.savefig('figs/Evs_magn_Test1.pdf')

plt.figure(3)
for i in range(len(spins)):
    plt.plot(temp, sampled_cv[i, :], label = 'L{}'.format(spins[i]))
#plt.axvline(2.269, color='k', linestyle='--')
plt.legend()
plt.ylabel(r'$C_v$')
plt.xlabel('Temperature')
plt.title('Expectation values of heat capacity versus temperature')
plt.savefig('figs/Evs_cv_Test1.pdf')

plt.figure(4)
for i in range(len(spins)):
    plt.plot(temp, sampled_suscept[i, :], label = 'L{}'.format(spins[i]))
#plt.axvline(2.269, color='k', linestyle='--')
plt.legend()
plt.ylabel(r'$\chi$')
plt.xlabel('Temperature')
plt.title('Expectation values of susceptibility versus temperature')
plt.savefig('figs/Evs_suscept_Test1.pdf')
