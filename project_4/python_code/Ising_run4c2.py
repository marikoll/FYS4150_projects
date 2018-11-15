#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 14 13:24:04 2018

@author: maritkollstuen
"""

import numpy as np
import time
import matplotlib.pyplot as plt

from Ising_model import MC_cutoff


spins       = 20
trials      = np.logspace(2, 7, 100, base = 10.0, dtype = np.dtype(np.int64)) #, np.dtype(np.int64))#[int(1e2), int(1e3), int(1e4), int(1e5), int(1e6), int(1e7)]

temp = [1.0, 2.4]

p = 0.1
trial_length = int(len(trials)*(1-p))
Temp = len(trials) - trial_length

sampled_energies = np.zeros((trial_length, len(temp)))
sampled_magnets = np.zeros((trial_length, len(temp)))
sampled_cv = np.zeros((trial_length, len(temp)))
sampled_suscept = np.zeros((trial_length, len(temp)))
sampled_absmagn = np.zeros((trial_length, len(temp)))
#sampled_count = np.zeros((int(9e5), len(temp)))



start = time.time()


for k in range(len(temp)):
    for i in range(len(trials)):
        grid = np.random.choice([-1,1],size=(spins, spins))
        energy_avg, magnet_avg, C_v, susceptibility, abs_magnet = MC_cutoff(grid, int(trials[i]), temp[k],  P = p)
        if i >Temp:
            sampled_energies[i - Temp, k] = energy_avg
            sampled_magnets[i - Temp, k] = np.abs(magnet_avg)
            sampled_cv[i - Temp, k] = C_v
            sampled_suscept[i - Temp, k] = susceptibility
            sampled_absmagn[i - Temp, k] = abs_magnet
#        sampled_count[:, k] = counter_list
        
        
#o_sampled_energies = np.zeros((trial_length, len(temp)))
#o_sampled_magnets = np.zeros((trial_length, len(temp)))
#o_sampled_cv = np.zeros((trial_length, len(temp)))
#o_sampled_suscept = np.zeros((trial_length, len(temp)))
#o_sampled_absmagn = np.zeros((trial_length, len(temp)))
##o_sampled_count = np.zeros((int(9e5), len(temp)))
#
#
#
#for k in range(len(temp)):
#    for i in range(len(trials)):
#        
#        grid = np.ones((spins, spins))
#        energy_avg, magnet_avg, C_v, susceptibility, abs_magnet, counter_list = MC(grid, trials[i], temp[k])
#        if i >Temp:
#            o_sampled_energies[i - Temp, k] = energy_avg
#            o_sampled_magnets[i - Temp, k] = np.abs(magnet_avg)
#            o_sampled_cv[i - Temp, k] = C_v
#            o_sampled_suscept[i - Temp, k] = susceptibility
#            o_sampled_absmagn[i - Temp, k] = abs_magnet
##            o_sampled_count[:, k] = counter_list
        

        
stop = time.time()
print('CPU: ', stop-start)  

## Expectation values for E and M



plt.figure()
plt.hist([sampled_energies[:,0]],
    bins = 50, density = 10)
plt.title('Energy distribution with T = 1')

plt.figure()
plt.hist([sampled_energies[:,1]],
    bins = 50, density = 10)
plt.title('Energy distribution with T = 2.4')



#plt.figure()
#ax1 = plt.subplot(211)
#ax1.set_title('Expectation values for energy and magnetization as functions of MC-cycles, T = 1', fontsize = 12)
#ax1.semilogx(trials[Temp+1:], sampled_energies[1:,0], 'r', label= 'T = 1.0')
#ax1.semilogx(trials[Temp+1:], sampled_energies[1:,1], 'b', label = 'T = 2.4')
#ax1.legend()
#ax1.set_ylabel(r'$\langle E \rangle$', fontsize = 10)
#ax1.grid()
#ax2 = plt.subplot(212)
#ax2.semilogx(trials[Temp+1:], sampled_absmagn[1:,0], 'r', label = 'T = 1.0')
#ax2.semilogx(trials[Temp+1:], sampled_absmagn[1:,1], 'b', label = 'T = 2.4')
#ax2.legend()
#ax2.set_ylabel(r'$\langle |M|\rangle$', fontsize = 10)
#ax2.set_xlabel('MC-cycles', fontsize = 10)
#ax2.grid()
##plt.savefig('figs/expectations_2020lattice_t1.pdf', bbox_inches = 'tight')
#plt.show()


#plt.figure()
#ax1 = plt.subplot(211)
#ax1.set_title('Expectation values for energy and magnetization as functions of MC-cycles, T = 2.4', fontsize = 12)
#ax1.plot(trials[Temp+1:], sampled_energies[1:,1], 'r--', label= 'ordered')
#ax1.plot(trials[Temp+1:], o_sampled_energies[1:,1], 'b--', label = 'random')
#ax1.legend()
#ax1.set_ylabel(r'$\langle E \rangle$', fontsize = 10)
#ax1.grid()
#ax2 = plt.subplot(212)
#ax2.plot(trials[Temp+1:], sampled_absmagn[1:,1], 'r--', label ='ordered')
#ax2.plot(trials[Temp+1:], o_sampled_absmagn[1:,1], 'b--', label = 'random')
#ax2.legend()
#ax2.set_ylabel(r'$\langle |M|\rangle$', fontsize = 10)
#ax2.set_xlabel('MC-cycles', fontsize = 10)
#ax2.grid()
##plt.savefig('figs/expectations_2020lattice_t2.pdf', bbox_inches = 'tight')
#plt.show()

### Number of accepted states: 
#t = np.linspace(1, trials[0], trials[0])
#plt.figure()
#plt.plot(t, sampled_count[:,0], label = 'Random spin config.')
#plt.plot(t, o_sampled_count[:,0], label = 'Ordered spin config.')
#plt.legend()
#plt.grid("on")
#plt.xlabel('MC-cycles', fontsize = 10)
#plt.ylabel('Number of accepted states', fontsize = 10)
#plt.title('Temperature = 1.0', fontsize = 10)
#plt.savefig('figs/n_accepted_states_t1.pdf', bbox_inches = 'tight')
#plt.figure()
#plt.plot(t, sampled_count[:,1], label = 'Random spin config.')
#plt.plot(t, o_sampled_count[:,1], label = 'Ordered spin config.')
#plt.legend()
#plt.grid("on")
#plt.xlabel('MC-cycles', fontsize = 10)
#plt.ylabel('Number of accepted states', fontsize = 10)
#plt.title('Temperature = 2.4', fontsize = 10)
#plt.savefig('figs/n_accepted_states_t2_4.pdf', bbox_inches = 'tight')