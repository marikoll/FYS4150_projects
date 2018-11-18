#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Checks and plots number of accepted states for a 20x20 lattice for temperature 
T = 1.0 and T = 2.4 with ordered (i.e all spins up) and unordered configuration
"""
import numpy as np
import time
import matplotlib.pyplot as plt

from Ising_model import MC

spins       = 20
trials      = int(1e4)

temp = [1.0, 2.4]

sampled_count = np.zeros((int(1e4), len(temp)))

start = time.time()

for k in range(len(temp)):
    grid = np.random.choice([-1,1],size=(spins, spins))
    energy_avg, magnet_avg, C_v, susceptibility, abs_magnet, counter_list = MC(grid, trials, temp[k])
    sampled_count[:, k] = counter_list

o_sampled_count = np.zeros((int(1e4), len(temp)))
for k in range(len(temp)):
    grid = np.ones((spins, spins))
    energy_avg, magnet_avg, C_v, susceptibility, abs_magnet, counter_list = MC(grid, trials, temp[k])
    o_sampled_count[:, k] = counter_list

stop = time.time()
print('CPU: ', stop-start)  
       
## Number of accepted states: 
t = np.linspace(1, trials, trials)
plt.figure()
plt.plot(t, sampled_count[:,0], 'r', label = 'Random spin config.')
plt.plot(t, o_sampled_count[:,0], 'b', label = 'Ordered spin config.')
plt.legend()
plt.grid("on")
plt.xlabel('MC-cycles', fontsize = 10)
plt.ylabel('Number of accepted states', fontsize = 10)
plt.title('Temperature = 1.0', fontsize = 15)
plt.savefig('figs/n_accepted_states_t1.pdf', bbox_inches = 'tight')
plt.figure()
plt.plot(t, sampled_count[:,1], 'r', label = 'Random spin config.')
plt.plot(t, o_sampled_count[:,1], 'b', label = 'Ordered spin config.')
plt.legend()
plt.grid("on")
plt.xlabel('MC-cycles', fontsize = 10)
plt.ylabel('Number of accepted states', fontsize = 10)
plt.title('Temperature = 2.4', fontsize = 15)
plt.savefig('figs/n_accepted_states_t2_4.pdf', bbox_inches = 'tight')