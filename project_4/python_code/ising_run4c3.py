#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Accepted states - task 4c)
"""
import numpy as np
import time
import matplotlib.pyplot as plt

from Ising_model import MC

spins       = 20
trials      = int(1e6)#np.arange(100, 1e6, 5000, np.dtype(np.int64))

temp = [1.0, 2.4]

sampled_count = np.zeros((int(1e6), len(temp)))

start = time.time()

for k in range(len(temp)):
    grid = np.random.choice([-1,1],size=(spins, spins))
    energy_avg, magnet_avg, C_v, susceptibility, abs_magnet, counter_list = MC(grid, trials, temp[k])
    sampled_count[:, k] = counter_list

o_sampled_count = np.zeros((int(1e6), len(temp)))
for k in range(len(temp)):
    grid = np.ones((spins, spins))
    energy_avg, magnet_avg, C_v, susceptibility, abs_magnet, counter_list = MC(grid, trials, temp[k])
    o_sampled_count[:, k] = counter_list

stop = time.time()
print('CPU: ', stop-start)  
       
## Number of accepted states: 
t = np.linspace(1, trials, trials)
plt.figure()
plt.plot(t, sampled_count[:,0], label = 'Random spin config.')
plt.plot(t, o_sampled_count[:,0], label = 'Ordered spin config.')
plt.legend()
plt.grid("on")
plt.xlabel('MC-cycles', fontsize = 10)
plt.ylabel('Number of accepted states', fontsize = 10)
plt.title('Temperature = 1.0', fontsize = 10)
plt.savefig('figs/n_accepted_states_t1.pdf', bbox_inches = 'tight')
plt.figure()
plt.plot(t, sampled_count[:,1], label = 'Random spin config.')
plt.plot(t, o_sampled_count[:,1], label = 'Ordered spin config.')
plt.legend()
plt.grid("on")
plt.xlabel('MC-cycles', fontsize = 10)
plt.ylabel('Number of accepted states', fontsize = 10)
plt.title('Temperature = 2.4', fontsize = 10)
plt.savefig('figs/n_accepted_states_t2_4.pdf', bbox_inches = 'tight')