#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Expectation values as funct. of MC-cycles
"""

import numpy as np
import time
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes
from mpl_toolkits.axes_grid1.inset_locator import mark_inset

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


T = np.linspace(0, trials, trials/10)



plt.figure(1)
fig, ax  = plt.subplots()
ax.plot(T, o_sampled_energies[::10,0], 'r', label= 'ordered configuration')
ax.plot(T, sampled_energies[::10,0], 'b', label= 'random configuration')
ax.legend(loc = 9)
ax.set_title('Expectation values for energy, T = 1.0', fontsize = 15)
ax.set_xlabel('MC-cycles', fontsize = 10)
ax.set_ylabel(r'$\langle E\rangle$', fontsize = 10)
axins = zoomed_inset_axes(ax, 19, loc=1)
axins.plot(T, o_sampled_energies[::10,0], 'r')
axins.plot(T, sampled_energies[::10,0], 'b')
x1, x2, y1, y2 = 0, 10000, -2.01, -1.95 # specify the limits
axins.set_xlim(x1, x2) # apply the x-limits
axins.set_ylim(y1, y2)
mark_inset(ax, axins, loc1=2, loc2=4, fc="none", ec="0.5")
plt.savefig('figs/expectations_2020lattice_temp1_energy.pdf', bbox_inches = 'tight')
plt.draw()
plt.show()

plt.figure(2)
fig, ax  = plt.subplots()
ax.plot(T, o_sampled_absmagn[::10,0], 'r',label= 'ordered configuration')
ax.plot(T, sampled_absmagn[::10,0], 'b', label= 'random configuration')
ax.legend(loc = 8)
ax.set_title('Expectation values for magnetization, T = 1.0', fontsize = 15)
ax.set_xlabel('MC-cycles', fontsize = 10)
ax.set_ylabel(r'$\langle |M|\rangle$', fontsize = 10)
axins = zoomed_inset_axes(ax, 11, loc=7)
axins.plot(T, o_sampled_absmagn[::10,0], 'r')
axins.plot(T, sampled_absmagn[::10,0], 'b')
x1, x2, y1, y2 = -10, 15000, 0.94, 1.01 # specify the limits
axins.set_xlim(x1, x2) # apply the x-limits
axins.set_ylim(y1, y2)
mark_inset(ax, axins, loc1=2, loc2=4, fc="none", ec="0.5")
plt.savefig('figs/expectations_2020lattice_temp1_magnet.pdf', bbox_inches = 'tight')
plt.draw()
plt.show()

plt.figure(3)
fig, ax  = plt.subplots()
ax.plot(T, o_sampled_energies[::10,1], 'r', label= 'ordered configuration')
ax.plot(T, sampled_energies[::10,1], 'b', label= 'random configuration')
ax.legend(loc = 9)
ax.set_title('Expectation values for energy, T = 2.4', fontsize = 15)
ax.set_xlabel('MC-cycles', fontsize = 10)
ax.set_ylabel(r'$\langle E\rangle$', fontsize = 10)
axins = zoomed_inset_axes(ax, 8, loc=7)
axins.plot(T, o_sampled_energies[::10,1], 'r')
axins.plot(T, sampled_energies[::10,1], 'b')
x1, x2, y1, y2 = 0, 55000, -1.3, -1.2 # specify the limits
axins.set_xlim(x1, x2) # apply the x-limits
axins.set_ylim(y1, y2)
mark_inset(ax, axins, loc1=2, loc2=4, fc="none", ec="0.5")
plt.savefig('figs/expectations_2020lattice_temp24_energy.pdf', bbox_inches = 'tight')
plt.draw()
plt.show()

plt.figure(4)
fig, ax  = plt.subplots()
ax.plot(T, o_sampled_absmagn[::10,1], 'r',label= 'ordered configuration')
ax.plot(T, sampled_absmagn[::10,1], 'b', label= 'random configuration')
ax.legend(loc = 9)
ax.set_title('Expectation values for magnetization, T = 2.4', fontsize = 15)
ax.set_xlabel('MC-cycles', fontsize = 10)
ax.set_ylabel(r'$\langle |M|\rangle$', fontsize = 10)
axins = zoomed_inset_axes(ax, 4, loc=7)
axins.plot(T, o_sampled_absmagn[::10,1], 'r')
axins.plot(T, sampled_absmagn[::10,1], 'b')
x1, x2, y1, y2 = -10, 170000, 0.4, 0.56 # specify the limits
axins.set_xlim(x1, x2) # apply the x-limits
axins.set_ylim(y1, y2)
mark_inset(ax, axins, loc1=2, loc2=4, fc="none", ec="0.5")
plt.savefig('figs/expectations_2020lattice_temp24_magnet.pdf', bbox_inches = 'tight')
plt.draw()
plt.show()

#plt.figure()
#ax1 = plt.subplot(211)
#ax1.set_title('Expectation values for energy and magnetization as functions of MC-cycles, T = 1.0', fontsize = 12)
#ax1.plot(T, o_sampled_energies[::10,0], 'r', label= 'ordered configuration')
#ax1.plot(T, sampled_energies[::10,0], 'b', label = 'random configuration')
#ax1.legend()
#ax1.set_ylabel(r'$\langle E \rangle$', fontsize = 10)
#ax1.grid()
#ax2 = plt.subplot(212)
#ax2.plot(T, o_sampled_absmagn[::10,0], 'r', label = 'ordered configuration')
#ax2.plot(T, sampled_absmagn[::10,0], 'b', label = 'random configuration')
#ax2.legend()
#ax2.set_ylabel(r'$\langle |M|\rangle$', fontsize = 10)
#ax2.set_xlabel('MC-cycles', fontsize = 10)
#ax2.grid()
#plt.savefig('figs/expectations_2020lattice_temp1.pdf', bbox_inches = 'tight')
#plt.show()
#
#
#
#plt.figure()
#ax1 = plt.subplot(211)
#ax1.set_title('Expectation values for energy and magnetization as functions of MC-cycles, T = 2.4', fontsize = 12)
#ax1.plot(T, o_sampled_energies[::10,1], 'r', label= 'ordered configuration')
#ax1.plot(T, sampled_energies[::10,1], 'b', label = 'random configuration')
#ax1.legend()
#ax1.set_ylabel(r'$\langle E \rangle$', fontsize = 10)
#ax1.grid()
#ax2 = plt.subplot(212)
#ax2.plot(T, o_sampled_absmagn[::10,1], 'r', label = 'ordered configuration')
#ax2.plot(T, sampled_absmagn[::10,1], 'b', label = 'random configuration')
#ax2.legend()
#ax2.set_ylabel(r'$\langle |M|\rangle$', fontsize = 10)
#ax2.set_xlabel('MC-cycles', fontsize = 10)
#ax2.grid()
#plt.savefig('figs/expectations_2020lattice_temp24.pdf', bbox_inches = 'tight')
#plt.show()
