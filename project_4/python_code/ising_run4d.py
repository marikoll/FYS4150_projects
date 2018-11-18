"""
Plots the energy distribution for a 20x20 lattice in 1e6 MC-cycles for 
temperature T = 1.0 and T = 2.4 and belonging percentage-cutoffs. 
"""


import numpy as np
import time
import matplotlib.pyplot as plt
import scipy.stats as st
import matplotlib.mlab as mlab

from Ising_model import MC_cutoff


spins       = 20
trials      = int(1e6)

temp = [1.0, 2.4]

p = [0.01, 0.1]



sampled_energies = np.zeros((trials, len(temp)))
sampled_variance = np.zeros(len(temp))

start = time.time()


for k, T in enumerate(temp):
    grid = np.random.choice([-1,1],size=(spins, spins))
    En,energy_var, energy_avg, magnet_avg, C_v, susceptibility, abs_magnet = MC_cutoff(grid, trials, T, P = p[k])
    sampled_energies[:, k] = En
    sampled_variance[k] = energy_var/spins**2

        
stop = time.time()
print('CPU: ', stop-start)  

sampled_energies_24 = sampled_energies[(int(trials*p[1])):,1]
sampled_energies_1 = sampled_energies[(int(trials*p[0])):,0]


plt.figure()
(mu, sigma2) = st.norm.fit([sampled_energies_24])
n, bins, patches = plt.hist([sampled_energies_24], 30, normed = 1, label = r'$\sigma^2 = {:.2f}$'.format(sampled_variance[1]))
y = mlab.normpdf( bins, mu, sigma2)
l = plt.plot(bins, y, 'r--', linewidth=2, label = r'$\sigma^2 = {:.2f}$'.format(sigma2))
plt.legend()
plt.xlabel('Energy per spin', fontsize = 10)
plt.ylabel('Normalized number of occurences', fontsize = 10)
plt.title('Energy distribution with T = 2.4\nMC-cutoff: {}%'.format(p[1]*100), fontsize = 15)
plt.savefig('figs/energy_dist_T24.pdf', bbox_inches = 'tight')
plt.show()

plt.figure()
(mu, sigma1) = st.norm.fit([sampled_energies_1])
n, bins, patches = plt.hist([sampled_energies_1], 30, normed = 1, label = r'$\sigma^2 = {:.2e}$'.format(sampled_variance[0]))
y = mlab.normpdf( bins, mu, sigma1)
l = plt.plot(bins, y, 'r--', linewidth=2, label = r'$\sigma^2 = {:.2e}$'.format(sigma1))
plt.legend()
plt.xlabel('Energy per spin', fontsize = 10)
plt.ylabel('Normalized number of occurences', fontsize = 10)
plt.title('Energy distribution with T = 1.0\nMC-cutoff: {}%'.format(p[0]*100), fontsize = 15)
plt.savefig('figs/energy_dist_T1.pdf', bbox_inches = 'tight')
plt.show()