
import numpy as np
import time
import matplotlib.pyplot as plt
import scipy.stats as st

from Ising_model import MC

spins       = 20
trials      = np.arange(100, 1e6, 1000, np.dtype(np.int64))

temp = [1.0, 2.4]

#p = 0.1
#trial_length = int(trials*(1-p))
#Temp = trials - trial_length



sampled_energies = np.zeros((len(trials), len(temp)))
sampled_magnets = np.zeros((len(trials), len(temp)))
sampled_cv = np.zeros((len(trials), len(temp)))
sampled_suscept = np.zeros((len(trials), len(temp)))
sampled_absmagn = np.zeros((len(trials), len(temp)))

start = time.time()


for k in range(len(temp)):
    for i in range(len(trials)):
        grid = np.random.choice([-1,1],size=(spins, spins))
        energy_avg, magnet_avg, C_v, susceptibility, abs_magnet, counter_list = MC(grid, trials[i], temp[k])
        sampled_energies[i, k] = energy_avg
        sampled_magnets[i, k] = np.abs(magnet_avg)
        sampled_cv[i, k] = C_v
        sampled_suscept[i, k] = susceptibility
        sampled_absmagn[i, k] = abs_magnet


        
stop = time.time()
print('CPU: ', stop-start)  



plt.figure()
plt.hist([sampled_energies[:,0]],
    bins = 50, density = 10)
plt.title('Energy distribution with T = 1')

plt.figure()
plt.hist([sampled_energies[:,1]],
    bins = 50, density = 10)
plt.title('Energy distribution with T = 2.4')


