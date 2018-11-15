
import numpy as np
import time
import matplotlib.pyplot as plt
import scipy.stats as st

from Ising_model import MC_cutoff


spins       = 20
trials      = int(1e6)# np.logspace(2, 7, 500, base = 10.0, dtype = np.dtype(np.int64)) #, np.dtype(np.int64))#[int(1e2), int(1e3), int(1e4), int(1e5), int(1e6), int(1e7)]

temp = [1.0, 2.4]

p = [0.0015, 0.0025]



sampled_energies = np.zeros((trials, len(temp)))


start = time.time()


for k, T in enumerate(temp):
    grid = np.random.choice([-1,1],size=(spins, spins))
    energy_avg, magnet_avg, C_v, susceptibility, abs_magnet,\
        sampling_starts_from = MC_cutoff(grid, trials, T, P = p[k])
    sampled_energies[:, k] = energy_avg

        
stop = time.time()
print('CPU: ', stop-start)  



plt.figure()
plt.hist([sampled_energies[:,0]],
    bins = 70, density = 10)
plt.title('Energy distribution with T = 1')

plt.figure()
plt.hist([sampled_energies[:,1]],
    bins = 70, density = 10)
plt.title('Energy distribution with T = 2.4')


