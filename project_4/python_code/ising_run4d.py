mport matplotlib
import numpy as np
import time
import matplotlib.pyplot as plt

from Ising_model import MC


spins       = 20
trials      = np.arange(100, 1e6, 1000, np.dtype(np.int64))#[int(1e2), int(1e3), int(1e4), int(1e5), int(1e6), int(1e7)]

temp = [1.0, 2.4]

sampled_energies = np.zeros((len(trials), len(temp)))
sampled_magnets = np.zeros((len(trials), len(temp)))
sampled_cv = np.zeros((len(trials), len(temp)))
sampled_suscept = np.zeros((len(trials), len(temp)))
sampled_absmagn = np.zeros((len(trials), len(temp)))

start = time.time()


for k in range(len(temp)):
    for i in range(len(trials)):
        grid = np.ones((spins, spins))
        energy_avg, magnet_avg, C_v, susceptibility, abs_magnet = MC(grid, trials[i], temp[k])
        sampled_energies[i, k] = energy_avg
        sampled_magnets[i, k] = np.abs(magnet_avg)
        sampled_cv[i, k] = C_v
        sampled_suscept[i, k] = susceptibility
        sampled_absmagn[i, k] = abs_magnet

o_sampled_energies = np.zeros((len(trials), len(temp)))
o_sampled_magnets = np.zeros((len(trials), len(temp)))
o_sampled_cv = np.zeros((len(trials), len(temp)))
o_sampled_suscept = np.zeros((len(trials), len(temp)))
o_sampled_absmagn = np.zeros((len(trials), len(temp)))



for k in range(len(temp)):
    for i in range(len(trials)):
        grid = np.ones((spins, spins))
        energy_avg, magnet_avg, C_v, susceptibility, abs_magnet = MC(grid, trials[i], temp[k], ordered = True)
        o_sampled_energies[i, k] = energy_avg
        o_sampled_magnets[i, k] = np.abs(magnet_avg)
        o_sampled_cv[i, k] = C_v
        o_sampled_suscept[i, k] = susceptibility
        o_sampled_absmagn[i, k] = abs_magnet


stop = time.time()
print('CPU: ', stop-start)



fig, axs = plt.subplots(1,2, tight_layout=True, figsize= (10,8))
axs[0].hist([sampled_energies[:,0] , o_sampled_energies[:,0]],
    bins = 30, density = 5)
axs[0].set_xlim(-2, -1.995, 3)
axs[0].set_ylabel('Probability')
axs[0].set_title('Energy distribution with T = 1')
axs[0].set_xlabel('E')
axs[0].grid()

axs[1].hist([sampled_magnets[:,0], o_sampled_magnets[:,0]],
    bins = 30, density = 1.5,)
axs[1].set_xlim(0.999, 1.00)
axs[1].legend(['random', 'ordered'], loc = 1)
axs[1].set_title('Magnetization distribution with T = 1')
axs[1].set_xlabel('M')
axs[1].grid()

plt.show()



plt.figure()
ax1 = plt.subplot(211)
ax1.set_title('Expectation values for energy and magnetization as functions of MC-cycles, T = 2.4', fontsize = 12)
ax1.plot(trials, sampled_energies[:,1], 'r-*', label= 'ordered')
ax1.plot(trials, o_sampled_energies[:,1], 'b--', label = 'random')
ax1.legend()
ax1.set_ylabel(r'$\langle E \rangle$', fontsize = 10)
ax1.grid()
ax2 = plt.subplot(212)
ax2.plot(trials, sampled_magnets[:,1], 'r-*', label ='ordered')
ax2.plot(trials, o_sampled_magnets[:,1], 'b--', label = 'random')
ax2.legend()
ax2.set_ylabel(r'$\langle M\rangle$', fontsize = 10)
ax2.set_xlabel('MC-cycles', fontsize = 10)
ax2.grid()
plt.savefig('expectations_2020lattice_t2.pdf', bbox_inches = 'tight')
plt.show()
