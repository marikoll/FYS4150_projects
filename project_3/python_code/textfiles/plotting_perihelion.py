""" file for plotting npy files with Mercury-Sun positions at perihelion"""


import numpy as np
import matplotlib.pyplot as plt


p_100 = np.load('precession.npy')
p_50 = np.load('perihelion_50yrs.npy')
r_0 = np.load('precession_not_relativistic.npy')
r_100 = np.load('precession_not_relativistic_100yrs.npy')

plt.figure()
plt.plot(np.arctan(p_100[:, 1]/p_100[:, 0])*180*3600/np.pi, 'b', label = 'corrected Newton, N = 1e8, dt = 1e-6')
plt.plot(np.arctan(p_50[:, 1]/p_50[:, 0])*180*3600/np.pi, 'g', label = 'corrected Newton, N = 1e9, dt = 1e-7')
plt.plot(np.arctan(r_100[:, 1]/r_100[:, 0])*180*3600/np.pi, 'r', label = 'pure Newtonian, N = 1e8, dt = 1e-6')
plt.plot(442, 43, 'ro', linewidth = 4, label = 'observed value')
plt.legend(loc = 'upper left', fontsize = 13)
plt.grid()
plt.xlabel('No. perihelions', fontsize = 12)
plt.ylabel('arc seconds', fontsize = 12)
plt.title('Perihelion precession of Mercury', fontsize = 15)
plt.savefig('perhelion_precession.pdf')
plt.show()
