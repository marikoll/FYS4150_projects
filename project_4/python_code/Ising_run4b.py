#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Expectation values for energy, magnetization, C_V and susceptibility for a 
2x2 lattice with various mc_cycles and temperature = 1.0. 
"""

import numpy as np
import time
from numba import prange

from Ising_model import MC


spins       = 2
mc_cycles   = [int(1e2), int(1e3), int(1e4), int(1e5), int(1e6), int(1e7)]
temp        = 1.0



start = time.time()


for i in prange(len(mc_cycles)):
    grid = np.ones((spins, spins))
    energy_avg, magnet_avg, C_v, susceptibility, abs_magnet, c= MC(grid, mc_cycles[i], temp)
    print('MC-cycles: {}\n<E>: {:.4f}\n<M>: {:.4f}\nC_v: {:.4f}\nchi: {:.4f}\n\n'.\
          format(mc_cycles[i], energy_avg[-1], abs_magnet[-1], C_v, susceptibility))


        
stop = time.time()
print('CPU: ', stop-start)  

