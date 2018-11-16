#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Nov 11 14:15:51 2018

@author: maritkollstuen
"""

import numpy as np
import time

from Ising_model import MC


spins       = 2
trials      = [int(1e2), int(1e3), int(1e4), int(1e5), int(1e6), int(1e7)]
temp = 1.0



start = time.time()


for i in range(len(trials)):
    grid = np.ones((spins, spins))
    energy_avg, magnet_avg, C_v, susceptibility, abs_magnet, c= MC(grid, trials[i], temp)
    print('MC-cycles: {}\n<E>: {:.4f}\n<M>: {:.4f}\nC_v: {:.4f}\nchi: {:.4f}\n\n'.\
          format(trials[i], energy_avg[-1], abs_magnet[-1], C_v, susceptibility))


        
stop = time.time()
print('CPU: ', stop-start)  

