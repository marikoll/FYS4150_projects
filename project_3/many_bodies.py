#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 18 16:34:04 2018

@author: maritkollstuen
"""

import numpy as np
import solar_system_class as ss
import matplotlib.pyplot as plt


massM= 3.3E23
posM =  np.array([-2.138801978256535E-01,-4.028538544608614E-01,-1.397398026440086E-02])
velM =  np.array([1.927259979627735E-02,-1.164174161133437E-02,-2.720041450765680E-03])*365


massE = 6E24
posE  = np.array([9.528047055398201E-01,3.053612869840809E-01,-9.272902073041313E-05])
velE    = np.array([-5.428888690270241E-03,1.636353485946535E-02,-4.491683144318728E-07])*365


massJ = 1.9E27
posJ   = np.array([-2.679859418467178E+00,-4.648870533678862E+00,7.922561878424454E-02])
velJ   = np.array([6.448353939132267E-03,-3.409518873130117E-03,-1.300875035507015E-04])*365


massSa = 5.5E26
posSa    = np.array([1.533542535074663E+00,-9.937711138425550E+00,1.117470354327505E-01])
velSa    = np.array([5.206732001657008E-03,8.336898877431193E-04,-2.220848255804162E-04])*365


massS     = 2E30
posS      = np.array([0,0,0])
velS      = -(massSa*velSa + massM*velM + massE*velE + massJ*velJ)/massS


Mercury = ss.planets(velM, posM, massM)
Earth = ss.planets(velE, posE, massE)
Jupiter = ss.planets(velJ, posJ, massJ)
Saturn = ss.planets(velSa, posSa, massSa)
Sun = ss.planets(velS, posS, massS)


system = ss.system(Sun, Earth, Jupiter, Mercury, Saturn)
solve = ss.solver(system, 80000,1/3650)
solve.velocity_verlet()

for i in range(len(system)):
    plt.plot(system[i].pos_vec[:,0],system[i].pos_vec[:,1])
plt.show()