#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 18 16:34:04 2018

@author: maritkollstuen
"""

import numpy as np
import solar_system_class as ss
import Earth_Jupiter_Sun as ejs
import Earth_Sun as es
import Mercury_Sun as ms
import matplotlib.pyplot as plt



massM   = 3.3E23
posM    =  np.array([-4.466974942983433E-02,-4.551297184815911E-01,-3.377443034243644E-02])
velM    =  np.array([2.235388667387174E-02,-1.255928387575858E-03,-2.154047309352793E-03])*365

posM    = np.array([0.3075,0,0])
velM    = np.array([0,12.44,0])

massV   = 4.9E24
posV    =  np.array([6.771183423096696E-01,2.635570892119448E-01,-3.564015770708658E-02])
velV    =  np.array([-7.233924240439758E-03,1.883053303306194E-02,6.755080544414482E-04])*365

massE   = 6E24
posE    = np.array([9.004267194046488E-01,4.329702891250327E-01,-9.309259935529284E-05])
velE    = np.array([-7.644048786784979E-03,1.548720517812966E-02,9.799447004349826E-08])*365

massMa  = 6.6E23
posMa   = np.array([1.386064640638050E+00,-7.195522876836125E-02,-3.574720875419098E-02])
velMa   = np.array([1.325375632613797E-03,1.516910673406354E-02,2.852845787253853E-04])*365

masstimes = 1
massJ   = 1.9E27 #*masstimes
posJ    = np.array([-2.628114780107304E+00,-4.675856649493857E+00,7.818020939743552E-02])
velJ    = np.array([6.487854045762937E-03,-3.337752963125740E-03,-1.312915365947336E-04])*365


massSa  = 5.5E26
posSa   = np.array([1.575176167525727E+00,-9.930952889722587E+00,1.099718612784312E-01])
velSa   = np.array([5.202593393091781E-03,8.555506892745189E-04,-2.216362261569332E-04])*365

massU   = 8.8E25
posU    = np.array([1.716586975644875E+01,1.001373868288497E+01,-1.851949198554858E-01])
velU    = np.array([-2.010719471096281E-03,3.213990583409865E-03,3.805302445255308E-05])*365

massN   = 1.03E26
posN    = np.array([2.892424054587586E+01,-7.707281387885470E+00,-5.078715497631403E-01])
velN    = np.array([7.870585771863339E-04,3.052076356482153E-03,-8.060978991971341E-05])*365

massP   = 1.31E22
posP    = np.array([1.165895866372055E+01,-3.157299883744404E+01,6.054513403127080E-03])
velP    = np.array([3.023518143470085E-03,4.354001925477489E-04,-9.190982697153275E-04])*365

 
massS   = 2E30

#velS    = massE*velE/massS
#
#Earth   = ejs.planets(velE, posE, massE)
#Sun     = ejs.planets(velS, posS, massS)          
#
#system = ejs.system(Earth, Sun)
#solver = ejs.solver(system, 8000, 1/3650)
#solver.velocity_verlet()
#
#for i in range(len(system)):
#    plt.plot(system[i].pos_vec[:,0],system[i].pos_vec[:,1])
#plt.show()
    

## 3c)

#dt = 10.0**np.arange(-1, -7, -1)# [1/3.65, 1/36.5, 1/365, 1/3650, 1/36500]
#planet = ['Sun', 'Earth']
#energy = []
#angular_momentum = []
## plt.figure()
#for i in dt:
#    #for method in ['.velocity_verlet()', '.euler_fwd()']:
#    posS    = np.array([0,0,0])
#    velS    = massE*velE/massS
#    Earth   = es.planets(velE, posE, massE)
#    Sun     = es.planets(velS, posS, massS)
#
#    systemES = es.system(Sun, Earth)
#
#    solveES  = es.solver(systemES, int(1/i), i)
#    solveES.velocity_verlet()
#    # deltas.append(np.linalg.norm(systemES[1].pos_vec[0,:]- systemES[1].pos_vec[-1,:]))
#    kin_energy = systemES[1].kinetic_energy(systemES[1].vel_vec)
#    pot_energy = systemES[1].potential_energy(systemES[1].pos_vec)
#    energy.append((kin_energy[-1] + pot_energy[-1])-(kin_energy[0] + pot_energy[0]))
##    for l in range(int(1/i)):
##        r = np.linalg.norm(systemES[1].pos_vec, axis = 1)
#    l1 = systemES[1].angular_momentum(systemES[1].vel_vec[0,:], systemES[1].pos_vec[0,:])
#    l2 = systemES[1].angular_momentum(systemES[1].vel_vec[-1,:], systemES[1].pos_vec[-1,:])
#    angular_momentum.append(l2 - l1)
#plt.semilogx(dt, np.linalg.norm(angular_momentum, axis = 1))
#plt.ticklabel_format(axis = 'y', style = 'sci', scilimits=(0,0))
#plt.grid()
#plt.xlabel('dt', fontsize = 10)
#plt.ylabel('L2 - L1', fontsize = 10)
#plt.title('L2 - L1 with Euler method', fontsize = 15)
#plt.savefig('Euler_angmom.pdf')
#    plt.semilogx(i, systemES[1].pos_vec[0,0]- systemES[1].pos_vec[-1,0])
#    plt.xlabel('dt')
#    plt.ylabel('position')
#    for k in range(len(systemES)):
#        
#        solveES  = es.solver(systemES, 8000, i)
#        solveES.velocity_verlet()
#        plt.plot(systemES[k].pos_vec[:,0],systemES[k].pos_vec[:,1], label = planet[k])
#        plt.title('Method Verlet, dt: {}'.format(i))
#        plt.axis('equal')
#    plt.legend()


#    pos_e = solveES.euler_fwd()
#    ax1 =plt.subplot(1, 2, 1)
#    ax1.plot(pos_v[:,0],pos_v[:,1])# , label = planets[i])
###    ax1.set_xlim([-20, 400])
###    ax1.set_ylim([-180, 50])
#    ax1.set_title('Velocity Verlet', fontsize  = 10)
#    ax1.set_ylabel('y [AU]', fontsize  = 10)
#    ax1.set_xlabel('x [AU]', fontsize  = 10)
##    ax1.legend()
#    ax2 = plt.subplot(1, 2, 2)
#    ax2.plot(pos_e[:,0],pos_e[:,1])#, label = planets[i])
#    ax2.set_title('Euler', fontsize  = 10)
#    ax2.set_xlabel('x [AU]', fontsize  = 10)
##    ax2.legend()




## 3f - TAR 100000 ÅR!            
#velS    = -(massSa*velSa + massM*velM + massE*velE + massJ*velJ + massN*velN \
#            + massP*velP + massU*velU + massV*velV + massMa*velMa)/massS
#    
#Mercury = ss.planets(velM, posM, massM)
#Earth   = ss.planets(velE, posE, massE)
#Jupiter = ss.planets(velJ, posJ, massJ)
#Saturn  = ss.planets(velSa, posSa, massSa)
#Sun     = ss.planets(velS, posS, massS)
#Venus   = ss.planets(velV, posV, massV)
#Mars    = ss.planets(velMa, posMa, massMa)
#Uranus  = ss.planets(velU, posU, massU)
#Neptun  = ss.planets(velN, posN, massN)
#Pluto   = ss.planets(velP, posP, massP)
#
#system  = ss.system(Sun, Mercury, Venus, Earth, Mars, Jupiter, Saturn, Uranus, Neptun, Pluto)
#solve   = ss.solver(system, 1000000,1/3650)
#solve.velocity_verlet()
#
#
#planets = ['Sun', 'Mercury', 'Venus', 'Earth', 'Mars', 'Jupiter', 'Saturn', 'Uranus', 'Neptun', 'Pluto']
#plt.figure()
#for i in range(len(system)):
#    plt.plot(system[i].pos_vec[:,0],system[i].pos_vec[:,1], label = planets[i])
#plt.legend()
#plt.title('Solar System', fontsize = 15)
#plt.xlabel('x [AU]', fontsize = 10)
#plt.ylabel('y [AU]', fontsize = 10)
#plt.savefig('solar_system.pdf')

            
            
## Task 3e): 

#velS1 = np.array([0,0,0])
#
#Earth1   = ejs.planets(velE, posE, massE)
#Jupiter1 = ejs.planets(velJ, posJ, massJ)
#Sun1     = ejs.planets(velS1, posS, massS)
#N = 100000
#dt = 1/3650
#system1 = ejs.system(Sun1, Earth1, Jupiter1)
#solve1  = ejs.solver(system1, N, dt)
#
#solve1.velocity_verlet()
#
#
#velS      = -(massSa*velSa + massM*velM + massE*velE + massJ*velJ)/massS
#
##Mercury = ss.planets(velM, posM, massM)
#Earth   = ss.planets(velE, posE, massE)
#Jupiter = ss.planets(velJ, posJ, massJ)
##Saturn  = ss.planets(velSa, posSa, massSa)
#Sun     = ss.planets(velS, posS, massS)
#
#
#system = ss.system(Sun, Earth, Jupiter)#, Mercury, Saturn)
#solve = ss.solver(system, 100000,1/3650)
#solve.velocity_verlet()
#
#planets = ['Sun', 'Earth', 'Jupiter']
#
#fig = plt.figure()
#fig.suptitle('Mass of Jupiter: {:.1e} kg'.format(massJ), fontsize = 15)
#for i in range(len(system1)):
#    ax1 =plt.subplot(1, 2, 1)
#    ax1.plot(system1[i].pos_vec[:,0],system1[i].pos_vec[:,1], label = planets[i])
##    ax1.set_xlim([-20, 400])
##    ax1.set_ylim([-180, 50])
#    ax1.set_title('Center of mass position at Sun', fontsize  = 10)
#    ax1.set_ylabel('y [AU]', fontsize  = 10)
#    ax1.set_xlabel('x [AU]', fontsize  = 10)
#    ax1.legend()
#    ax2 = plt.subplot(1, 2, 2)
#    ax2.plot(system[i].pos_vec[:,0],system[i].pos_vec[:,1], label = planets[i])
#    ax2.set_title('Center of mass position at origin', fontsize  = 10)
#    ax2.set_xlabel('x [AU]', fontsize  = 10)
#    ax2.legend()
#
#plt.savefig('task3e_massJ{}.pdf'.format(masstimes))
#plt.show()


## Task 3f)


dt = 1e-7
N = int(3e7)
p_vel = np.array(ms.velocity_verlet(posM, velM, N, dt))
# np.save("p_vel.txt", p_vel)

#planets = ['Sun', 'Mercury']
#
#fig = plt.figure()
##fig.suptitle('Mass of Jupiter: {:.1e} kg'.format(massJ), fontsize = 15)
#for i in range(len(system)):
##    ax1 =plt.subplot(1, 2, 1)
#    plt.plot(system[i].pos_vec[:,0],system[i].pos_vec[:,1], label = planets[i])
##    ax1.set_xlim([-20, 400])
##    ax1.set_ylim([-180, 50])
#    plt.title('Sun - Mercury', fontsize  = 10)
#    plt.ylabel('y [AU]', fontsize  = 10)
#    plt.xlabel('x [AU]', fontsize  = 10)
#    plt.legend()




