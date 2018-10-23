#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Script that reproduces the figures used in project 3 in FYS4150

The script imports the modules needed

Use script by uncommenting the sections separated by rows of #'s

To obtain information of modules, look up the specific module for explanation
"""

import numpy as np
import matplotlib.pyplot as plt


########## PLANETS ################################################
## POSITION AND VELOCITY TAKEN FROM NASA ON THE 19/10-18

## MERCURY
massM   = 3.3E23
posM    =  np.array([-4.466974942983433E-02,-4.551297184815911E-01,-3.377443034243644E-02])
velM    =  np.array([2.235388667387174E-02,-1.255928387575858E-03,-2.154047309352793E-03])*365

## VENUS
massV   = 4.9E24
posV    =  np.array([6.771183423096696E-01,2.635570892119448E-01,-3.564015770708658E-02])
velV    =  np.array([-7.233924240439758E-03,1.883053303306194E-02,6.755080544414482E-04])*365

## EARTH
massE   = 6E24
posE    = np.array([9.004267194046488E-01,4.329702891250327E-01,-9.309259935529284E-05])
velE    = np.array([-7.644048786784979E-03,1.548720517812966E-02,9.799447004349826E-08])*365

## MARS
massMa  = 6.6E23
posMa   = np.array([1.386064640638050E+00,-7.195522876836125E-02,-3.574720875419098E-02])
velMa   = np.array([1.325375632613797E-03,1.516910673406354E-02,2.852845787253853E-04])*365

## JUPITER
massJ   = 1.9E27
posJ    = np.array([-2.628114780107304E+00,-4.675856649493857E+00,7.818020939743552E-02])
velJ    = np.array([6.487854045762937E-03,-3.337752963125740E-03,-1.312915365947336E-04])*365

## SATURN
massSa  = 5.5E26
posSa   = np.array([1.575176167525727E+00,-9.930952889722587E+00,1.099718612784312E-01])
velSa   = np.array([5.202593393091781E-03,8.555506892745189E-04,-2.216362261569332E-04])*365

## URANUS
massU   = 8.8E25
posU    = np.array([1.716586975644875E+01,1.001373868288497E+01,-1.851949198554858E-01])
velU    = np.array([-2.010719471096281E-03,3.213990583409865E-03,3.805302445255308E-05])*365

## NETPUNE
massN   = 1.03E26
posN    = np.array([2.892424054587586E+01,-7.707281387885470E+00,-5.078715497631403E-01])
velN    = np.array([7.870585771863339E-04,3.052076356482153E-03,-8.060978991971341E-05])*365

## PLUTO
massP   = 1.31E22
posP    = np.array([1.165895866372055E+01,-3.157299883744404E+01,6.054513403127080E-03])
velP    = np.array([3.023518143470085E-03,4.354001925477489E-04,-9.190982697153275E-04])*365

## SUN 
massS   = 2E30
posS    = np.array([0, 0, 0])


########## TEST OF ALGORITHMS, TASK 3C ############################
################ Earth - Sun system ###############################

#import Earth_Sun as es
#
#velS    = massE*velE/massS
#
#Earth   = es.planets(velE, posE, massE)
#Sun     = es.planets(velS, posS, massS)          
#
#
#system_verlet = es.system(Earth, Sun)
#solver = es.solver(system_verlet, 2*int(1/10e-5), 10e-5)
#solver.velocity_verlet()
#system_euler = es.system(Earth, Sun)
#solver = es.solver(system_euler, 2*int(1/10e-5), 10e-5)
#solver.euler_fwd()
#
#plt.figure()
#for i in range(len(system_verlet)):
#    plt.plot(system_verlet[i].pos_vec[:,0],system_verlet[i].pos_vec[:,1], linestyle = '--', color = 'Red', label = 'Verlet')
#    plt.plot(system_euler[i].pos_vec[:,0],system_euler[i].pos_vec[:,1], linestyle = '-.', color = 'Blue', label = 'Euler')
#    plt.title('Orbits with Euler and Verlet method')
#    plt.legend()
#    plt.xlabel('x [AU]')
#    plt.ylabel('y [AU]')
#plt.show()

######### Energy/angular momentum conservation ####################

#import Earth_Sun as es
#
#dt = 10.0**np.arange(-1, -7, -1)
#planet = ['Sun', 'Earth']
#energy = []
#angular_momentum = []
#for i in dt:
#    posS    = np.array([0,0,0])
#    velS    = massE*velE/massS
#    Earth   = es.planets(velE, posE, massE)
#    Sun     = es.planets(velS, posS, massS)
#
#    systemES = es.system(Sun, Earth)
#
#    solveES  = es.solver(systemES, int(1/i), i)
#    solveES.velocity_verlet()
#    kin_energy = systemES[1].kinetic_energy(systemES[1].vel_vec)
#    pot_energy = systemES[1].potential_energy(systemES[1].pos_vec)
#    energy.append((kin_energy[-1] + pot_energy[-1])-(kin_energy[0] + pot_energy[0]))
#    l1 = systemES[1].angular_momentum(systemES[1].vel_vec[0,:], systemES[1].pos_vec[0,:])
#    l2 = systemES[1].angular_momentum(systemES[1].vel_vec[-1,:], systemES[1].pos_vec[-1,:])
#    angular_momentum.append(l2 - l1)
#
#methods = ['Euler', 'Verlet']
#
#
### Energy conservation: 
#for method in methods:
#    plt.figure()
#    plt.semilogx(dt, energy)
#    plt.ticklabel_format(axis = 'y', style = 'sci', scilimits=(0,0))
#    plt.grid()
#    plt.xlabel('dt', fontsize = 10)
#    plt.ylabel('E2 - E1', fontsize = 10)
#    plt.title('E2 - E1 with {} method'.format(method), fontsize = 15)
#    plt.savefig('figs/{}_energy.pdf'.format(method))
#    plt.show()
#
#
### Angular momentum: 
#for method in methods:
#    plt.figure()
#    plt.semilogx(dt, np.linalg.norm(angular_momentum, axis = 1))
#    plt.ticklabel_format(axis = 'y', style = 'sci', scilimits=(0,0))
#    plt.grid()
#    plt.xlabel('dt', fontsize = 10)
#    plt.ylabel('L2 - L1', fontsize = 10)
#    plt.title('L2 - L1 with {} method'.format(method), fontsize = 15)
#    plt.savefig('figs/{}_angmom.pdf'.format(method))
#    plt.show()


########## ESCAPE VELOCITY, TASK 3D ###############################

#import Earth_Sun as es
#
#velS      = np.array([0,0,0])
#Sun = es.planets(velS, posS, massS)
#
#plt.figure()
#
#posE = np.array([9.472338196836673E-01,3.216790877655131E-01,-9.310425816193890E-05])
#massE = 6E24
#
#escape_vel = (0.210945/10**3)*np.sqrt((2*6.67E-11*1.989E30)/(149.6E9))*0.84#AU/year
#velE = np.array([-5.428888690270241E-03,escape_vel,-4.491683144318728E-07])
#Earth = es.planets(velE, posE, massE)
#system = es.system(Sun, Earth)
#solve = es.solver(system,  800, 1/365)
#solve.velocity_verlet()
#
#plt.plot(Earth.pos_vec[:,0], Earth.pos_vec[:,1], label = '0.84 times escape velocity')
#
#escape_vel = (0.210945/10**3)*np.sqrt((2*6.67E-11*1.989E30)/(149.6E9))*0.95#AU/year
#velE = np.array([-5.428888690270241E-03,escape_vel,-4.491683144318728E-07])
#Earth = es.planets(velE, posE, massE)
#system = es.system(Sun, Earth)
#solve = es.solver(system,  800, 1/365)
#solve.velocity_verlet()
#
#plt.plot(Earth.pos_vec[:,0], Earth.pos_vec[:,1], label = '0.95 times escape velocity')
#plt.plot(Sun.pos_vec[:,0], Sun.pos_vec[:,1],'yo', linewidth=2, markersize=12)
#
#escape_vel = (0.210945/10**3)*np.sqrt((2*6.67E-11*1.989E30)/(149.6E9))*1.2#AU/year
#velE = np.array([-5.428888690270241E-03,escape_vel,-4.491683144318728E-07])
#Earth = es.planets(velE, posE, massE)
#system = es.system(Sun, Earth)
#solve = es.solver(system,  8000, 1/365)
#solve.velocity_verlet()
#
#plt.plot(Earth.pos_vec[:,0], Earth.pos_vec[:,1], label = '1.2 times escape velocity')
#
#plt.xlabel('x [AU]',fontsize = 10)
#plt.ylabel('y [AU]', fontsize = 10)
#plt.xlim(-6,2)
#plt.ylim(-5, 15)
#plt.title('Trajectory of Earth with changed initial velocity')
#plt.legend()
#plt.savefig('figs/escape_velocity.pdf')
#plt.show()

################# R^BETA ##########################################

#import beta_code as bc
#
#velS      = 0
#
#massE = 6E24
#posE = np.array([9.472338196836673E-01,3.216790877655131E-01,-9.310425816193890E-05])
#velE = np.array([-5.71E-03,1.62E-02, -3.03E-07]) * 365
#Earth = bc.planets(velE, posE, massE)
#Sun = bc.planets(velS, posS, massS)
#system = bc.system(Sun, Earth)
#plt.figure()
#solve = bc.solver(system,  8000, 1/3650, b = 3)
#solve.velocity_verlet()
#plt.plot(Sun.pos_vec[:,0], Sun.pos_vec[:,1], 'yo', linewidth = 2, label = 'Sun')
#plt.plot(Earth.pos_vec[:,0], Earth.pos_vec[:,1], label= 'beta = 2.0')
#
#solve = bc.solver(system,  8000, 1/3650, b = 3.5)
#solve.velocity_verlet()
#plt.plot(Earth.pos_vec[:,0], Earth.pos_vec[:,1],label =  'beta = 2.5')
#
#solve = bc.solver(system,  8000, 1/3650, b = 4)
#solve.velocity_verlet()
#plt.plot(Earth.pos_vec[:,0], Earth.pos_vec[:,1], label = 'beta = 3.0')
#
#plt.xlabel('x [AU]', fontsize = 10)
#plt.ylabel('y [AU]', fontsize = 10)
#plt.title('Trajectory of Earth with r^beta')
#plt.legend()
#plt.savefig('figs/verlet_beta_edition.pdf')
#plt.show()


############ SUN - EARTH - JUPITER, TASK 3E #######################


#import Earth_Jupiter_Sun as ejs
#
#masstimes= 1000 # change to adjust mass of Jupiter
#N = 100000
#dt = 1/3650
#
#
### Sun as Center of mass
#velS    = np.array([0,0,0])
#Earth   = ejs.planets(velE, posE, massE)
#Jupiter = ejs.planets(velJ, posJ, massJ*masstimes)
#Sun     = ejs.planets(velS, posS, massS)
#system  = ejs.system(Sun, Earth, Jupiter)
#solve   = ejs.solver(system, N, dt)
#
#solve.velocity_verlet()
#
#planets = ['Sun', 'Earth', 'Jupiter']
#
#plt.figure()
#plt.title('Mass of Jupiter: {:.1e} kg'.format(massJ), fontsize = 15)
#for i in range(len(system)):
#    plt.plot(system[i].pos_vec[:,0],system[i].pos_vec[:,1], label = planets[i])
#    plt.xlim([-20, 40]) # use when masstimes = 1000
#    plt.ylim([-25, 10])
#    plt.ylabel('y [AU]', fontsize  = 10)
#    plt.xlabel('x [AU]', fontsize  = 10)
#    plt.legend()
#
#
#plt.savefig('figs/threebody_massJ{}.pdf'.format(masstimes))
#plt.show()


################## ALL PLANETS, TASK 3F ###########################

#import solar_system_class as ss
#
#
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
#plt.savefig('figs/solar_system.pdf')

############# PRECESSION OF MERCURY, TASK 3F ######################

# Speed and position of Mercury at perihelion: 
#
posM    = np.array([0.3075,0,0])
velM    = np.array([0,12.44,0])

dt = 1e-7
N = int(100e7)
precession = np.array(ms.velocity_verlet(posM, velM, N, dt))
np.save("precession", precession)
#
#precession = np.load('textfiles/perihelion_50yrs.npy')
#
#
#
#plt.figure()
#plt.plot(np.arctan(precession[:,1]/precession[:,0])*180*3600/np.pi)
#plt.grid()
#plt.xlabel('Perihelion points')
#plt.ylabel('Arc seconds')
#plt.savefig('figs/precession.pdf')

