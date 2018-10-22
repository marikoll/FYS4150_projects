""" Extra code for escape velocity """




import numpy as np
import solar_system_class as ss
import matplotlib.pyplot as plt



massS     = 2E30
posS      = np.array([0,0,0])
velS      = 0
Sun = ss.planets(velS, posS, massS)

#print(Earth.vel_vec)

plt.figure()


#plt.plot(Earth.pos_vec[:,0], Earth.pos_vec[:,1], label = '0.8 times escape velocity')
posE = np.array([9.472338196836673E-01,3.216790877655131E-01,-9.310425816193890E-05])
massE = 6E24

escape_vel = (0.210945/10**3)*np.sqrt((2*6.67E-11*1.989E30)/(149.6E9))*0.84#AU/year
vel = [-5.428888690270241E-03,escape_vel,-4.491683144318728E-07]
velE = np.array(vel)
Earth = ss.planets(velE, posE, massE)
system = ss.system(Sun, Earth)
solve = ss.solver(system,  800, 1/365)
solve.velocity_verlet()

plt.plot(Earth.pos_vec[:,0], Earth.pos_vec[:,1], label = '0.84 times escape velocity')

escape_vel = (0.210945/10**3)*np.sqrt((2*6.67E-11*1.989E30)/(149.6E9))*0.95#AU/year
vel = [-5.428888690270241E-03,escape_vel,-4.491683144318728E-07]
velE = np.array(vel)
Earth = ss.planets(velE, posE, massE)
system = ss.system(Sun, Earth)
solve = ss.solver(system,  800, 1/365)
solve.velocity_verlet()

plt.plot(Earth.pos_vec[:,0], Earth.pos_vec[:,1], label = '0.95 times escape velocity')
plt.plot(Sun.pos_vec[:,0], Sun.pos_vec[:,1],'yo', linewidth=2, markersize=12)

escape_vel = (0.210945/10**3)*np.sqrt((2*6.67E-11*1.989E30)/(149.6E9))*1.2#AU/year
vel = [-5.428888690270241E-03,escape_vel,-4.491683144318728E-07]
velE = np.array(vel)
Earth = ss.planets(velE, posE, massE)
system = ss.system(Sun, Earth)
solve = ss.solver(system,  8000, 1/365)
solve.velocity_verlet()

plt.plot(Earth.pos_vec[:,0], Earth.pos_vec[:,1], label = '1.2 times escape velocity')

plt.xlabel('x [AU]',fontsize = 10)
plt.ylabel('y [AU]', fontsize = 10)
plt.xlim(-6,2)
plt.ylim(-5, 15)
plt.title('Trajectory of Earth with changed initial velocity')
plt.legend()
plt.savefig('escape_velocity.pdf')
plt.show()
