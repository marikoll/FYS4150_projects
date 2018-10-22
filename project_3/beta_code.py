


"""
To run, use python beta_code.py

Project 3,
changing the initial algorithm to include changes of the gravitational force.
The solver class now needs beta as input, and the acceleration can be changed easily.

At the end of the code, the class and functions are called and a plot is produced
with different values of beta.
"""


import numpy as np
import matplotlib.pyplot as plt
import time


class planets:
    def __init__(self, vel0, pos0, mass):
        self.vel0 = vel0        # AU/year
        self.pos0 = pos0        # position in astronomical units
        self.mass = mass/2.10E30
        self.G = 4*np.pi**2

    def kinetic_energy(self, vel):
        KE =  0.5*self.mass*np.linalg.norm(vel,axis=1)**2
        return KE

    def potential_energy(self, r, mass_of_sun=1):
        PE = -self.G*mass_of_sun*self.mass/np.linalg.norm(r, axis=1)#self.distance(pos)
        return PE

    def angular_momentum(self, vel, r):
        L = np.cross(r, vel*self.mass,axis=1)
        return L

    def make_pos_vec(self, N):
        self.pos_vec = np.zeros((N,3))
        self.pos_vec[0] = self.pos0

    def make_vel_vec(self,N):
        self.vel_vec = np.zeros((N,3))
        self.vel_vec[0] = self.vel0

    def make_acc_vec(self,N,acc0):
        self.acc_vec = np.zeros((N,3))
        self.acc_vec[0] = acc0



def system(*args):
    system = []
    for i in args:
        system.append(i)

    return system


class solver:
    def __init__(self, system, N, dt,b):
        self.dt = dt
        self.N = N
        self.system = system
        self.G = 4*np.pi**2
        self.b = b

        for i in self.system:
            i.make_pos_vec(N)
            i.make_vel_vec(N)
        for i in self.system:
            i.make_acc_vec(N, self.acceleration(i, 0, b))

    def acceleration(self, current, i, b):
        acc = np.zeros(3)

        for planet in self.system:
            if planet != current:
                r = current.pos_vec[i] - planet.pos_vec[i]
                acc -= r*self.G*planet.mass/(np.linalg.norm(r)**b)

        return acc


    def velocity_verlet(self):

        b = self.b
        dt = self.dt
        dt2 = dt**2


        start = time.time()

        for i in range(self.N-1):
            for planet in self.system:
                planet.pos_vec[i+1] = planet.pos_vec[i] + dt*planet.vel_vec[i]\
                + 0.5*dt2*planet.acc_vec[i]
            for planet in self.system:
                planet.acc_vec[i+1] = self.acceleration(planet, i+1, b)
            for planet in self.system:
                planet.vel_vec[i+1] = planet.vel_vec[i] + \
                0.5*dt*(planet.acc_vec[i] + planet.acc_vec[i+1])
        stop = time.time()
        print("Verlet: {:.3f} sec".format(stop-start))


if __name__ == "__main__":

    # As we need the Sun to be the mass center, its position is set in origo
    # We assume that the velocity of the Sun is negelectable compared to the Earth
    massS     = 2E30
    posS      = np.array([0,0,0])
    velS      = 0

    massE = 6E24
    posE  = np.array([9.528047055398201E-01,3.053612869840809E-01,-9.272902073041313E-05])
    posE = np.array([9.472338196836673E-01,3.216790877655131E-01,-9.310425816193890E-05])
    velE    = np.array([-7.644048786784979E-03,1.548720517812966E-02,9.799447004349826E-08])*365

    Earth = planets(velE, posE, massE)
    Sun = planets(velS, posS, massS)
    system = system(Sun, Earth)

    # Beta = 3
    plt.figure()
    solve = solver(system,  8000, 1/3650, 3)
    solve.velocity_verlet()
    plt.plot(Sun.pos_vec[:,0], Sun.pos_vec[:,1], 'yo', linewidth = 2, label = 'Sun')
    plt.plot(Earth.pos_vec[:,0], Earth.pos_vec[:,1], label= 'beta = 2.0')

    # Beta = 3.5
    solve = solver(system,  8000, 1/3650, 3.5)
    solve.velocity_verlet()
    plt.plot(Earth.pos_vec[:,0], Earth.pos_vec[:,1],label =  'beta = 2.5')

    # Beta = 4.0
    solve = solver(system,  8000, 1/3650, 4)
    solve.velocity_verlet()
    plt.plot(Earth.pos_vec[:,0], Earth.pos_vec[:,1], label = 'beta = 3.0')

    plt.xlabel('x [AU]', fontsize = 10)
    plt.ylabel('y [AU]', fontsize = 10)

    plt.xlim(-1,3)
    plt.ylim(-1.5, 1)
    plt.title('Trajectory of Earth with r^beta')
    plt.legend()
    plt.savefig('verlet_beta_edition.pdf')
    plt.show()
