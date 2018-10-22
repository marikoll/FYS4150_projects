


"""ekstra kode til oppgave om escape velocity


plotte beta """


import numpy as np
import matplotlib.pyplot as plt
import time


class planets:
    def __init__(self, vel0, pos0, mass):
        self.vel0 = vel0        # AU/year
        self.pos0 = pos0        # position in astronomical units
        self.mass = mass/2.10E30
        self.G = 4*np.pi**2
#        self.mass_of_sun = 1
#        self.pos_sun = [0,0,0]
#        self.vel_sun = [0,0,0]

#    def distance(self, pos):
#        self.r = np.linalg.norm(pos,axis=1)
#        return self.r

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
#
#        t_f = 1
#        t = 0
#        h = np.linspace(t, t_f, self.N*t_f)
#        dt = h[1]-h[0]
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

    def plot_escape(pos_v, dt_v, N):
        plt.figure()
        plt.plot(pos_vec[:,0], pos_vec[:,1])
        plt.legend('Sun', 'earth')
        plt.title('escape')
        plt.show()

if __name__ == "__main__":
    posE  = np.array([9.528047055398201E-01,3.053612869840809E-01,-9.272902073041313E-05])

    massM= 3.3E23
    posM =  np.array([-2.138801978256535E-01,-4.028538544608614E-01,-1.397398026440086E-02])
    velM =  np.array([1.927259979627735E-02,-1.164174161133437E-02,-2.720041450765680E-03])*365


    massJ = 1.9E27
    posJ   = np.array([-2.679859418467178E+00,-4.648870533678862E+00,7.922561878424454E-02])
    velJ   = np.array([6.448353939132267E-03,-3.409518873130117E-03,-1.300875035507015E-04])*365


    massSa = 5.5E26
    posSa    = np.array([1.533542535074663E+00,-9.937711138425550E+00,1.117470354327505E-01])
    velSa    = np.array([5.206732001657008E-03,8.336898877431193E-04,-2.220848255804162E-04])*365


    massS     = 2E30
    posS      = np.array([0,0,0])
    velS      = 0
    #-(massSa*velSa + massM*velM + massE*velE + massJ*velJ)/massS

    massE = 6E24
    posE = np.array([9.472338196836673E-01,3.216790877655131E-01,-9.310425816193890E-05])
    velE = np.array([-5.71E-03,1.62E-02, -3.03E-07]) * 365
    Earth = planets(velE, posE, massE)
    Sun = planets(velS, posS, massS)
    system = system(Sun, Earth)
    plt.figure()
    solve = solver(system,  8000, 1/3650, 3)
    solve.velocity_verlet()
    plt.plot(Sun.pos_vec[:,0], Sun.pos_vec[:,1], 'yo', linewidth = 2, label = 'Sun')
    plt.plot(Earth.pos_vec[:,0], Earth.pos_vec[:,1], label= 'beta = 2.0')

    solve = solver(system,  8000, 1/3650, 3.5)
    solve.velocity_verlet()
    plt.plot(Earth.pos_vec[:,0], Earth.pos_vec[:,1],label =  'beta = 2.5')

    solve = solver(system,  8000, 1/3650, 4)
    solve.velocity_verlet()
    plt.plot(Earth.pos_vec[:,0], Earth.pos_vec[:,1], label = 'beta = 3.0')

    plt.xlabel('x [AU]', fontsize = 10)
    plt.ylabel('y [AU]', fontsize = 10)
    plt.title('Trajectory of Earth with r^beta')
    plt.legend()
    plt.savefig('verlet_beta_edition.pdf')
    plt.show()
