import numpy as np
import numba


def periodic(i, limit, add):
    return (i + limit + add) % limit


@numba.njit(cache = True)
def initial_energy(spins, temp):
    E = 0
    M = 0
    num_spins = len(spins)
    
    for i in range(num_spins):
        for j in range(num_spins):
#            if (temp < 1.5): 
#                spins[i, j] = 1
            left = spins[i-1, j] if i>0 else spins[num_spins - 1, j]
            above = spins[i, j-1] if j>0 else spins[i, num_spins - 1]
            
            E -= spins[i,j]*(left+above)
            M += spins[i, j]
            
    return E, M


@numba.njit(cache=True)
def MC(spins, num_cycles, temperature, ordered = False):#, num_thermalization_steps=0):
    num_spins = len(spins)
    exp_values = np.zeros((int(num_cycles), 5))
    
    E, M = initial_energy(spins, temperature)
    counter_list = np.zeros(num_cycles)
    counter = 0
    for i in range(num_cycles):
        ix = np.random.randint(num_spins)
        iy = np.random.randint(num_spins)

        left = spins[ix - 1, iy] if ix > 0 else spins[num_spins - 1, iy]
        right = spins[ix + 1, iy] if ix < (num_spins - 1) else spins[0, iy]

        above = spins[ix, iy - 1] if iy > 0 else spins[ix, num_spins - 1]
        below = spins[ix, iy + 1] if iy < (num_spins - 1) else spins[ix, 0]

        delta_energy = (2 * spins[ix, iy] * (left + right + above + below))
        
 #       if i > int(num_cycles*0.1):
        if np.random.random() <= np.exp(-delta_energy / temperature):
            spins[ix, iy] *= -1.0
            E += delta_energy
            M += 2*spins[ix, iy]
            counter += 1
        
        exp_values[i,0] = E
        exp_values[i,1] = M
        exp_values[i,2] = E**2
        exp_values[i,3] = M**2
        exp_values[i,4] = np.abs(M)
        counter_list[i] = counter
        

    norm = 1/float(num_cycles)
    
    energy_avg = np.sum(exp_values[:,0])*norm
    magnet_avg = np.sum(exp_values[:,1])*norm
    energy2_avg = np.sum(exp_values[:,2])*norm
    magnet2_avg = np.sum(exp_values[:,3])*norm
    magnet_absavg = np.sum(exp_values[:,4])*norm
    
    energy_var = (energy2_avg - energy_avg**2)/(num_spins**2)#**2*temperature**2)
    magnet_var = (magnet2_avg - magnet_absavg**2)/(num_spins**2)#**2*temperature**2)
    
    energy_avg = energy_avg/num_spins**2
    magnet_avg = magnet_avg/num_spins**2
    C_v = energy_var/temperature**2
    susceptibility = magnet_var/temperature
    abs_magnet = magnet_absavg/num_spins**2

    return energy_avg, magnet_avg, C_v, susceptibility, abs_magnet, counter_list

    
if __name__ == "__main__": 
    pass
#    spins       = 2
#    trials      = int(1e1)#, int(1e3), int(1e4), int(1e5), int(1e6), int(1e7)]
#    temp = 1.0
#    grid = np.random.choice([-1,1],size=(spins, spins))#np.ones((spins, spins))
#    energy_avg, magnet_avg, C_v, susceptibility, abs_magnet, c= MC(grid, trials, temp)
#    sampled_energies = np.zeros(len(trials))
#    sampled_magnets = np.zeros(len(trials))
#    sampled_cv = np.zeros(len(trials))
#    sampled_suscept = np.zeros(len(trials))
#    sampled_absmagn = np.zeros(len(trials))
#    
#    
#    
#    for i in range(len(trials)):
#        grid = np.ones((spins, spins))
#        energy_avg, magnet_avg, C_v, susceptibility, abs_magnet, c= MC(grid, trials[i], temp)#, w)
#        sampled_energies[i] = energy_avg
#        sampled_magnets[i] = magnet_avg
#        sampled_cv[i] = C_v
#        sampled_suscept[i] = susceptibility
#        sampled_absmagn[i] = abs_magnet