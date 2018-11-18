"""
Two - dimensional Ising model 
"""


import numpy as np
import numba



@numba.njit(cache = True)
def initialization(spins):
    """
    Function that computes the initial energy and magnetization based on 
    the configuration of the input matrix. 
    
    Input: 
        spins   - <np.ndarray> matrix of sixe (spins x spins)
    
    Output: 
        E       - <float64> energy of initial configuration
        M       - <float64> magnetization of initial configuration
    """
    E = 0
    M = 0
    num_spins = len(spins)
    
    for i in range(num_spins):
        for j in range(num_spins):
            left = spins[i-1, j] if i>0 else spins[num_spins - 1, j]
            above = spins[i, j-1] if j>0 else spins[i, num_spins - 1]
            
            E -= spins[i,j]*(left+above)
            M += spins[i, j]
            
    return E, M


@numba.njit(cache=True)
def MC(spins, num_cycles, temperature):
    """
    Monte-Carlo algorithm for solving the Ising model with implemented Markov
    chain and Metropolis algorithm.
    
    Input: 
        spins       - <int> matrix of sixe (spins x spins)
        num_cycles  - <int> number of Monte-Carlo cycles
        temperature - <float> temperature of the system in units JT/k_B
        
    Output: 
        energy_avg      - <float64> expectation values of the energy
        magnet_avg      - <float64> expectation values of the  magnetization
        C_v             - <float> heat capacity
        susceptibility  - <float> susceptibility
        abs_magnet      - <float64> expectation values of the absolute 
                                    magnetization
        counter_list    - <float64> number of accepted states
    """
    num_spins = len(spins)

    exp_values = np.zeros((int(num_cycles), 5))
    
    E, M = initialization(spins)
    counter_list = np.zeros(num_cycles)
    counter = 0
    for i in range(num_cycles):
        for j in range(num_spins**2):
            ix = np.random.randint(num_spins)
            iy = np.random.randint(num_spins)
    
            left = spins[ix - 1, iy] if ix > 0 else spins[num_spins - 1, iy]
            right = spins[ix + 1, iy] if ix < (num_spins - 1) else spins[0, iy]
    
            above = spins[ix, iy - 1] if iy > 0 else spins[ix, num_spins - 1]
            below = spins[ix, iy + 1] if iy < (num_spins - 1) else spins[ix, 0]
    
            delta_energy = (2 * spins[ix, iy] * (left + right + above + below))
            
            # Metropolis: 
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
        
    energy_avg = np.cumsum(exp_values[:,0])/np.arange(1, num_cycles +1)
    magnet_avg = np.cumsum(exp_values[:,1])/np.arange(1, num_cycles +1)

    energy2_avg = np.cumsum(exp_values[:,2])/np.arange(1, num_cycles +1)
    magnet2_avg = np.cumsum(exp_values[:,3])/np.arange(1, num_cycles +1)
    magnet_absavg = np.cumsum(exp_values[:,4])/np.arange(1, num_cycles +1)
    energy_var = (energy2_avg[-1] - energy_avg[-1]**2)/(num_spins**2)
    magnet_var = (magnet2_avg[-1] - magnet_avg[-1]**2)/(num_spins**2)
    
    energy_avg = energy_avg/num_spins**2
    magnet_avg = magnet_avg/num_spins**2
    C_v = energy_var/temperature**2
    susceptibility = magnet_var/temperature
    abs_magnet = magnet_absavg/num_spins**2


    return energy_avg, magnet_avg, C_v, susceptibility, abs_magnet, counter_list


@numba.njit(cache=True)
def MC_cutoff(spins, num_cycles, temperature, P):
    """
    Monte-Carlo algorithm for solving the Ising model with implemented Markov
    chain and Metropolis algorithm. With cutoff percentage meaning that the 
    algorithm will not sample until the percentage of cycles are performed
    
    Input: 
        spins       - <int> matrix of sixe (spins x spins)
        num_cycles  - <int> number of Monte-Carlo cycles
        temperature - <float> temperature of the system in units JT/k_B
        P           - <float> cutoff percentage
        
    Output: 
        energy          - <float64> energy of the system
        energy_var      - <float> variance of the system
        energy_avg      - <float64> expectation values of the energy
        magnet_avg      - <float64> expectation values of the  magnetization
        C_v             - <float> heat capacity
        susceptibility  - <float> susceptibility
        abs_magnet      - <float64> expectation values of the absolute 
                                    magnetization
    """
    num_spins = len(spins)

    exp_values = np.zeros((int(num_cycles), 5))
    sampling_starts_from = int(num_cycles*(P))
    
    E, M = initialization(spins)

    for i in range(num_cycles):
        for j in range(num_spins**2):
            ix = np.random.randint(num_spins)
            iy = np.random.randint(num_spins)
    
            left = spins[ix - 1, iy] if ix > 0 else spins[num_spins - 1, iy]
            right = spins[ix + 1, iy] if ix < (num_spins - 1) else spins[0, iy]
    
            above = spins[ix, iy - 1] if iy > 0 else spins[ix, num_spins - 1]
            below = spins[ix, iy + 1] if iy < (num_spins - 1) else spins[ix, 0]
    
            delta_energy = (2 * spins[ix, iy] * (left + right + above + below))
    
            
            if np.random.random() <= np.exp(-delta_energy / temperature):
                spins[ix, iy] *= -1.0
                
    
                E += delta_energy
                M += 2*spins[ix, iy]
    
            if i >= sampling_starts_from:
                exp_values[i,0] = E
                exp_values[i,1] = M
                exp_values[i,2] = E**2
                exp_values[i,3] = M**2
                exp_values[i,4] = np.abs(M)

        
    energy = exp_values[:,0]/num_spins**2
    energy_avg = np.cumsum(exp_values[:,0])/np.arange(1, num_cycles +1)
    magnet_avg = np.cumsum(exp_values[:,1])/np.arange(1, num_cycles +1)
    energy2_avg = np.cumsum(exp_values[:,2])/np.arange(1, num_cycles +1)
    magnet2_avg = np.cumsum(exp_values[:,3])/np.arange(1, num_cycles +1)
    magnet_absavg = np.cumsum(exp_values[:,4])/np.arange(1, num_cycles +1)
    energy_var = (energy2_avg[-1] - energy_avg[-1]**2)/(num_spins**2)
    magnet_var = (magnet2_avg[-1] - magnet_avg[-1]**2)/(num_spins**2)
    
    energy_avg = energy_avg/num_spins**2
    magnet_avg = magnet_avg/num_spins**2
    C_v = energy_var/temperature**2
    susceptibility = magnet_var/temperature
    abs_magnet = magnet_absavg/num_spins**2

    return energy, energy_var, energy_avg, magnet_avg, C_v, susceptibility, abs_magnet

    
if __name__ == "__main__": 
    pass
