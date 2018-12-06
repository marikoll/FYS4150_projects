
import numpy as np
import matplotlib.pyplot as plt
from poisson_solver import poisson_jacobi_periodic


def periodic_matrix(n_rows, n_cols):
    A = np.zeros((n_rows, n_cols))
    for i in range(0,n_rows):
        for j in range(0,n_cols):
            if i == j:
                A[i, j] = 2.0
            elif abs(i-j) ==1:
                A[i,j] = -1.0
    A[0, n_cols-1] = -1.0
    A[n_rows -1, 0] = -1.0

    return A

def initialize(N_x, N_y, dx, dy):
    init_psi = np.zeros(N_x*N_y)
    init_zeta = np.zeros(N_x*N_y)

    for i in range(0, N_x):
        for j in range(0, N_y):
            x = i*dx
            init_psi[i*N_y + j] = np.sin(4.0*np.pi*x) 
            init_zeta[i*N_y + j] = -16.0*np.pi**2*np.sin(4.0*np.pi*x)

    return init_psi, init_zeta


def leapfrog(N_x,N_y, dx,dy, T, dt):
    psi_0, zeta_0 = initialize(N_x, N_y, dx, dy)

    alpha = dt/(2*dx)
    gamma =  dt/dx
    dx2 = dx**2


    multi = N_x*N_y
    psi_prev = np.zeros(multi)
    psi_curr = np.zeros(multi)
    zeta_prev = np.zeros(multi)
    zeta_pp = np.zeros(multi)
    zeta_curr = np.zeros(multi)
    
    rhs_poisson = np.zeros(multi-1)
    A = periodic_matrix(multi-1, multi-1)

    for i in range(0,N_x):
        for j in range(0,N_y):
            psi_prev[i*N_y + j] = psi_0[i*N_y +j]
            zeta_prev[i*N_y + j] = zeta_0[i*N_y +j]

    #initial Euler:
    for i in range(1, N_x-1):
        for j in range(0, N_y):
            zeta_prev[i*N_y +j] = zeta_0[i*N_y + j] - \
            alpha*(psi_0[(i+1)*N_y + j] - psi_0[(i-1)*N_y + j])
    for j in range(0, N_y):
        zeta_prev[0*N_y + j] = zeta_0[0*N_y + j] - alpha*(psi_0[1*N_y + j] \
                 - psi_0[(N_x -2)*N_y + j])
        zeta_prev[(N_x-1)*N_y + j] = zeta_prev[0*N_y + j]

    for i in range(0, N_x-1):
        rhs_poisson[i] = -dx2*zeta_prev[i]

    psi_prev = np.linalg.solve(A, rhs_poisson)

    data_out = np.zeros((N_x*N_y-1, int(float(T)/dt)+1))
    t = 0.0
    data_out[0,0] = t
    data_out[1:, 0] = psi_0[:-1]
    n = 0
    n2 = 1

    while t < T:
        #forward Euler:
        for i in range(1, N_x-2):
            for j in range(0,N_y-1):
                zeta_curr[i*N_y + j] = zeta_pp[i*N_y + j] \
                - gamma*(psi_prev[(i+1)*N_y + j] - psi_prev[(i-1)*N_y + j])
        for j in range(0, N_y):
            zeta_curr[0*N_y + j] = zeta_pp[0*N_y + j] - \
            gamma*(psi_prev[1*N_y + j] - psi_prev[(N_x -2)*N_y + j])
            zeta_curr[(N_x-1)*N_y + j] = zeta_curr[0*N_y + j]

        for i in range(0, N_x-1):
            rhs_poisson[i] = -dx2*zeta_curr[i]

        psi_curr = np.linalg.solve(A, rhs_poisson)
        #poisson_jacobi_periodic(zeta_curr, dx, dt, N_x. N_y, 50, psi_curr)

        for i in range(0, N_x-1):
            for j in range(0, N_y-1):

                psi_prev[i*N_y + j] = psi_curr[i*N_y + j]
                zeta_pp[i*N_y + j] = zeta_prev[i*N_y + j]
                zeta_prev[i*N_y + j] = zeta_curr[i*N_y + j]

        t += dt
        if (n % 50 == 0):
            data_out[0, n2] = t
            data_out[1:, n2] = psi_curr[:]
            
            n2 += 1

        n += 1

    return data_out


if __name__ == "__main__":

    T = 200
    dt = 0.01

    dx = 1.0/40
    dy = 1.0/40
    L = 1.0
    N_x = int(L/dx + 1)
    N_y = int(L/dy +1)


    outstuff2 = leapfrog(N_x, N_y, dx,dy, T, dt)

#    psiE_gauss = euler(init_psi_gauss, init_zeta_gauss, N, dx, T, dt)
#    psiLF_gauss = leapfrog(init_psi_gauss, init_zeta_gauss, N, dx, T, dt)

    #x = np.linspace(0, 1, N-1)

    #plt.figure()
    #plt.plot(x, outstuff[:,0], 'r-')
    #plt.plot(x, outstuff2[:,0], 'b-.')


#    plt.figure()
#    plt.plot(x, psiE_gauss[1:N-3], 'r-')
#    plt.plot(x, psiLF_gauss[1:N-3], 'b-.')
