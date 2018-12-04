
import numpy as np
import matplotlib.pyplot as plt



def assertion(init_psi, init_zeta,N_x, N_y):
    epsilon = 1e-10
    for j in range(0, N_y):

        if abs(init_psi[0 + j] - init_psi[(N_x -1)*N_y + j] > epsilon)
            print('psi_0: ', init_psi[0], 'psi_N:', init_psi[N_x-1])
            print('Error, initial condition does not satisfy BC')
    for i in range(0, N_x):
        if abs(init_psi[i*N_y + 0] - init_psi[(i*N_y)+(N_y -i)] > epsilon)
            print('psi_0: ', init_psi[0], 'psi_N:', init_psi[N_x-1])
            print('Error, initial condition does not satisfy BC')


    psi_0 = np.empty(len(init_psi))
    zeta_0 = np.empty(len(init_zeta))


    psi_0 = init_psi
    zeta_0 = init_zeta
    return psi_0, zeta_0

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


def leapfrog(init_psi, init_zeta, N_x,N_y, dx, T, dt):
    psi_0, zeta_0 = assertion(init_psi, init_zeta, N_x, N_y)

    alpha = dt/(2*dx)
    gamma =  dt/dx

    multi = N_x*N_y
    psi_prev = np.zeros(multi)
    psi_curr = np.zeros(multi)
    zeta_prev = np.zeros(multi)
    zeta_pp = np.zeros(multi)
    zeta_curr = np.zeros(multi)

    for i in range(0,N_x):
        for j in range(0,N_y):

            psi_prev[i*N_y + j] = psi_0[i*N_y +j]
            zeta_prev[i*N_y + j] = zeta_0[i*N_y +j]

    #initial Euler:
    for i in range(1, N_x-1):
        for j in range(0, N_y)
            zeta_prev[i*N_y +j] = zeta_0[i*N_y + j] \
            - alpha*(psi_0[(i+1)*N_y + j] - psi_0[(i-1)*N_y + j])
    for j in range(0, N_y):
        zeta_prev[0*N_y + j] = zeta_0[0*N_y + j] - \
        alpha*(psi_0[1*N_y + j] - psi_0[(N_x -2)*N_y + j])
        zeta_prev[(N_x-1)*N_y + j] = zeta_prev[0*N_y + j]


    SOME SOLVER HERE


    outstuff = np.zeros((N_x-1, int(float(T)/dt)+1))
    t = 0.0
    n = 2

    while t < T:
        #forward Euler:
        for i in range(1, N_x-1):
            for j in range(0,N_y):
                zeta_curr[i*N_y + j] = zeta_pp[i*N_y + j] /
                - gamma*(psi_prev[(i+1)*N_y + j] - psi_prev[(i-1)*N_y + j])
        for j in range(0, N_y):
            zeta_curr[0*N_y + j] = zeta_pp[0*N_y + j] - gamma*(psi_prev[1*N_y + j] - psi_prev[(N_x -2])*N_y + j)
            zeta_curr[(N_x-1)*N_y + j] = zeta_curr[0*N_y + j]

        SOME SOLVER HERE

        for i in range(0, N_x-1):
            for j in range(0, N_y0):

                psi_prev[i*N_y + j] = psi_curr[i*N_y + j]
                zeta_pp[i*N_y + j] = zeta_prev[i*N_y + j]
                zeta_prev[i*N_y + j] = zeta_curr[i*N_y + j]

        t += dt
        if (n % 50 == 0):
            outstuff[:, n] = psi_curr[:]

        n += 1

    return outstuff


if __name__ == "__main__":

    T = 150
    dt = 0.5

    dx = 1.0/40
    L = 1.0
    N = int(L/dx + 1)

    init_psi = np.zeros(N)
    init_zeta = np.zeros(N)

    init_psi_gauss = np.zeros(N)
    init_zeta_gauss = np.zeros(N)
    sigma = 0.1

    for i in range(0, N-1):
        x = i*dx
        init_psi[i] = np.sin(4.0*np.pi*x)
        init_zeta[i] = -16.0*np.pi**2*np.sin(4.0*np.pi*x)

#        init_psi_gauss[i] = np.exp(-((x-0.5)/sigma)**2)
#        init_zeta_gauss[i] = (4*((x-0.5)/sigma)**2) - (2/sigma**2)*(np.exp(-((x-0.5)/sigma)**2))
#

    outstuff= euler(init_psi, init_zeta, N, dx, T, dt)
    outstuff2 = leapfrog(init_psi, init_zeta, N, dx, T, dt)

#    psiE_gauss = euler(init_psi_gauss, init_zeta_gauss, N, dx, T, dt)
#    psiLF_gauss = leapfrog(init_psi_gauss, init_zeta_gauss, N, dx, T, dt)

    x = np.linspace(0, 1, N-1)

    plt.figure()
    plt.plot(x, outstuff[:,0], 'r-')
    plt.plot(x, outstuff2[:,0], 'b-.')


#    plt.figure()
#    plt.plot(x, psiE_gauss[1:N-3], 'r-')
#    plt.plot(x, psiLF_gauss[1:N-3], 'b-.')
