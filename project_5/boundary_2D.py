
import numpy as np
import matplotlib.pyplot as plt


def tridiag(b, y, N, soltn):

    b[0] = 2.0

    #forward substitution
    for i in range(1, N):
         b[i] = float(i + 2)/float(i + 1)
         y[i] = y[i] + (float(y[i-1])/float(b[i-1]))

    #backward substitution
    soltn[N-1] = float(y[N-1])/float(b[N-1])

    for i in range(N-2, 0, -1):
        soltn[i] = float(y[i] + soltn[i+1])/float(b[i])

    return soltn

def initialize(N_x, N_y, dx, dy):
    init_psi = np.zeros(N_y*N_x)
    init_zeta = np.zeros(N_y*N_x)

    for i in range(0, N-1):
        x = i*dx
        init_psi[i*N_y + j] = np.sin(4.0*np.pi*x)
        init_zeta[i*N_y + j] = -16.0*np.pi**2*np.sin(4.0*np.pi*x)

    return init_psi, init_zeta


def center(N_x, N_y, dy, dx, T, dt):
    psi_0, zeta_0 = initialize(N_x, N_y, dy, dx)
    alpha = dt/(2*dx)
    dx2 = dx**2
    bc_0 = 0.0
    bc_N = 0.0
    gamma = dt/dx

    psi_prev = np.zeros(N_x*N_y)
    psi_curr = np.zeros(N_x*N_y)
    zeta_prev = np.zeros(N_x*N_y)
    zeta_pp = np.zeros(N_x*N_y)
    zeta_curr = np.zeros(N_x*N_y)


    # Arrays for tridiagonal solver
    diag = np.ones(N_x-2)*(-1)
    rhs_diag = np.zeros(N_x-2)

    # initial condition and boundary conditions

    for i in range(0,N_x):
        for j in range(0, N_y):
            
    psi_prev  = psi_0
    zeta_pp = zeta_0

    psi_curr[0] = bc_0; psi_curr[N_x-1] = bc_N
    zeta_curr[0] = zeta_pp[0]; zeta_curr[N_x-1] = zeta_pp[N_x-1]


    #initial Euler:
    for i in range(1, N_x-1):
        zeta_prev[i] = zeta_0[i] - alpha*(psi_0[i+1] - psi_0[i-1])
#    for i in range(1, N_x-1):
#        rhs_diag[i-1] = -dx2*zeta_prev[i]
#
#
#    psi_prev = tridiag(diag, rhs_diag, N_x-2, psi_prev)
#    print(psi_prev[1:10])
    out_data = np.zeros((N_x, int(float(T)/dt)))
    t = 0.0
    n = 0

    while t < T:
        for i in range(1, N_x-1):
            zeta_curr[i] = zeta_pp[i] - gamma*(psi_prev[i+1] - psi_prev[i-1])
        for i in range(1,N_x-1):
            rhs_diag[i-1] = -dx2*zeta_curr[i]
        psi_curr = tridiag(diag, rhs_diag, N_x -2, psi_curr)

        for i in range(1, N_x -1):
            psi_prev[i] = psi_curr[i]
            zeta_pp[i] = zeta_prev[i]
            zeta_prev[i] = zeta_curr[i]
        t += dt
        if (n % 20 == 0):
            out_data[:, n] = psi_curr[:]
        n += 1

    return out_data
