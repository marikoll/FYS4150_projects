
import numpy as np
import matplotlib.pyplot as plt
from poisson_solver import poission_jacobi


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
    init_psi = np.zeros(N_x*N_y)
    init_zeta = np.zeros(N_x*N_y)

    for i in range(0, N_x):
        for j in range(0, N_y):
            x = i*dx
            init_psi[i*N_y + j] = np.sin(4.0*np.pi*x) 
            init_zeta[i*N_y + j] = -16.0*np.pi**2*np.sin(4.0*np.pi*x)

    return init_psi, init_zeta


def center(N_x, N_y, dy, dx, T, dt):
    psi_0, zeta_0 = initialize(N_x, N_y, dy, dx)
    alpha = dt/(2*dx)
    gamma = dt/dx
    
    psi_prev = np.zeros(N_x*N_y)
    psi_curr = np.zeros(N_x*N_y)
    zeta_prev = np.zeros(N_x*N_y)
    zeta_pp = np.zeros(N_x*N_y)
    zeta_curr = np.zeros(N_x*N_y)

    bc_0y = np.zeros(N_y)
    bc_Ny = np.zeros(N_y)
    
    bc_0x = np.zeros(N_x)
    bc_Nx = np.zeros(N_x)
    # Arrays for tridiagonal solver
#    diag = np.ones(N_x-2)*(-1)
#    rhs_diag = np.zeros(N_x-2)

    # initial condition and boundary conditions

    for i in range(0, N_x):
        for j in range(0, N_y):
            psi_prev[i*N_y +j] = psi_0[i*N_y +j]
            zeta_pp[i*N_y +j] = zeta_0[i*N_y +j]
            
            psi_curr[0 + j] = bc_0y[j]
            psi_curr[(N_x-1)*N_y + j] = bc_Ny[j]
            zeta_curr[0 + j] = zeta_pp[0 + j] 
            zeta_curr[(N_x-1)*N_y + j] = zeta_pp[(N_x-1)*N_y + j]
        psi_curr[i*N_y + 0] = bc_0x[i]
        psi_curr[i*N_y + N_y-1] = bc_Nx[i]
        zeta_curr[i*N_y + 0] = zeta_pp[i*N_y + 0] 
        zeta_curr[i*N_y + N_y -1] = zeta_pp[i*N_y + N_y -1]


    #initial Euler:
    for i in range(1, N_x-1):
        for j in range(1, N_y-1):
            zeta_prev[i*N_y + j] = zeta_0[i*N_y +j] - \
            alpha*(psi_0[(i+1)*N_y + j] - psi_0[(i-1)*N_y + j])
    
    
    psi_prev = poission_jacobi(zeta_prev, bc_0y, bc_Ny, bc_0x, bc_Nx, dx, dy, N_x, \
                               N_y, 50, psi_prev)
    
#    for i in range(1, N_x-1):
#        rhs_diag[i-1] = -dx2*zeta_prev[i]
#
#
#    psi_prev = tridiag(diag, rhs_diag, N_x-2, psi_prev)
#    print(psi_prev[1:10])
    data_out = np.zeros((int(N_x*N_y+1), int(float(T)/dt)+1))
    t = 0.0
    data_out[0,0] = t
    data_out[1:, 0] = psi_0[:]
    n = 0
    n2 = 1

    while t < T:
        for i in range(1, N_x-1):
            for j in range(1, N_y -1):
                zeta_curr[i*N_y + j] = zeta_pp[i*N_y + j] - \
                gamma*(psi_prev[(i+1)*N_y + j] - psi_prev[(i-1)*N_y + j])
        
        psi_curr = poission_jacobi(zeta_curr, bc_0y, bc_Ny, bc_0x, bc_Nx, dx, dy, N_x, \
                               N_y, 50, psi_curr)

        for i in range(1, N_x -1):
            for j in range(1, N_y -1):
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
    T = 100
    dt = 0.01

    dx = 1.0/40
    dy = 1.0/40
    L = 1.0
    N_x = int(L/dx + 1)
    N_y = int(L/dy +1)
    
    
    data_out = center(N_x, N_y, dy, dx, T, dt)
    
    
    
    
    
    
    
    
