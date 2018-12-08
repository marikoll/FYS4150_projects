
import numpy as np
import matplotlib.pyplot as plt
from poisson_solver import poisson_jacobi_matrix


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

def initialize(N_x, N_y, dx, dy, case):
    init_psi = np.zeros((N_x,N_y))
    init_zeta = np.zeros((N_x,N_y))
    
    if case == 'sine':
        for i in range(0, N_x):
            for j in range(0, N_y):
                x = i*dx
                y = i*dy
                init_psi[i, j] = np.sin(np.pi*y)*np.sin(4.0*np.pi*x) 
                init_zeta[i, j] = -17.0*np.pi**2*np.sin(4.0*np.pi*x)*np.sin(np.pi*y)
                if i==0 or j == 0:
                    init_psi[i, j] = 0.0
                    init_zeta[i, j] = 0.0
                if i==N_x-1 or j == N_y-1:
                    init_psi[i, j] = 0.0
                    init_zeta[i, j] = 0.0
    if case == 'gauss':
        for i in range(0, N_x):
            for j in range(0, N_y):
                x = i*dx
                y = i*dy
                sigma = 0.1
                init_psi[i] = np.exp(-((x-0.5)/sigma)**2)
                init_zeta[i] = 4.0*((x-0.5)/sigma**2)**2 - (2/sigma**2)*np.exp(-((x-0.5)/sigma)**2)
    
#    init_psi = init_psi.reshape(N_x*N_y)
#    init_zeta = init_zeta.reshape(N_x*N_y)
    return init_psi, init_zeta


def center(N_x, N_y, dy, dx, T, dt, case):
    psi_0, zeta_0 = initialize(N_x, N_y, dy, dx, case)
    alpha = dt/(2*dx)
    gamma = dt/dx
    
    psi_prev = psi_0
    psi_curr = np.zeros((N_x,N_y))
    zeta_prev = np.zeros((N_x,N_y))
    zeta_pp = zeta_0
    zeta_curr = np.zeros((N_x,N_y))

    bc_0y = np.zeros(N_y)
    bc_Ny = np.zeros(N_y)
    
    bc_0x = np.zeros(N_x)
    bc_Nx = np.zeros(N_x)

    for i in range(0, N_x):
        for j in range(0, N_y):
#            psi_prev[i*N_y +j] = psi_0[i*N_y +j]
#            zeta_pp[i*N_y +j] = zeta_0[i*N_y +j]
#            
            psi_curr[0, j] = bc_0y[j]
            psi_curr[N_x-1, j] = bc_Ny[j]
            zeta_curr[0, j] = zeta_pp[0, j] 
            zeta_curr[N_x-1, j] = zeta_pp[N_x-1, j]
        psi_curr[i, 0] = bc_0x[i]
        psi_curr[i, N_y-1] = bc_Nx[i]
        zeta_curr[i, 0] = zeta_pp[i, 0] 
        zeta_curr[i, N_y -1] = zeta_pp[i, N_y -1]


    #initial Euler:
    for i in range(1, N_x-1):
        for j in range(1, N_y-1):
            zeta_prev[i, j] = zeta_0[i,j] - alpha*(psi_0[i+1, j] - psi_0[i-1, j])
    
    
    psi_prev = poisson_jacobi_matrix(psi_prev, zeta_prev, N_x, dx)
    

    data_out = np.zeros((N_x, N_y, int(float(T)/dt)+1))
    t = 0.0
#    data_out[0,0, 0] = t
    data_out[:, :,0] = psi_0
    n = 0
    n2 = 1

    while t < T:
        for i in range(1, N_x-1):
            for j in range(1, N_y -1):
                zeta_curr[i, j] = zeta_pp[i, j] - gamma*(psi_prev[i+1, j] - psi_prev[i-1, j])
        
        psi_curr = poisson_jacobi_matrix(zeta_curr, psi_curr, N_x, dx)

        for i in range(1, N_x -1):
            for j in range(1, N_y -1):
                psi_prev[i, j] = psi_curr[i, j]
                zeta_pp[i, j] = zeta_prev[i, j]
                zeta_prev[i, j] = zeta_curr[i, j]
        t += dt
        if (n % 50 == 0):
            data_out[:, :, n2] = psi_curr
            
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
    
    
    data_out_sine = center(N_x, N_y, dy, dx, T, dt, 'sine')
    
#    t = data_out[0,:]
#    psi_sine = data_out_sine[1:, :200]
#    new_psi_sine = np.zeros((41, 41, 200))
#    for t in range(0, 200):
#        new_psi_sine[:,:, t] = psi_sine[:, t].reshape(41,41).transpose()
#    
#    x = np.linspace(0, 1, 41)
#    y = np.linspace(0, 1, 41)
#    
#    plt.style.use("ggplot")
#    fig = plt.figure(figsize = (9,7))
#    CS = plt.contourf(x, y, new_psi_sine[:,:,0], 20, cmap = plt.cm.RdBu_r)
#    plt.colorbar(CS, orientation = "vertical")
#    plt.title(r'Contour field of $\psi(x, y, 0)$ in bounded domain', fontsize = 15)
#    plt.xlabel('x', fontsize = 13)
#    plt.ylabel('y', fontsize = 13)
#    plt.savefig('figs/sine_boundary_2d.pdf', bbox_inches = 'tight')
    
    
#    data_out_gauss = center(N_x, N_y, dy, dx, T, dt, 'gauss')
#    
##    t = data_out[0,:]
#    psi_gauss = data_out_gauss[1:, :200]
#    new_psi_gauss = np.zeros((41, 41, 200))
#    for t in range(0, 200):
#        new_psi_gauss[:,:, t] = psi_gauss[:, t].reshape(41,41).transpose()
#    
#    x = np.linspace(0, 1, 41)
#    y = np.linspace(0, 1, 41)
#    
#    plt.style.use("ggplot")
#    fig = plt.figure(figsize = (9,7))
#    CS = plt.contourf(x, y, new_psi_gauss[:,:,0], 20, cmap = plt.cm.RdBu_r)
#    plt.colorbar(CS, orientation = "vertical")
#    plt.title(r'Contour field of $\psi(x, y, 0)$ in bounded domain', fontsize = 15)
#    plt.xlabel('x', fontsize = 13)
#    plt.ylabel('y', fontsize = 13)
#    plt.savefig('figs/gauss_boundary_2d.pdf', bbox_inches = 'tight')
#    
    
    
    
    
    
    
