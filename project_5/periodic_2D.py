
import numpy as np
import matplotlib.pyplot as plt
from poisson_solver import poisson_jacobi_periodic


def periodic_matrix(n_rows, n_cols):
    """
    Function to make periodic matrix to solve the 1D PDE in a periodic domain.
    The function returns a n_rows times n_cols matrix which is used in a matrix solver.

    Input:
        n_rows      <int>       number of rows in matrix
        n_cols      <int>       number of columns in matrix

    Output:
        A           <float64>   matrix with size a n_rows times n_cols
    """
    A = np.zeros((n_rows, n_cols))
    for i in range(0,n_rows):
        for j in range(0,n_cols):
            if i == j:
                A[i, j] = -2.0
            elif abs(i-j) ==1:
                A[i,j] = 1.0
    A[0, n_cols-1] = 1.0
    A[n_rows -1, 0] = 1.0

    return A

def initialize(N_x, N_y, dx, dy, case):
    """
    Function to set the inital values of the sine and the gaussian wave.
    This will be used in a function below to calculate the full solution.

    Input:
        N_x         <int>       number of points in x-direction
        N_y         <int>       number of points in y-direction
        dx          <float>     stepsize in x-direction
        dy          <float>     stepsize in y-direction
        case        <str>       select which case to initialize, sine or gaussian

    Output
        init_psi    <float64>   array with initial values for psi
        init_zeta   <float64>   array with inital values for zeta

    """
    init_psi = np.zeros(N_x*N_y)
    init_zeta = np.zeros(N_x*N_y)

    if case == 'sine':
        for i in range(0, N_x):
            for j in range(0, N_y):
                x = i*dx
#                y = i*dy
                init_psi[i*N_y + j] = np.sin(4.0*np.pi*x)
                init_zeta[i*N_y + j] = -16.0*np.pi**2*np.sin(4.0*np.pi*x)
    if case == 'gauss':
        for i in range(0, N_x):
            for j in range(0, N_y):
                x = i*dx
                sigma = 0.1
                init_psi[i] = np.exp(-((x-0.5)/sigma)**2)
                init_zeta[i] = 4.0*((x-0.5)/sigma**2)**2 - (2/sigma**2)*np.exp(-((x-0.5)/sigma)**2)

    return init_psi, init_zeta



def leapfrog(N_x,N_y, dx,dy, T, dt, case):
    """
    Function to use for time stepping using leapfrog (centered difference) to solve the PDE
    in question.
    Uses functions initialize in order to solve the PDE in question.
    poisson_jacobi_periodic is used to solve the PDE wrt x (Ax = y)


    Input:
        N_x     <int>       number of points in x-direction
        N_y     <int>       number of points in x-direction
        dx      <float>     stepsize in x-direction
        dy      <float>     stepsize in y-direction
        T       <int>       time limit
        dt      <float>     time stepsize
        case    <str>       select which case to initialize, sine or gaussian

    Output:
        out_data <float64>  array with time in first row and psi (solution)

    """
    psi_0, zeta_0 = initialize(N_x, N_y, dx, dy, case)

    alpha = dt/(2*dx)
    gamma =  dt/dx

    psi_prev = np.zeros(N_x*N_y)
    psi_curr = np.zeros(N_x*N_y)
    zeta_prev = np.zeros(N_x*N_y)
    zeta_pp = np.zeros(N_x*N_y)
    zeta_curr = np.zeros(N_x*N_y)

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

    data_out = np.zeros((N_x*N_y+1, int(float(T)/dt)+1))
    t = 0.0
    data_out[0,0] = t

    data_out[1:, 0] = psi_0[:]

    n = 0
    n2 = 1

    while t < T:

        for i in range(1, N_x-2):
            for j in range(0,N_y-1):
                zeta_curr[i*N_y + j] = zeta_pp[i*N_y + j] \
                - gamma*(psi_prev[(i+1)*N_y + j] - psi_prev[(i-1)*N_y + j])
        for j in range(0, N_y):
            zeta_curr[0*N_y + j] = zeta_pp[0*N_y + j] - \
            gamma*(psi_prev[1*N_y + j] - psi_prev[(N_x -2)*N_y + j])
            zeta_curr[(N_x-1)*N_y + j] = zeta_curr[0*N_y + j]

        psi_curr = poisson_jacobi_periodic(zeta_curr, dx, dt, N_x, N_y, 50, psi_curr)

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
        """
        Calling function with variables as set below.
        Plotting Hovmuller diagrams in 2D
        """

    T = 100
    dt = 0.01

    dx = 1.0/40
    dy = 1.0/40
    L = 1.0
    N_x = int(L/dx + 1)
    N_y = int(L/dy +1)

    data_out_sine = leapfrog(N_x, N_y, dy, dx, T, dt, 'sine')

    psi_sine = data_out_sine[1:, :200]
    new_psi_sine = np.zeros((41, 41, 200))
    for t in range(0, 200):
        new_psi_sine[:,:, t] = psi_sine[:, t].reshape(41,41).transpose()

    x = np.linspace(0, 1, 41)
    y = np.linspace(0, 1, 41)

    plt.figure(1, figsize = (10, 8))
    plt.subplot(2, 2, 1)
    plt.style.use("ggplot")
    CS = plt.contourf(x, y, new_psi_sine[:,:,0], 20, cmap = plt.cm.coolwarm)
    plt.colorbar(CS, orientation = "vertical")
    plt.title(r'$\psi(x, y, {:.2f})$'.format(psi_sine[0,0]), fontsize = 15)
    plt.ylabel('y', fontsize = 13)
    plt.subplot(2, 2, 2)
    plt.style.use("ggplot")
    CS = plt.contourf(x, y, new_psi_sine[:,:,2], 20, cmap = plt.cm.coolwarm)
    plt.colorbar(CS, orientation = "vertical")
    plt.title(r'$\psi(x, y, {:.2f})$'.format(psi_sine[0,2]), fontsize = 15)
    plt.subplot(2, 2, 3)
    plt.style.use("ggplot")
    CS = plt.contourf(x, y, new_psi_sine[:,:,4], 20, cmap = plt.cm.coolwarm)
    plt.colorbar(CS, orientation = "vertical")
    plt.title(r'$\psi(x, y, {:.2f})$'.format(psi_sine[0,4]), fontsize = 15)
    plt.xlabel('x', fontsize = 13)
    plt.ylabel('y', fontsize = 13)
    plt.subplot(2, 2, 4)
    plt.style.use("ggplot")
    CS = plt.contourf(x, y, new_psi_sine[:,:,6], 20, cmap = plt.cm.coolwarm)
    plt.colorbar(CS, orientation = "vertical")
    plt.title(r'$\psi(x, y, {:.2f})$'.format(psi_sine[0,6]), fontsize = 15)
    plt.xlabel('x', fontsize = 13)
    plt.savefig('figs/sine_periodic_2d.pdf', bbox_inches = 'tight')
