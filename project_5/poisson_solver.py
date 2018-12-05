import numpy as np


def poisson_jacobi_periodic(g, dx, dy, N_x, N_y, max_iter, f):
    dxdx = dx*dx
    dydy = dy*dy
    dxdxdydy = dxdx*dydy
    dxdx_plus_dydy_2 = 2*(dxdx*dydy)


    iter = 0
    diff = 1e20
    eps = 1e-6
    f_temp = np.empty(N_x*N_y)

    while (iter <= max_iter && abs(diff) > eps):
        diff = 0

        for i in range(0, N_x):
            for j in range(0,N_y):
                f_temp[i*N_y + j] = f[i*N_y + j]
        for i in range(0, N_x):
            for j in range(0,N_y):

                if i==0 & j ==0:
                    f[0*N_y+0] = (dydy*(f_temp[1*N_y +0] + f_temp[(N_x -2)*N_y + 0] )) + \
                    dxdx *(f_temp[(N_x -1)*N_y + 1] + f_temp[(N_x -1)*N_y + (N_y -2)]) -
                    dxdxdydy*g[(N_x-1)*N_y +0] / dxdx_plus_dydy_2
                elif i==0 && j == N_y-1:
                    f[0*N_y+(N_y -1)] = (dydy*(f_temp[1*N_y +(N_y -1)] + f_temp[(N_x -2)*N_y + (N_y -1)] )) + \
                    dxdx *(f_temp[(0*N_x +1)*N_y + 1] + f_temp[(0*N_y + (N_y -2)]) -
                    dxdxdydy*g[0*N_y +(N_y -1)]) / dxdx_plus_dydy_2
                elif i==N_x -1 && j == N_y-1:
                    f[(N_x -1)*N_y+(N_y -1)] = (dydy*(f_temp[1*N_y +(N_y -1)] + f_temp[(N_x -2)*N_y + (N_y -1)] )) + \
                    dxdx *(f_temp[(N_x -1)*N_y + 1] + f_temp[(N_x -1)*N_y + (N_y -2)]) -
                    dxdxdydy*g[(N_x-1)*N_y +(N_y -1)]) / dxdx_plus_dydy_2
                elif i==N_x -1 && j != N_y-1 && j != 0:
                    f[0*N_y+j] = (dydy*(f_temp[1*N_y + j] + f_temp[(N_x -2)*N_y + j)) + \
                    dxdx *(f_temp[(N_x -1)*N_y + 1] + f_temp[(N_x -1)*N_y + (N_y -2)]) -
                    dxdxdydy*g[(N_x-1)*N_y +(N_y -1)]) / dxdx_plus_dydy_2
                elif i !=N_x -1 && i != 0 && j == 0:
                    f[i*N_y + 0] = (dydy*(f_temp[(i+1)*N_y + 0] + f_temp[(N_x -2)*N_y + j)) + \
                    dxdx *(f_temp[(N_x -1)*N_y + 1] + f_temp[(N_x -1)*N_y + (N_y -2)]) -
                    dxdxdydy*g[(N_x-1)*N_y +(N_y -1)]) / dxdx_plus_dydy_2
