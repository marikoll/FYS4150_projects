import numpy as np

def poisson_jacobi(g, bc_0y, bc_1y, bc_0x, bc_1x, dx, dy, N_x, N_y, max_iter, f):
    dxdx = dx**2
    dydy = dy**2
    dxdxdydy = dxdx*dydy
    dxdx_pluss_dydy_2 = 2*(dxdx +dydy)   
    
    for j in range(0, N_y):
        f[0 + j] = bc_0y[j]
        f[(N_x - 1)*N_y + j] = bc_1y[j]
    for i in range(0, N_x):
        f[i*N_y + 0] = bc_0x[i]
        f[i*N_y + (N_y-1)] = bc_1x[i]
    
    iterations = 0
    diff = 1E20
    eps = 1E-6
    f_temp = np.zeros(N_x*N_y)
    
    while iterations < max_iter and abs(diff) > eps:
        diff = 0.0
        
        for i in range(0, N_x):
            for j in range(0, N_y):
                f_temp[i*N_y + j] = f[i*N_y + j]
        
        for i in range(1, N_x -1):
            for j in range(1, N_y -1):
                f[i*N_y + j] = (dydy*(f_temp[(i+1)*N_y + j] + f_temp[(i-1)*N_y + j]) \
                 + dxdx*(f_temp[i*N_y + (j+1)] + f_temp[i*N_y + (j-1)]) - \
                 dxdxdydy*g[i*N_y + j])/dxdx_pluss_dydy_2
                diff += f[i*N_y + j] - f_temp[i*N_y + j]
        iterations += 1
    return f

    
def poisson_jacobi_periodic(g, dx, dy, N_x, N_y, max_iter, f):
    dxdx = dx*dx
    dydy = dy*dy
    dxdxdydy = dxdx*dydy
    dxdx_plus_dydy_2 = 2*(dxdx+dydy)


    iterations = 0
    diff = 1e20
    eps = 1e-6
    f_temp = np.empty(N_x*N_y)

    while (iterations <= max_iter and abs(diff) > eps):
        diff = 0.0

        for i in range(0, N_x):
            for j in range(0,N_y):
                f_temp[i*N_y + j] = f[i*N_y + j]
                
        for i in range(0, N_x):
            for j in range(0,N_y):
                # Corner, x = 0 and y = 0
                if i==0 and j ==0:
                    f[0*N_y+0] = (dydy*(f_temp[1*N_y +0] + f_temp[(N_x -2)*N_y + 0] ) + \
                    dxdx*(f_temp[0*N_y + 1] + f_temp[0*N_y + (N_y -2)]) - \
                    dxdxdydy*g[0*N_y +0]) / dxdx_plus_dydy_2
                     
                # Corner, x = 1 and y = 0     
                elif i==N_x-1 and j == 0:
                    f[(N_x-1)*N_y+0] = (dydy*(f_temp[1*N_y +0] + f_temp[(N_x -2)*N_y + (N_y +0)] ) + \
                    dxdx *(f_temp[(N_x -1)*N_y + 1] + f_temp[(N_x-1)*N_y + (N_y -2)]) - \
                    dxdxdydy*g[(N_x-1)*N_y +0] )/ dxdx_plus_dydy_2
                     
                # Corner, x = 0 and y = 1     
                elif i==0 and j == N_y-1:
                    f[0*N_y+(N_y -1)] = (dydy*(f_temp[1*N_y +(N_y -1)] + f_temp[(N_x -2)*N_y + (N_y -1)] ) + \
                    dxdx *(f_temp[0*N_y + 1] + f_temp[0*N_y + (N_y -2)]) - \
                    dxdxdydy*g[0*N_y +(N_y -1)] )/ dxdx_plus_dydy_2
                     
                # Corner, x = 1 and y = 1
                elif i==N_x-1 and j == N_y-1:
                    f[(N_x-1)*N_y+(N_y -1)] = (dydy*(f_temp[1*N_y +(N_y -1)] + f_temp[(N_x -2)*N_y + (N_y -1)] ) + \
                    dxdx *(f_temp[(N_x-1)*N_y + 1] + f_temp[(N_x-1)*N_y + (N_y -2)]) - \
                    dxdxdydy*g[(N_x-1)*N_y +(N_y -1)] )/ dxdx_plus_dydy_2
                       
                # Wall, x = 0    
                elif i==0 and j != 0 and j != N_y-1:
                    f[0*N_y+j] = (dydy*(f_temp[1*N_y + j] + f_temp[(N_x -2)*N_y + j]) + \
                    dxdx *(f_temp[0*N_y + j +1] + f_temp[0*N_y + (j -1)]) - \
                    dxdxdydy*g[0*N_y +j] )/ dxdx_plus_dydy_2
                     
                # Wall, x = 1    
                elif i ==N_x-1 and j != 0 and j != N_y -1:
                    f[(N_x-1)*N_y + j] = (dydy*(f_temp[1*N_y + j] + f_temp[(N_x -2)*N_y + j]) + \
                    dxdx *(f_temp[(N_x -1)*N_y +j + 1] + f_temp[(N_x-1)*N_y + j-1]) - \
                    dxdxdydy*g[(N_x-1)*N_y +j] )/ dxdx_plus_dydy_2
                
                # Wall, y = 0
                elif j == 0 and i != 0 and i != N_x-1:
                    f[i*N_y + 0] = (dydy*(f_temp[(i+1)*N_y + 0] + f_temp[(i-1)*N_y + 0]) + \
                    dxdx *(f_temp[i*N_y + 1] + f_temp[i*N_y + (N_y -2)]) - \
                    dxdxdydy*g[i*N_y +0]) / dxdx_plus_dydy_2
                     
                # Wall, y = 1
                elif j == N_y-1 and i != 0 and i != N_x-1:
                    f[i*N_y + (N_y-1)] = (dydy*(f_temp[(i+1)*N_y + (N_y-1)] + f_temp[(i-1)*N_y + (N_y-1)]) + \
                    dxdx *(f_temp[i*N_y + 1] + f_temp[i*N_y + (N_y -2)]) - \
                    dxdxdydy*g[i*N_y + (N_y-1)]) / dxdx_plus_dydy_2
                     
                else:
                    f[i*N_y +j] = (dydy*(f_temp[(i+1)*N_y + j] + f_temp[(i-1)*N_y + j]) + \
                    dxdx *(f_temp[i*N_y + (j+1)] + f_temp[i*N_y + (j-1)]) - \
                    dxdxdydy*g[i*N_y + j]) / dxdx_plus_dydy_2
                diff += f[i*N_y + j] - f_temp[i*N_y + j]
        iterations += 1
                
    return f
