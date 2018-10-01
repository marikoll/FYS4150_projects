#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Package containing functions for making a tridiogonal matrix (A) and calculate 
eigenvalues(of A) by an analytic equation aswell as eigenvalues and the
corresponding eigenvectors by Jacobi's method of rotation. 

Functions: 
    A = make_matrix(n, task, rho_max = None, case = None)            
    eigenval = analytical_eig(A)
    diagonal(A),R, AMax, teller = jacobi(A, epsilon = 1.0E-8)
    
"""
import numpy as np
from numpy import identity, diagonal
import math as m
import time
from numba import jit
import matplotlib.pyplot as plt


def make_matrix(n, task, rho_max = None, case = None):
    """
    Creates a tridiogongal (nxn) Toepliz - matrix
    
    input:
        n       - size of matrix
        task    - choose between task = 1 (task 2b), task = 2(task 2d),
                  task = 3 (task 2e)
        case    - choose between case = 1 (omega1), case = 2 (omega2),
                  case = 3 (omega3), case = 4 (omega4)
                      
                  Default = None
                      
                  USE IN TASK 2e
        rho_max - <int> maximum value of rho
            
                  Default = None
                      
                  USE IN TASK 2d AND 2e
    output: 
        A       - tridiagonal Toepliz matriz
    """
    d = np.ones(n)
    if task == 1:               # Make diagonal elements according to task 2b
        h = 1/float(n)
        d *= 2/h**2
    elif task == 2:
        for i in range(n):      # Make diagonal elements according to task 2d
            h = rho_max/float(n)
            v = np.empty(n)
            v[i] = ((i+1)*h)**2
            d[i] = 2/h**2 + v[i]
    elif task == 3:             # Make diagonal elements according to task 2e
        for i in range(n):
            omega = np.array([0.01, 0.5, 1, 5])
            h = rho_max/float(n)
            v = np.empty(n)
            v[i] = omega[case]**2*(i*h)**2 + 1/(i+1*h)
            d[i] = 2/h**2 + v[i]
    a = -1/h**2
    A = np.zeros((n, n))
    for i in range(n):
        for j in range(n):
            if i == j:
                A[i, j] = d[i]
        for j in range(n-1):
            if i == (j+1):
                A[j, i] = a
                A[i, j] = a
    return A

def analytical_eig(A):
    """
    Function for finding the analytical eigenvaliues for a tridiagonal Toepliz
    matric implemented on the buckling beam problem (task 2b)
    
    Input: 
        A       - tridiagonal Toepliz matriz
    Output: 
        eigenval- eigenvalues of A 
    
    """
    n = len(A)
    h = 1/float(n)
    d = 2/float(h)**2
    a = -1/float(h)**2
    eigenval = np.empty(n)
    for j in range(1,n+1):
        eigenval[j-1] = d + 2*a*np.cos((j*np.pi)/(float(n)+1)) # Analytic solution
        
    return eigenval


def jacobi(A, epsilon = 1.0E-8):
    """
    Implements Jacobi's method of rotation to find the eigenvalues and 
    eigenvectors of a tridiagonal, symmetric Toepliz matrix of size nxn
    
    Input: 
        A           - Tridiagonal Toepliz matriz
        epsilon     - Tolerance
                      Default = 1.0E-8
    Output: 
        diagonal(A) - Eigenvalues on the diagonal of the transformed matrix A
        R           - Eigenvectors
        AMax        - Maximum value of the last iteration (used for testing)
        iterations  - Number of iterations it took until all non-diagonal 
                      elements became zero
    """
    t1 = time.time()
    
    @jit(debug = True)
    def maxElem(A): 
        """
        Finds the largest elements on the non-diagonals of A 
        
        Input: 
            A           - Tridiagonal Toepliz matriz
        Output: 
            Amax        - Maximum value on non-diagonal
            k, l        - position of AMax
        """
        n = len(A)
        AMax = 0.0
        for i in range(n):
            for j in range(i+1,n):
                if abs(A[i,j]) >= AMax:
                    AMax = abs(A[i,j])
                    k = i;l = j
        return AMax, k, l

    @jit(debug = True)
    def rotate(A, R, k, l):
        """
        Rotates A  
        
        Input: 
            A           - Tridiagonal Toepliz matriz
            R           - nxn sized identity matrix  
            k, l        - position of largest elements on non-diagonal
        """
        n = len(A)
                
        
        tau = (A[l,l] - A[k,k])/(2*A[k,l])
        if tau > 0:
            t = 1.0/(tau +m.sqrt(1.0 + tau**2))
        else:
            t = -1.0/(-tau +m.sqrt(1.0 + tau**2))
        c = 1/m.sqrt(1 + t**2)
        s = c*t
        

        a_kk = A[k,k]
        a_ll = A[l,l]
        A[k,k] = c**2*a_kk - 2.0*c*s*A[k,l] + s**2*a_ll
        A[l,l] = s**2*a_kk + 2.0*c*s*A[k,l] + c**2*a_ll
        A[k,l] = 0.0
        A[l,k] = 0.0
        
        for i in range(n):
            if i != k and i != l:
                a_ik = A[i,k]
                a_il = A[i,l]
                A[i,k] = c*a_ik - s*a_il
                A[k,i] = A[i,k]
                A[i,l] = c*a_il + s*a_ik
                A[l,i] = A[i,l]
            r_ik = R[i,k]
            r_il = R[i,l]
            R[i,k] = c*r_ik - s*r_il
            R[i,l] = c*r_il + s*r_ik
    n = len(A)
    maxRot = n**3
    R = identity(n)*1.0
    iterations = 0
    for i in range(maxRot):
        if i == maxRot-1:
            print('Warning: max iterations reached')
        AMax, k, l = maxElem(A)
        iterations +=1
        if AMax < epsilon:
            
            t2 = time.time()
            time_tot = (t2-t1)
            return diagonal(A),R, AMax, iterations, time_tot 
        rotate(A,R,k,l)


if __name__ == "__main__":
    
    list1 = [10, 100, 200,250, 300, 350, 400]
    list2 = [5, 6, 7, 8]
#    outfile = open('result.txt', 'w')
#    for i in range(len(list1)):
#        for j in range(len(list2)):
#            A= make_matrix(list1[i],2,list2[j])
#            eigval,eigvec, Amax , teller, time_tot = jacobi(A)
#            eigval = sorted(eigval)
#            outfile.write('n: {}   Lambda_1: {:.5f}   Lambda_2: {:.5f}   Lambda_3: {:.5f}   Lambda_4: {:.5f}   CPU time: {:.5f} sec    iterations: {}\n'.format(list1[i], eigval[0], eigval[1], eigval[2], eigval[3], time_tot, teller))
#    outfile.close()
    outfile = open('result.txt', 'w')
    for i in range(len(list1)):
        A= make_matrix(list1[i],2,5)
        eigval,eigvec, Amax , teller, time_tot = jacobi(A)
        eigval = sorted(eigval)
        outfile.write('n: {}   Lambda_1: {:.5f}   Lambda_2: {:.5f}   Lambda_3: {:.5f}   Lambda_4: {:.5f}   CPU time: {:.5f} sec    iterations: {}\n'.format(list1[i], eigval[0], eigval[1], eigval[2], eigval[3], time_tot, teller))
    outfile.close()
    
#     A = make_matrix(10, 2, 9)
##     eigval_np, teigvec_np = np.linalg.eig(A)
#     eigval, eigvec, Amax, teller, time_tot = jacobi(A)
#     np.testing.assert_allclose(sorted(eigval), sorted(eigval_np), rtol=1e-08, atol=0)
#    eigenval = analytical_eig(A)
    
#    omegas = [0, 1, 2, 3]
#    n = 10
#    eigvals = np.zeros((len(omegas), len(range(n))))
#    plt.figure()
#    
#    for j in range(n):
#        A = make_matrix(10, 3, 9, omegas[1])
#        eigvals, eigvec, amax, t, time_tot = jacobi(A) 
#        plt.plot(eigvec[:,i], label = 'Omega = {}'.format(omegas[1]))
#    #plt.legend(fontsize = 10)
#    plt.show()
    
#    A = make_matrix(100, 3, 9, 0)
#    eigvals, eigvec, amax, t, time_tot = jacobi(A)
#    plt.plot(eigvec[:,19])
#    plt.show()
    
