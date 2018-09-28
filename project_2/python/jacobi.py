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
    if task == 1:
        h = 1/float(n)
        d *= 2/h**2
    elif task == 2:
        for i in range(n):
            h = rho_max/float(n)
            v = np.empty(n)
            v[i] = ((i+1)*h)**2
            d[i] = 2/h**2 + v[i]
    elif task == 3:
        for i in range(n):
            omega = np.array([0.01, 0.5, 1, 5])
            h = rho_max/float(n)
            v = np.empty(n)
            v[i] = omega[case-1]**2*(i*h)**2 + 1/(i+1*h)
            d[i] = 2/h**2 + v[i]
    a = -1/h**2
    A = np.zeros((n, n))
    for i in range(n):
        for j in range(n):
            if i == j:
                A[i, j] = d[i-1]
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
        eigenval[j-1] = d + 2*a*np.cos((j*np.pi)/(float(n)+1)) #exact
        
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
        n = len(A)
        
        if A[k,l] !=0.0:
            tau = (A[l,l] - A[k,k])/(2*A[k,l])
            if tau > 0:
                t = 1.0/(tau +m.sqrt(1.0 + tau**2))
            else:
                t = -1.0/(-tau +m.sqrt(1.0 + tau**2))
            c = 1/m.sqrt(1 + t**2)
            s = c*t
        else:
            c = 1.0
            s = 0.0
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
            print('advarsel')
        AMax, k, l = maxElem(A)
        iterations +=1
        if AMax < epsilon:
            
            t2 = time.time()
            print("Time used: {:.5f} sec".format(t2-t1))
            return diagonal(A),R, AMax, iterations
        rotate(A,R,k,l)
    

if __name__ == "__main__":
    A= make_matrix(10, 3, 1,3)
    eigenval = analytical_eig(A)
    eigval,eigvec, Amax , teller = jacobi(A)
