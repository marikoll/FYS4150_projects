#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Unit test for the jacobi method applied on a tridiagonal symmetric Toepliz 
matrix. 

Functions: 
    test_jacobi_eigvec(n)
    test_jacobi_orthogonality(n)
"""

from jacobi import jacobi
from jacobi import make_matrix
from jacobi import analytical_eig
import numpy as np

def test_jacobi_eigvec(n):
    """
    Tests if the Jacobi method yields the same eigenvalues as numpy.linalg.eig
    and the analytical eigenvalues
    
    Input: 
        n - size of matrix
    """
    A = make_matrix(n,1)
    eigval = analytical_eig(A)
    eigval_np, teigvec_np = np.linalg.eig(A)
    eigval_J, eigvec_J, a, r, t= jacobi(A)
    np.testing.assert_allclose(sorted(eigval), sorted(eigval_np), rtol=1e-08, atol=0) #rtol - relative tolerance, atol - absolute tolerance
    np.testing.assert_allclose(sorted(eigval), sorted(eigval_J), rtol=1e-08, atol=0)
    print('test_jacobi_eigvec passed')
     

def test_jacobi_orthogonality(n):
    """
    Tests if the orthogonality is preserved for the eigenvectors calculated by 
    the Jacobi method
    
    Input: 
        n - size of matrix
    """
    A = make_matrix(n,1)
    eigval_np, teigvec_np = np.linalg.eig(A)
    eigval_J, eigvec_J, a, r, t = jacobi(A)
    err = 1E-5
    assert 1 - np.dot(eigvec_J[0], eigvec_J[0]) <= err
    assert np.dot(eigvec_J[0], eigvec_J[1]) <= err
    assert 1- np.dot(eigvec_J[3], eigvec_J[3]) <= err
    assert np.dot(eigvec_J[3], eigvec_J[1]) <= err
    
    print('test_jacobi_eigvec passed')
    
def test_largest_elem(n):
    """
    Tests if the largest element on the non-diagonal that was picked on the 
    last iteration in the Jacobi method was indeed the largest element on the
    rotated matrix
    
    Input: 
        n - size of matrix
    """
    A = make_matrix(n,1)
    eigval_J, eigvec_J,AMax,r, t = jacobi(A)
    maxelem = 0.0
    for i in range(n):
        for j in range(i+1,n):
            if abs(A[i,j]) >= maxelem:
                maxelem = abs(A[i,j])
    assert maxelem == AMax
    print('test_largest_elem passed')
    
    
if __name__ == "__main__":
    test_jacobi_eigvec(4)
    test_jacobi_orthogonality(4)
    test_largest_elem(4)
