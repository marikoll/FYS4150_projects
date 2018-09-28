#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 27 17:45:16 2018

@author: maritkollstuen
"""

from jacobi import jacobi
from jacobi import make_matrix
from jacobi import analytical_eig
import numpy as np

def test_jacobi_eigvec(n):
    A = make_matrix(n,1)
    eigval = analytical_eig(A)
    eigval_np, teigvec_np = np.linalg.eig(A)
    eigval_J, eigvec_J,a,r = jacobi(A)
    np.testing.assert_allclose(sorted(eigval), sorted(eigval_np), rtol=1e-08, atol=0) #rtol - relative tolerance, atol - absolute tolerance
    np.testing.assert_allclose(sorted(eigval), sorted(eigval_J), rtol=1e-08, atol=0)
    print('test_jacobi_eigvec passed')
     

def test_jacobi_orthogonality(n):
    A = make_matrix(n,1)
    eigval_np, teigvec_np = np.linalg.eig(A)
    eigval_J, eigvec_J,a,r = jacobi(A)
    err = 1E-5
    assert 1 - np.dot(eigvec_J[0], eigvec_J[0]) <= err
    assert np.dot(eigvec_J[0], eigvec_J[1]) <= err
    assert 1- np.dot(eigvec_J[3], eigvec_J[3]) <= err
    assert np.dot(eigvec_J[3], eigvec_J[1]) <= err
    
    print('test_jacobi_eigvec passed')
    
if __name__ == "__main__":
    test_jacobi_eigvec(100)
    test_jacobi_orthogonality(100)
