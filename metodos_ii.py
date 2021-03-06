# -*- coding: utf-8 -*-
"""METODOS_II.ipynb

Automatically generated by Colaboratory.

Original file is located at
    https://colab.research.google.com/drive/1-TbqRlBLshgezDLoRCOxRjIpE1YJv7pg

## Algoritmo para resolver un sistema $Ax = b$ por el método de iteración de Jacobi

#### Implementacion matricial (vectorizada)
"""



import numpy as np
from numpy import linalg as LA

def Jacobi(A,b,nmax,tol):
  n=len(b)
  # extract diagonal as a vector
  D=np.diag(A)
  # substract diagonal as a matrix (call "diag" again)
  LU = A-np.diag(D)
  x=np.zeros(n)
  k=0
  error=10
  while ( k < nmax and tol < error):
    k +=1
    xnew=(b-np.dot(LU,x))/D
    error=LA.norm(xnew-x)
    print("solucion x para la iteracion k= %2d con error %5.4f"%(k,error))
    print(" ",xnew)
    x=xnew
  if k==nmax:
    print("maximo numero de iteraciones fue alcanzado. ")
    print("probablemente no hay convegencia ")

A = np.array([[5.,-2.,3.],
              [-3.,8.,1.],
              [-3.,-1.,-6.]])

b = np.array([-1.,4.,0.])



nmax=1000
tol=0.0001

Jacobi(A,b,nmax,tol)

