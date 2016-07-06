# -*- coding: utf-8 -*-
"""
SU2 interpolation library test
"""
import scipy as sp;
from scipy import linalg

def interp2d(qx, qy, qz):
    Vandermonde = sp.zeros((4,4))
    Vandermonde[:,0] = 1
    Vandermonde[:,1] = qx
    Vandermonde[:,2] = qy
    Vandermonde[:,3] = qx*qy
    Vandermonde = Vandermonde.T
    Vinv = sp.linalg.inv(Vandermonde)
    def SU2(x,y):
        RHS = sp.array([1,x,y,x*y])
        b = sp.dot(Vinv,RHS)
        return sp.dot(b,qz.T)
    SU2 = sp.vectorize(SU2)
    return SU2

