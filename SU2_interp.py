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
    Vinv = sp.linalg.inv(Vandermonde)
    
    print 'Vandermonde\n', Vandermonde
    print
    print 'Vandermonde inverse \n', Vinv
    Vinv = inverse(Vandermonde, 4)
    V22 = sp.copy(Vinv.T)
    print 'Transpose official'
    print V22
    for i in range(3):
        for j in range(i+1,4):
            d = Vinv[i,j]
            Vinv[i,j]= Vinv[j,i]
            Vinv[j,i]= d
    print 'Index ranspose\n', Vinv
    print 'Check transpose\n', Vinv-V22
    
    
    def SU2(x,y):
        RHS = sp.array([1,x,y,x*y])
        b = sp.dot(Vinv,RHS)
        return sp.dot(b,qz.T)
    SU2 = sp.vectorize(SU2)
    return SU2

def inverse(V,size):
    temp=sp.zeros([size,2*size])
    max_val = 0;
    
    
    for i in range(size):
        for j in range(size):
            temp[i,j]=V[i,j];
        temp[i,size+i]=1;
            
    #Pivot the rows 
    for i in range(size):
        max_val = temp[i,i];
        for j in range(i,size):
            if abs(temp[j,j]) > max_val:
                for k in range(size):
                    d = temp[i,k]
                    temp[i,k] = temp[j,k]
                    temp[j,k] = d
                                
    for k in range(size-1):
        if temp[k,k] != 0:
            for i in range(k+1,size):           
                c = temp[i,k]/temp[k,k]
                for j in range(2*size):
                    temp[i,j] = temp[i,j] - temp[k,j]*c
                    
                        
    for k in range(size-1,-1,-1):
        if temp[k,k] != 0:
            for i in range(k-1,0,-1):            
                    c = temp[i,k]/temp[k,k]
                    for j in range(2*size):
                        temp[i,j] = temp[i,j] - temp[k,j]*c             
            
    for i in range(size):
        c = temp[i,i]
        for j in range(size):
            temp[i,j+size] = temp[i,j+size]/c
    
    print 'Identity check'
    print sp.dot(temp[:,4:],V)

    return temp[:,4:]
    
    
    
    