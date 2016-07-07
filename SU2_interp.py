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
    print 'Vandermonde inverse official \n', Vinv    
    Vinv = inverse(Vandermonde, 4)
    print 'Vandermonde inverse Gauss \n', Vinv
    V22 = sp.copy(Vinv.T)
    print 'Identity check'
    print sp.dot(Vinv,Vandermonde)
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
    for k in range(size-1):
        max_idx = k
        max_val = temp[k,k];
        for j in range(k,size):
            if abs(temp[j,k]) > max_val:
                max_idx = j;
                max_val = temp[j,k];
        for j in range(2*size):
                d = temp[k,j]
                temp[k,j] = temp[max_idx,j]
                temp[max_idx,j] = d                                
        for i in range(k+1,size):           
            c = temp[i,k]/temp[k,k]
            for j in range(0,2*size):
                temp[i,j] = temp[i,j] - temp[k,j]*c
    
    print "Reduced Echelon"                    
    print temp[:,:4]
    print "Reduced Echelon X"                    
    print temp[:,4:]
                        
    for k in range(size-1,0,-1):
        if temp[k,k] != 0:
            for i in range(k-1,-1,-1):            
                    c = temp[i,k]/temp[k,k]
                    for j in range(0,2*size):
                        temp[i,j] = temp[i,j] - temp[k,j]*c             
    
    print "Reduced Reduced Echelon"                    
    print temp[:,:4]
    print "Reduced Reduced Echelon X"                    
    print temp[:,4:]
        
        
    for i in range(size):
        c = temp[i,i]
        for j in range(size):
            temp[i,j+size] = temp[i,j+size]/c
    
    

    return temp[:,4:]
    
    
    
    