"""
Simple(est) verification of the interpolation implementation and skweness test
"""
import scipy as sp
import matplotlib
from mpl_toolkits.mplot3d import axes3d
from scipy import interpolate
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter
from matplotlib import cm
from scipy import linalg
import SU2_interp
from shutil import copyfile
import os 
import sys

#Z = lambda x: sp.cos(x[0]-0.02) + sp.cos(x[1]+0.01) 
Z = lambda x: (x[0]) + (x[1])
InternalAngle = sp.deg2rad(80);
 
quad_x = sp.array([10.0, 11.0, 11+sp.cos(InternalAngle), 10+sp.cos(InternalAngle) ])
quad_y = sp.array([10.0, 10.0, 10+sp.sin(InternalAngle), 10+sp.sin(InternalAngle) ])
quad_z = Z([quad_x, quad_y])

x_samples = sp.zeros((100,100))
y_samples = sp.zeros((100,100))
for i in range(100):
    y_samples[:,i]= sp.linspace(10,10+sp.sin(InternalAngle),100)
    x_lim = 10+(y_samples[i,0]-10)/sp.tan(InternalAngle)
    x_samples[i,:] = sp.linspace(x_lim,1+x_lim,100)
    
#The correct values 
z_samples = Z([x_samples, y_samples])
#The SciPy values
interp_points_x = sp.copy(x_samples).reshape((100*100))
interp_points_y = sp.copy(y_samples).reshape((100*100))
interp_points   = sp.column_stack((interp_points_x,interp_points_y))
scipy_interp = sp.interpolate.griddata(sp.column_stack((quad_x,quad_y)),\
quad_z,interp_points,\
            method='linear') 
scipy_interp = scipy_interp.reshape((100,100))
scipy_error  = (scipy_interp - z_samples)/z_samples
#The SU2 values
SU2_Z  = SU2_interp.interp2d(quad_x,quad_y,quad_z)
su2_interp = SU2_Z(interp_points_x,interp_points_y)
su2_interp = su2_interp.reshape((100,100))
su2_error  = (su2_interp - z_samples)/z_samples

def plot(i,zz):
    plt.figure(i, figsize=(10,10))
    plt.plot(sp.hstack((quad_x,quad_x[0])),sp.hstack((quad_y,quad_y[0])), '-g')
    plt.plot(quad_x[0],quad_y[0], 'ro')
    plt.axis('equal')
    plt.grid('on')
    plt.xlim((9,12))
    plt.ylim((9,12))
    #plt.contourf(x_samples,y_samples,z_samples,100, interpolation=None)
    plt.contourf(x_samples,y_samples,abs(zz),100, interpolation=None)
    plt.colorbar()


plot(0,z_samples)
plot(1,scipy_error)
plot(2,su2_error)
plt.figure(3)
plt.plot(sp.sort(abs(scipy_error.reshape(100*100))), label='Scipy')
plt.plot(sp.sort(abs(su2_error.reshape(100*100))),label='SU2')
plt.legend()
Test2 = sp.matrix(
[[1, 0, 0, 0, 1, 0, 0, 0],\
 [1, 0.63290000000000646, -3941.25, -2494.4171250000254, 0, 1, 0, 0],\
 [1, 0, 1801.7300000000978, 0, 0, 0, 1, 0],\
 [1, 0.63290000000000646, -2168.0899999999674, -1372.1841609999933, 0, 0, 0, 1]])
print 'Test2'
print sp.around(Test2,2)
print 'Test2 inverse'
print sp.linalg.inv(Test2[:,:4])
print 'Elimination'
GT2I =  SU2_interp.inverse(Test2,4)
print 'Gauss inverse'
print GT2I
print 'Identity test'
print sp.dot(GT2I,Test2[:,:4])



            

