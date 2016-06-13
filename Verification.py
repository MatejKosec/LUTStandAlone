# -*- coding: utf-8 -*-
"""
Created on Mon Jun 13 22:20:07 2016
A brief routine to create a set of random points in a rho, p plane and query
the LUT for interpolated properties through all SetTD_State functions. 

@author: matej kosec
"""

import scipy as sp
import matplotlib.pyplot as plt
import os 
import sys


#Create a set of random pressures and densities, the min/max range must be
#congruent with the range of the LUT table
rho_min = 59
pressure_min = 115000
rho_max = 100
pressure_max = 300500 
nrho = 20
np   = 20
rho = rho_max + sp.rand(nrho,np)*(rho_max-rho_min)
pressure = pressure_max + sp.rand(nrho,np)*(pressure_max-pressure_min)
plt.plot(rho,pressure)

os.system('Debug/TableReader>log')
with open('log') as l:
    l= l.readlines()
print l



"""
Test SetTDState_rhoe
"""


"""
Test SetTDState_Prho
"""

"""
Test SetTDState_PT
"""

"""
Test SetTDState_rhoT
"""

"""
Test SetTDState_Ps
"""

"""
Test 	SetTDState_hs
"""


