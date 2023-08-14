# -*- coding: utf-8 -*-
"""
Created on Tue Dec 14 13:57:19 2021

@author: Student
"""

import Crystal_Field_Calcuations as cef
import numpy as np
import matplotlib.pyplot as plt

###################################################################################
np.set_printoptions(suppress=True)

ion = 'Cr3'
L = 3
S = 3/2
Z = 2

Cr3 = cef.LS(L,S)

Splus = Cr3.Sp(L,S)
Sminus = Cr3.Sm(L,S)

Sz = np.diag(np.arange(-S,S+1,1))

Sp = np.roll(np.diag(np.diag(np.sqrt(S*(S+1) - Sz*(Sz+1)))),1)
Sm = np.roll(np.diag(np.diag(np.sqrt(S*(S+1) - Sz*(Sz-1)))),-1)

Sxyz = Sz + (1/2)*(Sp + Sm)

Sij = Sz*Sz + (1/2)*(np.dot(Sp,Sm) + np.dot(Sm,Sp))
SxyzSxyz = np.dot(Sxyz,Sxyz)

####################### Supporting Information for Air Stable and Layer 
####################### Dependent Ferromagnetism in Atomically Thin van der Waals CrPS4 
"""
H_spin = J_1{sum(S_i*S_j)} + J_2{sum(S_i*S_j)} + J_c{sum(S_i*S_j)} + K{sum(S_i^2)} - g*u_B*H*{sum(S_i)}
"""

J1 = -2.39 # meV
J2 = -2.39 # meV
#J3 = -0.51
Jc = 0.1530 # meV
K = -0.00417 # meV
mu_B = (5.7883818012e-5)*1000 #meV
g_e = 2.0023193
B = 0.1

H_spin = J1*(SxyzSxyz) +\
    J2*(SxyzSxyz) +\
    Jc*(SxyzSxyz) -\
    K*(Sz**2) +\
    mu_B*(g_e*Sz)*B


E_val,E_val_excitation,H_vt = Cr3.Diag(H_spin,printfunction=True)


Sz_m = np.zeros(4)
for k in range(4):
    Sz_m[k] = 2*np.dot(np.conjugate(H_vt[k]),np.dot(Sz,H_vt[k]))
print(Sz_m)

