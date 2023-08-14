# -*- coding: utf-8 -*-
"""
Created on Wed Dec  8 18:12:29 2021

@author: Student
"""
import Crystal_Field_Calcuations as cef
import numpy as np
np.set_printoptions(suppress=True)

ion = 'Cu2'
L = 2
S = 1/2
Z = 2
Cu2_oct = cef.LS(L,S)

Cu =  np.array([14.688220,    1.615607,    5.571410])
N1 =   np.array([14.688220,    3.658199,    5.571410])
O1 =   np.array([14.121525,    1.623502,    3.708554])
O2 =   np.array([12.349230,    1.676360,    6.042306])
N2 =   np.array([14.688220,   -0.431790,    5.571410])
O3 =   np.array([15.254915,    1.623502,    7.434267])
O4 =   np.array([17.027211,    1.676360,    5.100515])



N1_d = N1 - Cu
N2_d = N2 - Cu
O1_d = O1 - Cu
O2_d = O2 - Cu
O3_d = O3 - Cu
O4_d = O4 - Cu

d = np.array([N1_d,N2_d,O1_d,O2_d,O3_d,O4_d])

B = Cu2_oct.PC(ion,L,S,d,Z)

rho = 0.2
Qeff = .2
B = Cu2_oct.SOM(ion,L,S,d,Z,rho,Qeff)

B = B/1000


O20 = Cu2_oct.Olm(L,S,2,0)
O21 = Cu2_oct.Olm(L,S,2,1)
O22 = Cu2_oct.Olm(L,S,2,2)
O40 = Cu2_oct.Olm(L,S,4,0)
O41 = Cu2_oct.Olm(L,S,4,1)
O42 = Cu2_oct.Olm(L,S,4,2)
O43 = Cu2_oct.Olm(L,S,4,3)
O44 = Cu2_oct.Olm(L,S,4,4)


Hcf = B[0]*O20 + B[1]*O21 + B[2]*O22 + B[3]*O40 + B[4]*O41 + B[5]*O42 + B[6]*O43 + B[7]*O44
Ecf_val,Ecf_val_excitation,H_cf_vt = Cu2_oct.Diag(Hcf)