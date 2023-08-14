# -*- coding: utf-8 -*-
"""
Created on Fri Feb 12 16:03:00 2021

@author: chemuser
"""
import crystalfield_test as cef
import numpy as np
import matplotlib.pyplot as plt

np.set_printoptions(suppress=True)


ion = 'Co2'
L = 1
S = 3/2
Z = 2
Co2_tet = cef.LS(L,S)




#################################### Crystal Field ##########################################
print('---------------- From literature ---------------')



O20 = Co2_tet.Olm(L,S,2,0)
#print('O02: ' + str(O02))
#print(O02_diag)

O40 = Co2_tet.Olm(L,S,4,0)
#print('O04: ' + str(O04))
#print(O04_diag)

O43 = Co2_tet.Olm(L,S,4,3)



B20 = 0.1896 #meV
B40 = 3.121 #meV
B43 = -0.00183 #meV

print('B02 (meV): ' + str(B20))
print('B04 (meV): ' + str(B40))
print('B34 (meV): ' + str(B43))

print('---------------- Crystal Field ---------------')

Hcf = B20*O20 + B40*O40 + B43*O43
Ecf_val,Ecf_val_excitation,H_cf_vt = Co2_tet.Diag(Hcf)




#################################### Spin-Orbit coupling ##########################################
print('---------------- Crystal Field with Spin-Orbit coupling ---------------')


SO = 19 #meV

SO_matrix = Co2_tet.SO(ion,L,S,SO)

Hcf_so = Hcf + SO_matrix

Ecf_so_val,Ecf_so_val_excitation,H_cf_so_vt = Co2_tet.Diag(Hcf_so)

print('---------------- Crystal Field with Spin-Orbit coupling and Molecular Field ---------------')

Hm = 2  #meV

Hg = Co2_tet.Molecular_Field_Sz(ion,L,S,Hm)

Hcf_so_g = Hcf_so + Hg

Ecf_so_g_val,Ecf_so_g_val_excitation,H_cf_so_g_vt = Co2_tet.Diag(Hcf_so_g)