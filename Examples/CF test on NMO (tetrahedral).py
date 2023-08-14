# -*- coding: utf-8 -*-
"""
Created on Sun Sep  6 22:01:34 2020

@author: brian
"""
import numpy as np
import crystalfield_test as cef


np.set_printoptions(suppress=True)

ion = 'Ni2'
L = 3
S = 1
Z = 2
Ni2p = cef.LS(L,S)


Ni = np.array([-2/3, 2/3, 0.05458])

O1 = np.array([-1.02182, 0.48909, -0.02633])

O2 = np.array([-0.48909, 1.02182, -0.02633])

O3 = np.array([-0.48909, 0.48909, -0.02633])

O4 = np.array([-2/3, 2/3, 0.2552])


O1_d = O1 - Ni
O2_d = O2 - Ni
O3_d = O3 - Ni
O4_d = O4 - Ni


d = np.array([O1_d,O2_d,O3_d,O4_d])
print('---------------- Ligand positions ---------------')
print(d)

B = Ni2p.PC(ion,L,S,d,Z)


#################################### Crystal Field ##########################################
print('---------------- What I am using ---------------')

O20 = Ni2p.Olm(L,S,2,0)
#print('O02: ' + str(O02))
#print(O02_diag)

O40 = Ni2p.Olm(L,S,4,0)
#print('O04: ' + str(O04))
#print(O04_diag)

O43 = Ni2p.Olm(L,S,4,3)



B02 = 0.6290083342348753
B04 = 1.3566225534052305
B34 = 4.9839923321314785

print('B02 (meV): ' + str(B02))
print('B04 (meV): ' + str(B04))
print('B34 (meV): ' + str(B34))

Hcf = B02*O20 + B04*O40 + B34*O43
Ecf = np.linalg.eigvals(Hcf)
Ecf1 = np.sort(Ecf)
Ecf_sol = Ecf1 - Ecf1[0]
print('Crystal Field: ' + str(Ecf_sol))


#################################### Spin-Orbit coupling ##########################################



SO_matrix = Ni2p.SO(ion,L,S)

Hcf_so = Hcf + SO_matrix
Ecf_so = np.linalg.eigvals(Hcf_so)
Ecf_so = np.sort(Ecf_so)
Ecf_so_sol = Ecf_so - Ecf_so[0]
print('Crystal Field with SO coupling: ' + str(Ecf_so_sol))

#################################### Spin-Spin interaction ##########################################

"""
p = 1.1/8.065548


Hss = Ni2p.SS(ion,L,S,p)
Hcf_so_ss = Hcf_so + Hss
"""

#################################### Molecular Field ##########################################


Hg = Ni2p.Molecular_Field_Sz(ion,L,S)

Hcf_so_g = Hcf_so + Hg

Ecf_so_g = np.linalg.eigh(Hcf_so_g)
Ecf_so_g_val = Ecf_so_g[0]
H_vt = Ecf_so_g[1].T

E_final = Ecf_so_g_val-Ecf_so_g_val[0]
print('Crystal Field with SO coupling and molecular field: ' + str(E_final))


#################################### Magnetic Moment ##########################################


print('\n')
print('--------- Magnetic Properties --------')




print('----------------- Lz ---------------')
Lz = Ni2p.Lz(L,S)
Lz_m = np.zeros((Ni2p.LS_degeneracy))
for k in range(Ni2p.LS_degeneracy):
    Lz_m[k] = np.dot(np.conjugate(H_vt[k]),np.dot(Lz,H_vt[k]))
print(Lz_m)
print('----------------- Sz ---------------')
Sz = Ni2p.Sz(L,S)
Sz_m = np.zeros((Ni2p.LS_degeneracy))
for k in range(Ni2p.LS_degeneracy):
    Sz_m[k] = np.dot(np.conjugate(H_vt[k]),np.dot(Sz,H_vt[k]))
print(Sz_m)
print('----------------- Mz ---------------')
Mz = Lz_m + 2*Sz_m
print(Mz)



print("----------- Crystal Field with SO coupling and magnetic field to induce nondegeneracy -----------")

mu_B = 5.7883818012e-2 # meV
Lx = Ni2p.Lx(L,S)
g_e = 2.0023193
Sx = Ni2p.Sx(L,S)
B = 0.1


Hcf_so_g_mag = Hcf_so_g + mu_B*(Lx + g_e*Sx)*B
Ecf_so_g_val,H_cef_SO_eigS = np.linalg.eigh(Hcf_so_g_mag)
E_final = Ecf_so_g_val - Ecf_so_g_val[0]
H_vt = H_cef_SO_eigS.T
print(E_final)


print('--------------- Lx ---------------')
Lp = Ni2p.Lp(L,S)
Lm = Ni2p.Lm(L,S)
Lp_m = np.zeros((Ni2p.LS_degeneracy))
for k in range(Ni2p.LS_degeneracy):
    Lp_m[k] = np.dot(np.conjugate(H_vt[k]),np.dot(Lp+Lm,H_vt[k]))
Lx_m = 0.5*Lp_m
print(Lx_m)



print('--------------- Sx ---------------')
Sp = Ni2p.Sp(L,S)
Sm = Ni2p.Sm(L,S)
Sp_m = np.zeros((Ni2p.LS_degeneracy))
for k in range(Ni2p.LS_degeneracy):
    Sp_m[k] = np.dot(np.conjugate(H_vt[k]),np.dot(Sp+Sm,H_vt[k]))
Sx_m = 0.5*Sp_m
print(Sx_m)

print('--------------- Mx ---------------')

Mx = Lx_m + g_e*Sx_m
print(Mx)


Lz = Ni2p.Lx(L,S)
Sz = Ni2p.Sx(L,S)

print("----------- Crystal Field with SO coupling and magnetic field to induce nondegeneracy -----------")

mu_B = 5.7883818012e-2 # meV
Ly = Ni2p.Ly(L,S)
g_e = 2.0023193
Sy = Ni2p.Sy(L,S)
B = 0.1

Hcf_so_g_mag = Hcf_so_g + mu_B*(Ly + g_e*Sy)*B
Ecf_so_g_val,H_cef_SO_eigS = np.linalg.eigh(Hcf_so_g_mag)
E_final = Ecf_so_g_val - Ecf_so_g_val[0]
H_vt = H_cef_SO_eigS.T
print(E_final)

print('--------------- Ly ---------------')
Ly = Ni2p.Ly(L,S)
Ly_m = np.zeros((Ni2p.LS_degeneracy))
for k in range(Ni2p.LS_degeneracy):
    Ly_m[k] = np.dot(np.conjugate(H_vt[k]),np.dot(Ly,H_vt[k]))
print(Ly_m)



print('--------------- Sy ---------------')
Sy = Ni2p.Sy(L,S)
Sy_m = np.zeros((Ni2p.LS_degeneracy))
for k in range(Ni2p.LS_degeneracy):
    Sy_m[k] = np.dot(np.conjugate(H_vt[k]),np.dot(Sy,H_vt[k]))
print(Sy_m)

print('--------------- My ---------------')

My = Ly_m + g_e*Sy_m
print(My)
