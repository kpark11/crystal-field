# -*- coding: utf-8 -*-
"""
Created on Sun Sep  6 22:01:34 2020

@author: brian
"""
import numpy as np
import Crystal_Field_Calcuations as cef


np.set_printoptions(suppress=True)

ion = 'Ni2'
L = 3
S = 1
Z = 2
Ni2p = cef.LS(L,S)


Ni_pos_oct = np.array([1/3, 2/3, 0.6172])
O_pos1_oct = np.array([0.02182, 0.51091, 0.47367])
O_pos2_oct = np.array([0.48909, 0.97818, 0.47367])
O_pos3_oct = np.array([0.48909, 0.51091, 0.47367])
O_pos4_oct = np.array([0.16834, 0.33668, 0.7422])
O_pos5_oct = np.array([0.66332, 0.83166, 0.7422])
O_pos6_oct = np.array([0.16834, 0.83166, 0.7422])

O1_oct = O_pos1_oct - Ni_pos_oct
O2_oct = O_pos2_oct - Ni_pos_oct
O3_oct = O_pos3_oct - Ni_pos_oct
O4_oct = O_pos4_oct - Ni_pos_oct
O5_oct = O_pos5_oct - Ni_pos_oct
O6_oct = O_pos6_oct - Ni_pos_oct

d = np.array([O1_oct,O2_oct,O3_oct,O4_oct,O5_oct,O6_oct])
print('---------------- Ligand positions ---------------')
print(d)

B = Ni2p.PC(ion,L,S,d,Z)


#################################### Crystal Field ##########################################
print('---------------- Crystal Field ---------------')

O20 = Ni2p.Olm(L,S,2,0)
#print('O02: ' + str(O02))
#print(O02_diag)

O40 = Ni2p.Olm(L,S,4,0)
#print('O04: ' + str(O04))
#print(O04_diag)

O43 = Ni2p.Olm(L,S,4,3)


Dq = 806.5548 # cm-1
A2 = -33.3 # cm-1
A4 = -399.4 # cm-1

B02 = (A2/3)/8.065548
B04 = (A4/12 + Dq/18)/8.065548
B34 = (-20*Dq/(9*np.sqrt(2)))/8.065548

print('B02 (meV): ' + str(B02))
print('B04 (meV): ' + str(B04))
print('B34 (meV): ' + str(B34))

Hcf = B02*O20 + B04*O40 + B34*O43
Ecf_val,Ecf_val_excitation,H_cf_vt = Ni2p.Diag(Hcf,printfunction=True)



#################################### Spin-Orbit coupling ##########################################

SO_matrix = Ni2p.SO(ion,L,S)

Hcf_so = Hcf + SO_matrix
Ecf_so_val,Ecf_so_val_excitation,H_cf_so_vt = Ni2p.Diag(Hcf_so,printfunction=True)

#################################### Spin-Spin interaction ##########################################

"""
p = 1.1/8.065548


Hss = Ni2p.SS(ion,L,S,p)
Hcf_so_ss = Hcf_so + Hss
"""

#################################### Molecular Field ##########################################


Hg = Ni2p.Molecular_Field_Sz(ion,L,S)

Hcf_so_g = Hcf_so + Hg

Ecf_so_g_val,Ecf_so_g_val_excitation,H_cf_so_g_vt = Ni2p.Diag(Hcf_so_g,printfunction=True)


#################################### Magnetic Moment ##########################################


print('\n')
print('--------- Magnetic Properties --------')




MagMoment_L_z,MagMoment_S_z,MagMoment_M_z = Ni2p.MagMoment_z(L,S,Hcf_so)

MagMoment_L_x,MagMoment_S_x,MagMoment_M_x = Ni2p.MagMoment_x(L,S,Hcf_so)

MagMoment_L_y,MagMoment_S_y,MagMoment_M_y = Ni2p.MagMoment_y(L,S,Hcf_so)



"""

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


"""
