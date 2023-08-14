# -*- coding: utf-8 -*-
"""
Created on Sun Sep  6 22:01:34 2020

@author: brian
"""
import numpy as np
import Crystal_Field_Calcuations as cef
import matplotlib.pyplot as plt

np.set_printoptions(suppress=True)

ion = 'Fe2'
L = 2
S = 2
Z = 2
Fe2p = cef.LS(L,S)


Fe = np.array([-1/3, 1/3, 0.01301])

O1 = np.array([-0.48830, 0.48830, -0.13710])

O2 = np.array([-0.48830, 0.02340, -0.13710])

O3 = np.array([-0.02340, 0.48830, -0.13710])

O4 = np.array([-1/6, 1/6, 0.13440])

O5 = np.array([-0.66580, 1/6, 0.13440])

O6 = np.array([-1/6, 0.66580, 0.13440])

O1_d = O1 - Fe
O2_d = O2 - Fe
O3_d = O3 - Fe
O4_d = O4 - Fe
O5_d = O5 - Fe
O6_d = O6 - Fe


d = np.array([O1_d,O2_d,O3_d,O4_d,O5_d,O6_d])
print('---------------- Ligand positions ---------------')
print(d)

B = Fe2p.PC(ion,L,S,d,Z)

#################################### Crystal Field ##########################################

print('---------------- From literature ---------------')
O20 = Fe2p.Olm(L,S,2,0)
#print('O02: ' + str(O02))
#print(O02_diag)

O40 = Fe2p.Olm(L,S,4,0)
#print('O04: ' + str(O04))
#print(O04_diag)

O43 = Fe2p.Olm(L,S,4,3)


Dq = 685 # meV
A2 = -596.0 # meV
A4 = -332.0 # meV

B02 = (A2/3)/8.065548
B04 = (A4/12 + Dq/18)/8.065548
B34 = (-20*Dq/(9*np.sqrt(2)))/8.065548

print('B02 (meV): ' + str(B02))
print('B04 (meV): ' + str(B04))
print('B34 (meV): ' + str(B34))

Hcf = B02*O20 + B04*O40 + B34*O43
Ecf = np.linalg.eigvals(Hcf)
Ecf1 = np.sort(Ecf)
Ecf_sol = Ecf1 - Ecf1[0]
print('Crystal Field: ' + str(Ecf_sol))


#################################### Spin-Orbit coupling ##########################################


SO = -116.5/8.065548

SO_matrix = Fe2p.SO(ion,L,S,SO)

Hcf_so = Hcf + SO_matrix
Ecf_so = np.linalg.eigvals(Hcf_so)
Ecf_so = np.sort(Ecf_so)
Ecf_so_sol = Ecf_so - Ecf_so[0]
print('Crystal Field with SO coupling: ' + str(Ecf_so_sol))

#################################### Spin-Spin interaction ##########################################


p = 1.1


Hss = Fe2p.SS(ion,L,S,p)
Hcf_so_ss = Hcf_so + Hss


#################################### Molecular Field ##########################################

Hm = 13.6/8.065548

Hg = Fe2p.Molecular_Field_Sz(ion,L,S,Hm)

Hcf_so_g = Hcf_so + Hg

Ecf_so_g = np.linalg.eigh(Hcf_so_g)
Ecf_so_g_val = Ecf_so_g[0]
H_vt = Ecf_so_g[1].T

E_final = Ecf_so_g_val-Ecf_so_g_val[0]
print('Crystal Field with SO coupling and molecular field: ' + str(E_final))


#################################### Magnetic Moment ##########################################
print('--------- Magnetic Properties --------')




print('----------------- Lz ---------------')
Lz = Fe2p.Lz(L,S)
Lz_m = np.zeros((Fe2p.LS_degeneracy))
for k in range(Fe2p.LS_degeneracy):
    Lz_m[k] = np.dot(np.conjugate(H_vt[k]),np.dot(Lz,H_vt[k]))
print(Lz_m)
print('----------------- Sz ---------------')
Sz = Fe2p.Sz(L,S)
Sz_m = np.zeros((Fe2p.LS_degeneracy))
for k in range(Fe2p.LS_degeneracy):
    Sz_m[k] = np.dot(np.conjugate(H_vt[k]),np.dot(Sz,H_vt[k]))
print(Sz_m)
print('----------------- Mz ---------------')
Mz = Lz_m + 2*Sz_m
print(Mz)



print("----------- Crystal Field with SO coupling and magnetic field to induce nondegeneracy -----------")

mu_B = 5.7883818012e-2 # meV
Lx = Fe2p.Lx(L,S)
g_e = 2.0023193
Sx = Fe2p.Sx(L,S)
B = 0.1


Hcf_so_ss_g_mag = Hcf_so_g + mu_B*(Lx + g_e*Sx)*B
Ecf_so_ss_g_val,H_cef_SO_eigS = np.linalg.eigh(Hcf_so_ss_g_mag)
E_final = Ecf_so_ss_g_val - Ecf_so_ss_g_val[0]
H_vt = H_cef_SO_eigS.T
print(E_final)


print('--------------- Lx ---------------')
Lp = Fe2p.Lp(L,S)
Lm = Fe2p.Lm(L,S)
Lp_m = np.zeros((Fe2p.LS_degeneracy))
for k in range(Fe2p.LS_degeneracy):
    Lp_m[k] = np.dot(np.conjugate(H_vt[k]),np.dot(Lp+Lm,H_vt[k]))
Lx_m = 0.5*Lp_m
print(Lx_m)



print('--------------- Sx ---------------')
Sp = Fe2p.Sp(L,S)
Sm = Fe2p.Sm(L,S)
Sp_m = np.zeros((Fe2p.LS_degeneracy))
for k in range(Fe2p.LS_degeneracy):
    Sp_m[k] = np.dot(np.conjugate(H_vt[k]),np.dot(Sp+Sm,H_vt[k]))
Sx_m = 0.5*Sp_m
print(Sx_m)

print('--------------- Mx ---------------')

Mx = Lx_m + g_e*Sx_m
print(Mx)


Lz = Fe2p.Lx(L,S)
Sz = Fe2p.Sx(L,S)

print("----------- Crystal Field with SO coupling and magnetic field to induce nondegeneracy -----------")

mu_B = 5.7883818012e-2 # meV
Ly = Fe2p.Ly(L,S)
g_e = 2.0023193
Sy = Fe2p.Sy(L,S)
B = 0.1

Hcf_so_ss_g_mag = Hcf_so_g + mu_B*(Ly + g_e*Sy)*B
Ecf_so_ss_g_val,H_cef_SO_eigS = np.linalg.eigh(Hcf_so_ss_g_mag)
E_final = Ecf_so_ss_g_val - Ecf_so_ss_g_val[0]
H_vt = H_cef_SO_eigS.T
print(E_final)

print('--------------- Ly ---------------')
Ly = Fe2p.Ly(L,S)
Ly_m = np.zeros((Fe2p.LS_degeneracy))
for k in range(Fe2p.LS_degeneracy):
    Ly_m[k] = np.dot(np.conjugate(H_vt[k]),np.dot(Ly,H_vt[k]))
print(Ly_m)



print('--------------- Sy ---------------')
Sy = Fe2p.Sy(L,S)
Sy_m = np.zeros((Fe2p.LS_degeneracy))
for k in range(Fe2p.LS_degeneracy):
    Sy_m[k] = np.dot(np.conjugate(H_vt[k]),np.dot(Sy,H_vt[k]))
print(Sy_m)

print('--------------- My ---------------')

My = Ly_m + g_e*Sy_m
print(My)


################## specific heat ###############

C_schottky,T = Fe2p.SpecificHeat(Ecf_sol)
C_schottky,T = Fe2p.SpecificHeat(Ecf_so_sol)
C_schottky,T = Fe2p.SpecificHeat(E_final)


################## magnetic entropy ###############

entropy,T = Fe2p.Entropy(Ecf_sol)
entropy,T = Fe2p.Entropy(Ecf_so_sol)
entropy,T = Fe2p.Entropy(E_final)


############# Magnetisation ###############
Fields = 5 # T


Zeeman = Fe2p.Zeeman_z(L,S,Hcf_so_g,Fields)
Zeeman = Fe2p.Zeeman_x(L,S,Hcf_so_g,Fields)
Zeeman = Fe2p.Zeeman_y(L,S,Hcf_so_g,Fields)

############# Magnetic susceptibility ###############

X_z,T = Fe2p.Susceptibility_z(L,S,Hcf_so,SO)
X_x,T = Fe2p.Susceptibility_x(L,S,Hcf_so,SO)
X_y,T = Fe2p.Susceptibility_y(L,S,Hcf_so,SO)

plt.plot(T,X_z,'bo',T,X_x,'ro',T,X_y,'yo')
#plt.ylim(0,0.2)
plt.xlabel('Temperature (K)')
plt.ylabel('Magnetic Susceptibility ($\mu_B$)')
plt.title('Magnetic Susceptibility')
#plt.legend('')
plt.grid(alpha=.4,linestyle='--')
plt.show()

plt.plot(T,1/X_z,'bo',T,1/X_x,'ro',T,1/X_y,'yo')
#plt.ylim(0,250)
plt.xlabel('Temperature (K)')
plt.ylabel('Inverse Magnetic Susceptibility ($\mu_B$)')
plt.title('Inverse Magnetic Susceptibility')
#plt.legend('')
plt.grid(alpha=.4,linestyle='--')
plt.show()



