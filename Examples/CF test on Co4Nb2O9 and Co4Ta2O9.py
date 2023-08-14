# -*- coding: utf-8 -*-
"""
Created on Fri Feb 12 16:03:00 2021

@author: chemuser
"""
import Crystal_Field_Calcuations as cef
import numpy as np
import matplotlib.pyplot as plt

##################   CNO    ########################



ion = 'Co2'
L = 3
S = 3/2
Z = 2
Co2_oct = cef.LS(L,S)

###################### Co1 ##############################


O1  =  np.array([1.642877,    3.056128,   12.947560])
O2  =  np.array([ -0.881323,    4.375178,   12.947560])
O3  =  np.array([ -0.761554,    1.529632,   12.947560])
O4  =  np.array([  1.098873,    4.480469,   10.609275])
O5  =  np.array([ -1.842836,    3.191886,   10.609275])
O6  =  np.array([  0.743964,    1.288583,   10.609275])
Co_oct = np.array([-0.000000,    2.986979,   11.415014])

O1_d = O1 - Co_oct
O2_d = O2 - Co_oct
O3_d = O3 - Co_oct
O4_d = O4 - Co_oct
O5_d = O5 - Co_oct
O6_d = O6 - Co_oct


d = np.array([O1_d,O2_d,O3_d,O4_d,O5_d,O6_d])
print('---------------- Ligand positions ---------------')
print(d)

B = Co2_oct.PC(ion,L,S,d,Z)

"""
rho = 0.13
Qeff = 0.7
B = Co2_oct.SOM(ion,L,S,d,Z,rho,Qeff)
"""

# For Co2+ Q is usually 0.6 - 1 while rho is usually 0.12 - 0.16

#################################### Crystal Field ##########################################
print('---------------- From literature ---------------')



O20 = Co2_oct.Olm(L,S,2,0)
#print('O02: ' + str(O02))
#print(O02_diag)

O40 = Co2_oct.Olm(L,S,4,0)
#print('O04: ' + str(O04))
#print(O04_diag)

O43 = Co2_oct.Olm(L,S,4,3)

B20 = B[0] #meV
B40 = B[3] #meV
B43 = B[6] #meV

#B20 = 0.1896 #meV
#B40 = 3.121 #meV
#B43 = -0.00183 #meV

print('B02 (meV): ' + str(B20))
print('B04 (meV): ' + str(B40))
print('B34 (meV): ' + str(B43))

print('---------------- Crystal Field cubic ---------------')

Hcf = B20*O20 + B40*O40 #+ B43*O43
Ecf_val,Ecf_val_excitation,H_cf_vt = Co2_oct.Diag(Hcf,printfunction=True)


print('---------------- Crystal Field non-cubic ---------------')

Hcf = B20*O20 + B40*O40 + B43*O43
Ecf_val,Ecf_val_excitation,H_cf_vt = Co2_oct.Diag(Hcf,printfunction=True)




#################################### Spin-Orbit coupling ##########################################
print('---------------- Crystal Field with Spin-Orbit coupling ---------------')


SO_matrix = Co2_oct.SO(ion,L,S)

Hcf_so = Hcf + SO_matrix

Ecf_so_val,Ecf_so_val_excitation,H_cf_so_vt = Co2_oct.Diag(Hcf_so,printfunction=True)

print('---------------- Crystal Field with Spin-Orbit coupling and Molecular Field ---------------')


Hg_x = Co2_oct.Molecular_Field_Sx(ion,L,S)

Hg_y = Co2_oct.Molecular_Field_Sy(ion,L,S)


Hcf_so_g = Hcf_so + Hg_x + Hg_y

Ecf_so_g_val,Ecf_so_g_val_excitation,H_cf_so_g_vt = Co2_oct.Diag(Hcf_so_g,printfunction=True)



#################################### Magnetic Moment ##########################################
print('--------- Magnetic Properties --------')




MagMoment_L_z,MagMoment_S_z,MagMoment_M_z = Co2_oct.MagMoment_z(L,S,Hcf_so)

MagMoment_L_x,MagMoment_S_x,MagMoment_M_x = Co2_oct.MagMoment_x(L,S,Hcf_so)

MagMoment_L_y,MagMoment_S_y,MagMoment_M_y = Co2_oct.MagMoment_y(L,S,Hcf_so)


################## specific heat ###############

C_schottky,T = Co2_oct.SpecificHeat(Ecf_val_excitation)
C_schottky,T = Co2_oct.SpecificHeat(Ecf_so_val_excitation)
C_schottky,T = Co2_oct.SpecificHeat(Ecf_so_g_val_excitation)


################## magnetic entropy ###############

entropy,T = Co2_oct.Entropy(Ecf_val_excitation)
entropy,T = Co2_oct.Entropy(Ecf_so_val_excitation)
entropy,T = Co2_oct.Entropy(Ecf_so_g_val_excitation)


############# Magnetic susceptibility ###############

X_z,T = Co2_oct.Susceptibility_z(L,S,Hcf_so,SO)
X_x,T = Co2_oct.Susceptibility_x(L,S,Hcf_so,SO)
X_y,T = Co2_oct.Susceptibility_y(L,S,Hcf_so,SO)

plt.plot(T,X_z,'b.',T,X_x,'r.',T,X_y,'y.')
plt.xlim(0,300)
#plt.ylim(0,2)
plt.xlabel('Temperature (K)')
plt.ylabel('Magnetic Susceptibility ($\mu_B$)')
plt.title('Magnetic Susceptibility')
#plt.legend('')
plt.grid(alpha=.4,linestyle='--')
plt.show()


