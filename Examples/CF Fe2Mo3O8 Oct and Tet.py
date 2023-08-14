# -*- coding: utf-8 -*-
"""
Created on Sun Sep  6 22:01:34 2020

@author: brian
"""
import numpy as np
import Crystal_Field_Calculations as cef
import matplotlib.pyplot as plt


np.set_printoptions(suppress=True)  

print("Fe2Mo3O8 Point Charge Modeling:")
a = 5.7732 # Angstrom
b = a # Angstrom
c = 10.0542 # Angstrom


ion = 'Fe2'
L = 2
S = 2
Z = 2
Fe2_tet = cef.LS(L,S)


Fe_tet = np.array([0.000000,    3.333159,    9.570593])
O1 = np.array([1.544620,    2.441372,    8.675769])
O2 = np.array([-1.544620,    2.441372,    8.675769])
O3 = np.array([0.000000,    3.333159,   11.532168])
O4 = np.array([0.000000,    5.116732,    8.675769])


O1_d = O1 - Fe_tet
O2_d = O2 - Fe_tet
O3_d = O3 - Fe_tet
O4_d = O4 - Fe_tet


d = np.array([O1_d,O2_d,O3_d,O4_d])
B = Fe2_tet.PC(ion,L,S,d,Z)


#rho = 0.1
#Qeff = 0.7
#B = Fe2_tet.SOM(ion,L,S,d,Z,rho,Qeff)

B02 = B[0]
#B02_oct = B02
B04 = B[3]
#B04_oct = B04
B34 = B[6]
#B34_oct = B34


#################################### Crystal Field ##########################################
print('---------------- From literature ---------------')





Dq = 302.8 # cm-1
A2 = -33.3 # cm-1
A4 = -399.4 # cm-1

B02 = (A2/3)
B04 = (A4/12 + Dq/18)
B34 = (-20*Dq/(9*np.sqrt(2)))





O20 = Fe2_tet.Olm(L,S,2,0)
#print('O02: ' + str(O02))
#print(O02_diag)

O40 = Fe2_tet.Olm(L,S,4,0)
#print('O04: ' + str(O04))
#print(O04_diag)

O43 = Fe2_tet.Olm(L,S,4,3)


print('B02 (meV): ' + str(B02))
print('B04 (meV): ' + str(B04))
print('B34 (meV): ' + str(B34))



print('---------------- Crystal Field ---------------')

Hcf_tet = B02*O20 + B04*O40 + B34*O43
Ecf_val_tet,Ecf_val_excitation_tet,H_cf_vt_tet = Fe2_tet.Diag(Hcf_tet,printfunction=True)



#################################### Spin-Orbit coupling ##########################################
print('---------------- Crystal Field with Spin-Orbit coupling ---------------')


SO = -116.5

SO_matrix = Fe2_tet.SO(ion,L,S,SO)

Hcf_so_tet = Hcf_tet + SO_matrix

Ecf_so_val_tet,Ecf_so_val_excitation_tet,H_cf_so_vt_tet = Fe2_tet.Diag(Hcf_so_tet,printfunction=True)

#################################### Spin-Spin interaction ##########################################
print('---------------- Crystal Field with Spin-Orbit coupling, spin-spin interaction, and Molecular Field ---------------')


p = 1.1


Hss_tet = Fe2_tet.SS(ion,L,S,p)
Hcf_so_ss_tet = Hcf_so_tet + Hss_tet

Ecf_so_ss_val_tet,Ecf_so_ss_val_excitation_tet,H_cf_so_ss_vt_tet = Fe2_tet.Diag(Hcf_so_ss_tet,printfunction=True)

#################################### Molecular Field ##########################################


print('---------------- Crystal Field with Spin-Orbit coupling and Molecular Field ---------------')

Hm = 18.2*8.065548

Hg_tet = Fe2_tet.Molecular_Field_Sz(ion,L,S,Hm)

Hcf_so_ss_g_tet = Hcf_so_ss_tet + Hg_tet

Ecf_so_ss_g_val_tet,Ecf_so_ss_g_val_excitation_tet,H_cf_so_ss_g_vt_tet = Fe2_tet.Diag(Hcf_so_ss_g_tet,printfunction=True)
#################################### Magnetic Moment ##########################################
"""
print('\n')



MagMoment_L_z,MagMoment_S_z,MagMoment_M_z = Fe2_tet.MagMoment_z(L,S,Hcf_so)

MagMoment_L_x,MagMoment_S_x,MagMoment_M_x = Fe2_tet.MagMoment_x(L,S,Hcf_so)

MagMoment_L_y,MagMoment_S_y,MagMoment_M_y = Fe2_tet.MagMoment_y(L,S,Hcf_so)




################## specific heat ###############

C_schottky,T = Fe2_tet.SpecificHeat(Ecf_val_excitation)
C_schottky,T = Fe2_tet.SpecificHeat(Ecf_so_val_excitation)
C_schottky,T = Fe2_tet.SpecificHeat(Ecf_so_g_val_excitation)


################## magnetic entropy ###############

entropy,T = Fe2_tet.Entropy(Ecf_val_excitation)
entropy,T = Fe2_tet.Entropy(Ecf_so_val_excitation)
entropy,T = Fe2_tet.Entropy(Ecf_so_g_val_excitation)


############# Magnetisation ###############

Fields = 5 # T


Zeeman = Fe2_tet.Zeeman_z(L,S,Hcf_so_g,Fields)
Zeeman = Fe2_tet.Zeeman_x(L,S,Hcf_so_g,Fields)
Zeeman = Fe2_tet.Zeeman_y(L,S,Hcf_so_g,Fields)


############# Magnetic susceptibility ###############

X_z,T = Fe2_tet.Susceptibility_z(L,S,Hcf_so,SO)
X_x,T = Fe2_tet.Susceptibility_x(L,S,Hcf_so,SO)
X_y,T = Fe2_tet.Susceptibility_y(L,S,Hcf_so,SO)






"""














# -*- coding: utf-8 -*-
"""
Created on Sun Sep  6 22:01:34 2020

@author: brian
"""
import numpy as np
import Crystal_Field_Calculations as cef
import matplotlib.pyplot as plt

np.set_printoptions(suppress=True)

ion = 'Fe2'
L = 2
S = 2
Z = 2
Fe2_oct = cef.LS(L,S)


Fe_oct = np.array([0.000000,    3.333159,    5.157905])
O1 = np.array([1.341980,    2.558366,    3.648669])
O2 = np.array([0.000000,    4.882744,    3.648669])
O3 = np.array([-1.341980,   2.558366,    3.648669])
O4 = np.array([-1.439547,    4.164282,    6.378385])
O5 = np.array([0.000000,    1.670912,    6.378385])
O6 = np.array([1.439547,    4.164282,    6.378385])


O1_d = O1 - Fe_oct
O2_d = O2 - Fe_oct
O3_d = O3 - Fe_oct
O4_d = O4 - Fe_oct
O5_d = O5 - Fe_oct
O6_d = O6 - Fe_oct


d = np.array([O1_d,O2_d,O3_d,O4_d,O5_d,O6_d])

B = Fe2_oct.PC(ion,L,S,d,Z)


#rho = 0.1
#Qeff = 0.8
#B = Fe2_oct.SOM(ion,L,S,d,Z,rho,Qeff)

B02 = B[0]
#B02_oct = B02
B04 = B[3]
#B04_oct = B04
B34 = B[6]
#B34_oct = B34

print('B02 (meV): ' + str(B02))
print('B04 (meV): ' + str(B04))
print('B34 (meV): ' + str(B34))

#################################### Crystal Field ##########################################

print('---------------- From literature ---------------')
O20 = Fe2_oct.Olm(L,S,2,0)
#print('O02: ' + str(O02))
#print(O02_diag)

O40 = Fe2_oct.Olm(L,S,4,0)
#print('O04: ' + str(O04))
#print(O04_diag)

O43 = Fe2_oct.Olm(L,S,4,3)


Dq = -685 # cm-1
A2 = -596.0 # cm-1
A4 = -332.0 # cm-1

B02 = (A2/3)
B04 = (A4/12 + Dq/18)
B34 = (-20*Dq/(9*np.sqrt(2)))

print('B02 (meV): ' + str(B02))
print('B04 (meV): ' + str(B04))
print('B34 (meV): ' + str(B34))

print('---------------- Crystal Field ---------------')

Hcf_oct = B02*O20 + B04*O40 + B34*O43
Ecf_val_oct,Ecf_val_excitation_oct,H_cf_vt_oct = Fe2_oct.Diag(Hcf_oct,printfunction=True)


#################################### Spin-Orbit coupling ##########################################
print('---------------- Crystal Field with Spin-Orbit coupling ---------------')


SO = -116.5

SO_matrix = Fe2_oct.SO(ion,L,S,SO)

Hcf_so_oct = Hcf_oct + SO_matrix

Ecf_so_val_oct,Ecf_so_val_excitation_oct,H_cf_so_vt_oct = Fe2_oct.Diag(Hcf_so_oct,printfunction=True)

#################################### Spin-Spin interaction ##########################################
print('---------------- Crystal Field with Spin-Orbit coupling, spin-spin interaction, and Molecular Field ---------------')


p = 1.1


Hss_oct = Fe2_oct.SS(ion,L,S,p)

Hcf_so_ss_oct = Hcf_so_oct + Hss_oct

Ecf_so_ss_val_oct,Ecf_so_ss_val_excitation_oct,H_cf_so_ss_vt_oct = Fe2_oct.Diag(Hcf_so_ss_oct,printfunction=True)

#################################### Molecular Field ##########################################

print('---------------- Crystal Field with Spin-Orbit coupling and Molecular Field ---------------')


Hm = 13.6*8.065548


Hg_oct = Fe2_oct.Molecular_Field_Sz(ion,L,S,Hm)

Hcf_so_ss_g_oct = Hcf_so_ss_oct + Hg_oct

Ecf_so_ss_g_val_oct,Ecf_so_ss_g_val_excitation_oct,H_cf_so_ss_g_vt_oct = Fe2_oct.Diag(Hcf_so_ss_g_oct,printfunction=True)







"""
#################################### Magnetic Moment ##########################################
print('--------- Magnetic Properties --------')




MagMoment_L_z,MagMoment_S_z,MagMoment_M_z = Fe2_oct.MagMoment_z(L,S,Hcf_so)

MagMoment_L_x,MagMoment_S_x,MagMoment_M_x = Fe2_oct.MagMoment_x(L,S,Hcf_so)

MagMoment_L_y,MagMoment_S_y,MagMoment_M_y = Fe2_oct.MagMoment_y(L,S,Hcf_so)



################## specific heat ###############

C_schottky,T = Fe2_oct.SpecificHeat(Ecf_val_excitation)
C_schottky,T = Fe2_oct.SpecificHeat(Ecf_so_val_excitation)
C_schottky,T = Fe2_oct.SpecificHeat(Ecf_so_g_val_excitation)


################## magnetic entropy ###############

entropy,T = Fe2_oct.Entropy(Ecf_val_excitation)
entropy,T = Fe2_oct.Entropy(Ecf_so_val_excitation)
entropy,T = Fe2_oct.Entropy(Ecf_so_g_val_excitation)


############# Zeeman effect ###############

Fields = 5 # T


Zeeman = Fe2_oct.Zeeman_z(L,S,Hcf_so_g,Fields)
Zeeman = Fe2_oct.Zeeman_x(L,S,Hcf_so_g,Fields)
Zeeman = Fe2_oct.Zeeman_y(L,S,Hcf_so_g,Fields)


############# Magnetic susceptibility ###############

X_z,T = Fe2_oct.Susceptibility_z(L,S,Hcf_so,SO)
X_x,T = Fe2_oct.Susceptibility_x(L,S,Hcf_so,SO)
X_y,T = Fe2_oct.Susceptibility_y(L,S,Hcf_so,SO)

X_z,T = Fe2_oct.Susceptibility_z(L,S,Hcf_so,SO)
X_x,T = Fe2_oct.Susceptibility_x(L,S,Hcf_so,SO)
X_y,T = Fe2_oct.Susceptibility_y(L,S,Hcf_so,SO)

plt.plot(T,X_z,'b.',T,X_x,'r.',T,X_y,'y.')
plt.xlim(0,300)
#plt.ylim(0,2)
plt.xlabel('Temperature (K)')
plt.ylabel('Magnetic Susceptibility ($\mu_B$)')
plt.title('Magnetic Susceptibility')
#plt.legend('')
plt.grid(alpha=.4,linestyle='--')
plt.show()



"""


