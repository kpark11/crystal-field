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


Fe = np.array([1/3, 2/3, -0.0481])

O1 = np.array([-0.0234, 0.4883, -0.1371])

O2 = np.array([0.5117, 1.0234, -0.1371])

O3 = np.array([0.5117, 0.4883, -0.1371])

O4 = np.array([1/3, 2/3, 0.147])


O1_d = O1 - Fe
O2_d = O2 - Fe
O3_d = O3 - Fe
O4_d = O4 - Fe


d = np.array([O1_d,O2_d,O3_d,O4_d])
print('---------------- Ligand positions ---------------')
print(d)


rho = 0.1
Qeff = 0.7
B = Fe2p.PC(ion,L,S,d,Z)
B = Fe2p.SOM(ion,L,S,d,Z,rho,Qeff)


#################################### Crystal Field ##########################################
print('---------------- From literature ---------------')



O20 = Fe2p.Olm(L,S,2,0)
#print('O02: ' + str(O02))
#print(O02_diag)

O40 = Fe2p.Olm(L,S,4,0)
#print('O04: ' + str(O04))
#print(O04_diag)

O43 = Fe2p.Olm(L,S,4,3)


Dq = 302.8 # cm-1
A2 = -33.3 # cm-1
A4 = -399.4 # cm-1

B02 = (A2/3)/8.065548
B04 = (A4/12 + Dq/18)/8.065548
B34 = (-20*Dq/(9*np.sqrt(2)))/8.065548

print('B02 (meV): ' + str(B02))
print('B04 (meV): ' + str(B04))
print('B34 (meV): ' + str(B34))

print('---------------- Crystal Field ---------------')

Hcf = B02*O20 + B04*O40 + B34*O43
Ecf_val,Ecf_val_excitation,H_cf_vt = Fe2p.Diag(Hcf)



#################################### Spin-Orbit coupling ##########################################
print('---------------- Crystal Field with Spin-Orbit coupling ---------------')


SO = -116.5/8.065548

SO_matrix = Fe2p.SO(ion,L,S,SO)

Hcf_so = Hcf + SO_matrix

Ecf_so_val,Ecf_so_val_excitation,H_cf_so_vt = Fe2p.Diag(Hcf_so)

#################################### Spin-Spin interaction ##########################################


p = 1.1


Hss = Fe2p.SS(ion,L,S,p)
Hcf_so_ss = Hcf_so + Hss


#################################### Molecular Field ##########################################


print('---------------- Crystal Field with Spin-Orbit coupling and Molecular Field ---------------')

Hm = 18.2/8.065548

Hg = Fe2p.Molecular_Field_Sz(ion,L,S,Hm)

Hcf_so_g = Hcf_so + Hg

Ecf_so_g_val,Ecf_so_g_val_excitation,H_cf_so_g_vt = Fe2p.Diag(Hcf_so_g)
#################################### Magnetic Moment ##########################################

print('\n')



MagMoment_L_z,MagMoment_S_z,MagMoment_M_z = Fe2p.MagMoment_z(L,S,Hcf_so)

MagMoment_L_x,MagMoment_S_x,MagMoment_M_x = Fe2p.MagMoment_x(L,S,Hcf_so)

MagMoment_L_y,MagMoment_S_y,MagMoment_M_y = Fe2p.MagMoment_y(L,S,Hcf_so)


################## specific heat ###############

C_schottky,T = Fe2p.SpecificHeat(Ecf_val_excitation)
C_schottky,T = Fe2p.SpecificHeat(Ecf_so_val_excitation)
C_schottky,T = Fe2p.SpecificHeat(Ecf_so_g_val_excitation)


################## magnetic entropy ###############

entropy,T = Fe2p.Entropy(Ecf_val_excitation)
entropy,T = Fe2p.Entropy(Ecf_so_val_excitation)
entropy,T = Fe2p.Entropy(Ecf_so_g_val_excitation)


############# Magnetisation ###############
Fields = 5 # T


Zeeman = Fe2p.Zeeman_z(L,S,Hcf_so_g,Fields)
Zeeman = Fe2p.Zeeman_x(L,S,Hcf_so_g,Fields)
Zeeman = Fe2p.Zeeman_y(L,S,Hcf_so_g,Fields)

############# Magnetic susceptibility ###############

X_z,T = Fe2p.Susceptibility_z(L,S,Hcf_so,SO)
X_x,T = Fe2p.Susceptibility_x(L,S,Hcf_so,SO)
X_y,T = Fe2p.Susceptibility_y(L,S,Hcf_so,SO)

plt.plot(T,X_z,'b.',T,X_x,'r.',T,X_y,'y.')
plt.xlim(0,300)
#plt.ylim(0,2)
plt.xlabel('Temperature (K)')
plt.ylabel('Magnetic Susceptibility ($\mu_B$)')
plt.title('Magnetic Susceptibility')
#plt.legend('')
plt.grid(alpha=.4,linestyle='--')
plt.show()

plt.plot(T,1/X_z,'bo',T,1/X_x,'ro',T,1/X_y,'yo')
plt.xlim(0,300)
plt.ylim(0,60)
plt.xlabel('Temperature (K)')
plt.ylabel('Inverse Magnetic Susceptibility ($\mu_B$)')
plt.title('Inverse Magnetic Susceptibility')
#plt.legend('')
plt.grid(alpha=.4,linestyle='--')
plt.show()

        
X_z_g,T = Fe2p.Susceptibility_z(L,S,Hcf_so_g)
X_x_g,T = Fe2p.Susceptibility_x(L,S,Hcf_so_g)
X_y_g,T = Fe2p.Susceptibility_y(L,S,Hcf_so_g)

plt.plot(T,X_z_g,'b.',T,X_x_g,'r.',T,X_y_g,'y.')
plt.xlim(0,300)
#plt.ylim(0,2)
plt.xlabel('Temperature (K)')
plt.ylabel('Magnetic Susceptibility ($\mu_B$)')
plt.title('Magnetic Susceptibility')
#plt.legend('')
plt.grid(alpha=.4,linestyle='--')
plt.show()

plt.plot(T,1/X_z_g,'bo',T,1/X_x_g,'ro',T,1/X_y_g,'yo')
plt.xlim(0,300)
#plt.ylim(0,60)
plt.xlabel('Temperature (K)')
plt.ylabel('Inverse Magnetic Susceptibility ($\mu_B$)')
plt.title('Inverse Magnetic Susceptibility')
#plt.legend('')
plt.grid(alpha=.4,linestyle='--')
plt.show()


