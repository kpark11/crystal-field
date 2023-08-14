# -*- coding: utf-8 -*-
"""
Created on Sun Sep  6 22:01:34 2020

@author: brian
"""
import numpy as np
import Crystal_Field_Calcuations as cef
import matplotlib.pyplot as plt
from scipy import optimize

######## Literature Yang et al, Physica B 370 (2005). #######
#  B20 = v - 2*sqrt(2)*v'
#  B40 = (4/3)*v + 2*sqrt(2)*v' - 14*Dq
#  B43 = -(7/sqrt(70))*((2/3)*v + sqrt(2)*v' + 20*Dq)



np.set_printoptions(suppress=True)

ion = 'Ni2'
L = 3
S = 1
Z = 2
Ni2_oct = cef.LS(L,S)




Ni_pos_oct = np.array([2.874350,    1.659507,    1.156248])
O_pos1_oct = np.array([4.297096,    0.838084,    2.389448])
O_pos2_oct = np.array([2.874350,    3.302352,    2.389448])
O_pos3_oct = np.array([1.451604,    0.838084,    2.389448])
O_pos4_oct = np.array([1.531252,    2.434945,   -0.259761])
O_pos5_oct = np.array([2.874350,    0.108631,   -0.259761])
O_pos6_oct = np.array([4.217448,    2.434945,   -0.259761])



O1_oct = O_pos1_oct - Ni_pos_oct
O2_oct = O_pos2_oct - Ni_pos_oct
O3_oct = O_pos3_oct - Ni_pos_oct
O4_oct = O_pos4_oct - Ni_pos_oct
O5_oct = O_pos5_oct - Ni_pos_oct
O6_oct = O_pos6_oct - Ni_pos_oct

d_oct = np.array([O1_oct,O2_oct,O3_oct,O4_oct,O5_oct,O6_oct])

B_oct = Ni2_oct.PC(ion,L,S,d_oct,Z)

#################################### Crystal Field ##########################################
print('---------------- From literature ---------------')
##### French study



"""
Dq_oct = 875.5 # cm-1
A2_oct = -97.3 # cm-1
A4_oct = -579.3 # cm-1


B02_oct = (A2_oct/3)/8.065548
B04_oct = (A4_oct/12 + Dq_oct/18)/8.065548
B34_oct = (-20*Dq_oct/(9*np.sqrt(2)))/8.065548
"""
##################  Maybe #####################
Dq_oct = 110.6 #meV
B34_oct = (-20*Dq_oct/(9*np.sqrt(2)))/20
B04_oct = -np.sqrt(2)*(B34_oct)/(9*np.sqrt(40))
B02_oct = 200 #4.31 # meV


B02_oct = 15.51/8.065548
B04_oct = -0.38/8.065548
B34_oct = -9.28/8.065548

print("B34_oct: " + str(B34_oct))
print("B04_oct: " + str(B04_oct))
print("B02_oct: " + str(B02_oct))
###########3 To play around


#B04_oct = 0.4158937711201375  
#B34_oct = -8.68955666658135  
#B02_oct = 125  
"""
### From Radwanski

B04_oct = 21/8.065548  # in K, but 1 K = 1 cm-1, which is why i am converting to meV the same way
B02_oct = 50/8.065548
Dq_oct = 60*B04_oct
B34_oct = (-20*Dq_oct/(9*np.sqrt(2)))/8.065548


### From me
Dq_oct = 806.5 # cm-1
B04_oct = Dq_oct/(60*8.065548)
B02_oct = 50/8.065548

B34_oct = (-20*Dq_oct/(9*np.sqrt(2)))/8.065548
"""


print('---------------- Crystal Field Parameters ---------------')
print('B02 (meV): ' + str(B02_oct))
print('B04 (meV): ' + str(B04_oct))
print('B34 (meV): ' + str(B34_oct))

#################################### Crystal Field ##########################################
O20 = Ni2_oct.Olm(L,S,2,0)
O21 = Ni2_oct.Olm(L,S,2,1)
O22 = Ni2_oct.Olm(L,S,2,2)
O40 = Ni2_oct.Olm(L,S,4,0)
O41 = Ni2_oct.Olm(L,S,4,1)
O42 = Ni2_oct.Olm(L,S,4,2)
O43 = Ni2_oct.Olm(L,S,4,3)
O44 = Ni2_oct.Olm(L,S,4,4)
O2m1 = Ni2_oct.Olm(L,S,2,-1)
O2m2 = Ni2_oct.Olm(L,S,2,-2)
O4m1 = Ni2_oct.Olm(L,S,4,-1)
O4m2 = Ni2_oct.Olm(L,S,4,-2)
O4m3 = Ni2_oct.Olm(L,S,4,-3)
O4m4 = Ni2_oct.Olm(L,S,4,-4)

#B_oct = np.array([ 3.47822028,  2.13070696, -7.52546629])
#print(B_oct)


Hcf_oct = B02_oct*O20 + B04_oct*O40 + B34_oct*O43
#Hcf_oct = B_oct[0]*O20 + B_oct[1]*O40 + B_oct[2]*O43
Ecf_val_oct,Ecf_val_excitation_oct,H_cf_vt_oct = Ni2_oct.Diag(Hcf_oct,printfunction=True)



#################################### Spin-Orbit coupling ##########################################
print('---------------- Crystal Field with Spin-Orbit coupling ---------------')

SO_oct = 41

SO_oct = -303.38/8.065548

SO_matrix_oct = Ni2_oct.SO(ion,L,S,SO_oct)

Hcf_so_oct = Hcf_oct + SO_matrix_oct

Ecf_so_val_oct,Ecf_so_val_excitation_oct,H_cf_so_vt_oct = Ni2_oct.Diag(Hcf_so_oct,printfunction=True)

#################################### Spin-Spin interaction ##########################################
"""

p = 1.1


Hss = Fe2p.SS(ion,L,S,p)
Hcf_so_ss = Hcf_so + Hss
"""
#################################### Molecular Field ##########################################


print('---------------- Crystal Field with Spin-Orbit coupling and Molecular Field ---------------')

Hm=23.3/8.065548

Hg_oct = Ni2_oct.Molecular_Field_Sz(ion,L,S,Hm)

Hcf_so_g_oct = Hcf_so_oct + Hg_oct

Ecf_so_g_val_oct,Ecf_so_g_val_excitation_oct,H_cf_so_g_vt_oct = Ni2_oct.Diag(Hcf_so_g_oct,printfunction=True)

#################################### Magnetic Moment ##########################################

print('\n')



MagMoment_L_z_oct,MagMoment_S_z_oct,MagMoment_M_z_oct = Ni2_oct.MagMoment_z(L,S,Hcf_so_g_oct)

MagMoment_L_x_oct,MagMoment_S_x_oct,MagMoment_M_x_oct = Ni2_oct.MagMoment_x(L,S,Hcf_so_g_oct)

MagMoment_L_y_oct,MagMoment_S_y_oct,MagMoment_M_y_oct = Ni2_oct.MagMoment_y(L,S,Hcf_so_g_oct)

"""
################## specific heat ###############

C_schottky_oct,T = Ni2_oct.SpecificHeat(Ecf_val_excitation_oct)
C_schottky_oct,T = Ni2_oct.SpecificHeat(Ecf_so_val_excitation_oct)


################## magnetic entropy ###############

entropy_oct,T = Ni2_oct.Entropy(Ecf_val_excitation_oct)
entropy_oct,T = Ni2_oct.Entropy(Ecf_so_val_excitation_oct)


############# Magnetisation ###############
#Fields = 5 # T


#Zeeman = Ni2_oct.Zeeman_z(L,S,Hcf_so_g,Fields)
#Zeeman = Ni2_oct.Zeeman_x(L,S,Hcf_so_g,Fields)
#Zeeman = Ni2_oct.Zeeman_y(L,S,Hcf_so_g,Fields)

############# Magnetic susceptibility ###############


X_z_oct,T = Ni2_oct.Susceptibility_z(L,S,Hcf_so_oct,SO_oct)
X_x_oct,T = Ni2_oct.Susceptibility_x(L,S,Hcf_so_oct,SO_oct)
X_y_oct,T = Ni2_oct.Susceptibility_y(L,S,Hcf_so_oct,SO_oct)


plt.plot(T,X_z_oct,'b.',T,X_x_oct,'r.',T,X_y_oct,'y.')
plt.xlim(0,300)
#plt.ylim(0,2)
plt.xlabel('Temperature (K)')
plt.ylabel('Magnetic Susceptibility ($\mu_B$)')
plt.title('Magnetic Susceptibility')
#plt.legend('')
plt.grid(alpha=.4,linestyle='--')
plt.show()




Exp_oct = np.array([0.00001,0.00001,0.00001,0.00001,0.00001,0.00001,1.06,1.06,1.06,1.06,1.06,1.06,1.141,1.141,1.141,1.494,1.494,1.494,1.566,1.566,1.566])*1000
def Fit_B(B,E_exp):
        #Qeff = X[0]
        #rho = X[1]
        #SO = B[8]
        #SOM_2 = rho*((2/(1-rho))**(3))*(Qeff/Z)*B[0:3]
        #SOM_4 = rho*((2/(1-rho))**(5))*(Qeff/Z)*B[3:8]
        #SOM = np.array([])
        #SOM = np.append(SOM,SOM_2)
        #SOM = np.append(SOM,SOM_4)
        #SO_matrix = Cr3_oct.SO(ion,L,S,SO)
        Hcf_oct = B[0]*O20 + B[1]*O40 + B[2]*O43
        E = np.linalg.eigh(Hcf_oct)
        E_val = E[0]
        #E_state = E[1].T
        E_excitation = E_val-E_val[0]
        
        return np.sum((E_excitation[6:]/E_exp[6:] - 1)**2)

result = optimize.minimize(Fit_B, B_oct, args = (Exp_oct), method = 'Powell',tol=1e-8)
B_oct = result.x
print(result)

"""




#######################################################################################################
#######################################################################################################
#######################################################################################################
#######################################################################################################
#######################################################################################################
#######################################################################################################
#######################################################################################################
#######################################################################################################


ion = 'Ni2'
L = 3
S = 1
Z = 2
Ni2_tet = cef.LS(L,S)





Ni_tet = np.array([2.874350,    1.659507,    5.471264])
O1 = np.array([2.874350,    1.659507,    7.450501])
O2 = np.array([1.343098,    2.543576,    4.673039])
O3 = np.array([2.874350,   -0.108631,    4.673039])
O4 = np.array([4.405602,    2.543576,    4.673039])



O1_d = O1 - Ni_tet
O2_d = O2 - Ni_tet
O3_d = O3 - Ni_tet
O4_d = O4 - Ni_tet


d_tet = np.array([O1_d,O2_d,O3_d,O4_d])

B_tet = Ni2_tet.PC(ion,L,S,d_tet,Z)

"""
####### From Allen Scheie


NiPosTet = np.array([0.33333,  0.66667, -0.05200])

OposTet = np.array([[0.33333,  0.66667,  0.14610],
	[-0.02400,  0.48800, -0.13410],
	[0.51200,  0.48800, -0.13410],
	[0.51200,  1.02400, -0.13410]])

O1_d = OposTet[0] - NiPosTet
O2_d = OposTet[1] - NiPosTet
O3_d = OposTet[2] - NiPosTet
O4_d = OposTet[3] - NiPosTet


d_tet = np.array([O1_d,O2_d,O3_d,O4_d])

B_tet = Ni2_tet.PC(ion,L,S,d_tet,Z)

"""

#rho = 0.11
#Qeff = 0.9
#B_tet = Ni2_tet.SOM(ion,L,S,d_tet,Z,rho,Qeff)

#################################### Crystal Field ##########################################



O20 = Ni2_tet.Olm(L,S,2,0)
#print('O02: ' + str(O02))
#print(O02_diag)

O40 = Ni2_tet.Olm(L,S,4,0)
#print('O04: ' + str(O04))
#print(O04_diag)

O43 = Ni2_tet.Olm(L,S,4,3)


B02_tet = B_tet[0]
#B02_tet  = B02
B04_tet = B_tet[3]
#B04_tet = B04
B34_tet = B_tet[6]
#B34_tet = B34



#################################### Crystal Field ##########################################
print('---------------- From literature ---------------')
##### French study



"""
Dq_tet = 241.9 # cm-1
A2_tet = -5.43 # cm-1
A4_tet = -500.3 # cm-1

B02_tet = (A2_tet/3)/8.065548
B04_tet = (A4_tet/12 + Dq_tet/18)/8.065548
B34_tet = (-20*Dq_tet/(9*np.sqrt(2)))/8.065548


print(B02_tet)
print(B04_tet)
print(B34_tet)


### From Radwanski

Dq_tet = 63
B04_tet = Dq_tet/60
B02_tet = 50/8.065548
B34_tet = (-20*Dq_tet/(9*np.sqrt(2)))/8.065548


### From me

Dq_tet = 25
B04_tet = Dq_tet/60
B02_tet = 20/8.065548
B34_tet = (-20*Dq_tet/(9*np.sqrt(2)))/8.065548
"""
############ this might be it
Dq_tet = 25 # meV
B34_tet = (-20*Dq_tet/(9*np.sqrt(2)))
B04_tet = -np.sqrt(2)*(B34_tet)/(9*np.sqrt(40))
B02_tet = 28 #4.31 # meV

##########3 To play around

#Dq_tet = 15 # meV


B02_tet = 7.88/8.065548
B04_tet = 0.065/8.065548
B34_tet = 9.66/8.065548


print('---------------- Crystal Field Parameters ---------------')
print('B02 (meV): ' + str(B02_tet))
print('B04 (meV): ' + str(B04_tet))
print('B34 (meV): ' + str(B34_tet))

print('---------------- Crystal Field ---------------')

#B_tet = np.array([15.40479649,  0.30189252, 10.35892682])
#print(B_tet)


Hcf_tet = B02_tet*O20 + B04_tet*O40 + B34_tet*O43
#Hcf_tet = B2m2_tet*O2m2 + B2m1_tet*O2m1 + B20_tet*O20 + B21_tet*O21 + B22_tet*O22
#Hcf_tet = B4m4_tet*O4m4 + B4m3_tet*O4m3 + B4m2_tet*O4m2 + B4m1_tet*O4m1 + B40_tet*O40 + Hcf_tet
#Hcf_tet = B44_tet*O44 + B43_tet*O43 + B42_tet*O42 + B41_tet*O41 + Hcf_tet
Ecf_val_tet,Ecf_val_excitation_tet,H_cf_vt_tet = Ni2_tet.Diag(Hcf_tet,printfunction=True)



#################################### Spin-Orbit coupling ##########################################
print('---------------- Crystal Field with Spin-Orbit coupling ---------------')


SO_tet = -332.53/8.065548

SO_matrix_tet = Ni2_tet.SO(ion,L,S,SO_tet)

Hcf_so_tet = Hcf_tet + SO_matrix_tet

Ecf_so_val_tet,Ecf_so_val_excitation_tet,H_cf_so_vt_tet = Ni2_tet.Diag(Hcf_so_tet,printfunction=True)

#################################### Spin-Spin interaction ##########################################
"""

p = 1.1


Hss = Fe2p.SS(ion,L,S,p)
Hcf_so_ss = Hcf_so + Hss
"""

#################################### Molecular Field ##########################################


print('---------------- Crystal Field with Spin-Orbit coupling and Molecular Field ---------------')

Hm = 10

Hm=13.24/8.065548


Hg_tet = Ni2_tet.Molecular_Field_Sz(ion,L,S,Hm)

Hcf_so_g_tet = Hcf_so_tet + Hg_tet

Ecf_so_g_val_tet,Ecf_so_g_val_excitation_tet,H_cf_so_g_vt_tet = Ni2_tet.Diag(Hcf_so_g_tet,printfunction=True)


#################################### Magnetic Moment ##########################################

print('\n')



MagMoment_L_z_tet,MagMoment_S_z_tet,MagMoment_M_z_tet = Ni2_tet.MagMoment_z(L,S,Hcf_so_g_tet)

MagMoment_L_x_tet,MagMoment_S_x_tet,MagMoment_M_x_tet = Ni2_tet.MagMoment_x(L,S,Hcf_so_g_tet)

MagMoment_L_y_tet,MagMoment_S_y_tet,MagMoment_M_y_tet = Ni2_tet.MagMoment_y(L,S,Hcf_so_g_tet)

'''
################## specific heat ###############4

C_schottky_tet,T = Ni2_tet.SpecificHeat(Ecf_val_excitation_tet)
C_schottky_tet,T = Ni2_tet.SpecificHeat(Ecf_so_val_excitation_tet)


################## magnetic entropy ###############

entropy_tet,T = Ni2_tet.Entropy(Ecf_val_excitation_tet)
entropy_tet,T = Ni2_tet.Entropy(Ecf_so_val_excitation_tet)


############# Magnetisation ###############
#Fields = 5 # T


#Zeeman = Ni2_tet.Zeeman_z(L,S,Hcf_so_g,Fields)
#Zeeman = Ni2_tet.Zeeman_x(L,S,Hcf_so_g,Fields)
#Zeeman = Ni2_tet.Zeeman_y(L,S,Hcf_so_g,Fields)

############# Magnetic susceptibility ###############
X_z_tet,T = Ni2_tet.Susceptibility_z(L,S,Hcf_so_tet,SO_tet)
X_x_tet,T = Ni2_tet.Susceptibility_x(L,S,Hcf_so_tet,SO_tet)
X_y_tet,T = Ni2_tet.Susceptibility_y(L,S,Hcf_so_tet,SO_tet)




plt.plot(T,X_z_tet,'b.',T,X_x_tet,'r.',T,X_y_tet,'y.')
plt.xlim(0,300)
#plt.ylim(0,2)
plt.xlabel('Temperature (K)')
plt.ylabel('Magnetic Susceptibility ($\mu_B$)')
plt.title('Magnetic Susceptibility')
#plt.legend('')
plt.grid(alpha=.4,linestyle='--')
plt.show()





Exp_tet = np.array([0,0,0,0.0192,0.0192,0.0192,0.0192,0.0192,0.0192,0.2,0.2,0.2,0.255,0.255,0.255,0.554,0.554,0.554,0.64,0.64,0.64])*1000
def Fit_B(B,E_exp):
        #Qeff = X[0]
        #rho = X[1]
        #SO = B[8]
        #SOM_2 = rho*((2/(1-rho))**(3))*(Qeff/Z)*B[0:3]
        #SOM_4 = rho*((2/(1-rho))**(5))*(Qeff/Z)*B[3:8]
        #SOM = np.array([])
        #SOM = np.append(SOM,SOM_2)
        #SOM = np.append(SOM,SOM_4)
        #SO_matrix = Cr3_oct.SO(ion,L,S,SO)
        Hcf_oct = B[0]*O20 + B[1]*O40 + B[2]*O43
        E = np.linalg.eigh(Hcf_oct)
        E_val = E[0]
        #E_state = E[1].T
        E_excitation = E_val-E_val[0]
        
        return np.sum((E_excitation[3:]/E_exp[3:] - 1)**2)

result = optimize.minimize(Fit_B, B_tet, args = (Exp_tet), method = 'Powell',tol=1e-12  )
B_tet = result.x
print(result)

'''
