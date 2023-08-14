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

ion = 'Mn2'
L = 0
S = 5/2
Z = 2
Mn2_oct = cef.LS(L,S)



Mn = np.array([3.034000,    8.750595,    0.000000])
S1 = np.array([0.987336,    8.797884,    1.614580])
S2 = np.array([4.021336,    6.965167,    1.614580])
S3 = np.array([4.089697,   10.508700,    1.613166])
S4 = np.array([2.046664,    6.965167,   -1.614580])
S5 = np.array([5.080664,    8.797884,   -1.614580])
S6 = np.array([1.978303,   10.508700,   -1.613166])


S1_d = S1 - Mn
S2_d = S2 - Mn
S3_d = S3 - Mn
S4_d = S4 - Mn
S5_d = S5 - Mn
S6_d = S6 - Mn

d = np.array([S1_d,S2_d,S3_d,S4_d,S5_d,S6_d])

#print(d)

########################################### Crystal Field Splitting ########################################


B = Mn2_oct.PC(ion,L,S,d,Z)

rho = 0.3
Qeff = 1.2
B = Mn2_oct.SOM(ion,L,S,d,Z,rho,Qeff)


O20 = Mn2_oct.Olm(L,S,2,0)
O21 = Mn2_oct.Olm(L,S,2,1)
O22 = Mn2_oct.Olm(L,S,2,2)
O40 = Mn2_oct.Olm(L,S,4,0)
O41 = Mn2_oct.Olm(L,S,4,1)
O42 = Mn2_oct.Olm(L,S,4,2)
O43 = Mn2_oct.Olm(L,S,4,3)
O44 = Mn2_oct.Olm(L,S,4,4)



############### From Literature ###########

#B[0] = -350.4/8.065548
B[0] = -450.4/8.065548
#B[2] = -140.2/8.065548
B[2] = -100.2/8.065548
B[3] = -6.7/8.065548
B[5] = 71.8/8.065548
B[7] = 50.9/8.065548

print("#################################################")

print("This is the value from literature: ")
print("B20: " + str(B[0]))
print("B22: " + str(B[2]))
print("B40: " + str(B[3]))
print("B42: " + str(B[5]))
print("B44: " + str(B[7]))

print("#################################################")


#Hcf = B[0]*O20 + B[1]*O21 + B[2]*O22 + B[3]*O40 + B[4]*O41 + B[5]*O42 + B[6]*O43 + B[7]*O44
Hcf = B[0]*O20 + B[2]*O22 + B[3]*O40 + B[5]*O42 + B[7]*O44


Ecf_val,Ecf_val_excitation,H_cf_vt = Mn2_oct.Diag(Hcf,printfunction=True)


####################################### Spin-Orbit Coupling #######################################

SO = 31 #meV

SO_matrix = Mn2_oct.SO(ion,L,S)

Hcf_so = Hcf + SO_matrix

Ecf_so_val,Ecf_so_val_excitation,H_cf_so_vt = Mn2_oct.Diag(Hcf_so,printfunction=True)



#################################### Magnetic Moment ##########################################
print('--------- Magnetic Properties --------')




MagMoment_L_z,MagMoment_S_z,MagMoment_M_z = Mn2_oct.MagMoment_z(L,S,Hcf_so)

MagMoment_L_x,MagMoment_S_x,MagMoment_M_x = Mn2_oct.MagMoment_x(L,S,Hcf_so)

MagMoment_L_y,MagMoment_S_y,MagMoment_M_y = Mn2_oct.MagMoment_y(L,S,Hcf_so)



################## specific heat ###############

C_schottky,T = Mn2_oct.SpecificHeat(Ecf_val_excitation)
C_schottky,T = Mn2_oct.SpecificHeat(Ecf_so_val_excitation)


################## magnetic entropy ###############

entropy,T = Mn2_oct.Entropy(Ecf_val_excitation)
entropy,T = Mn2_oct.Entropy(Ecf_so_val_excitation)


############# Zeeman effect ###############
"""
Fields = 25 # T


Zeeman = Cr3_oct.Zeeman_z(L,S,Hcf_so,Fields,savefunction=False)
Zeeman = Cr3_oct.Zeeman_x(L,S,Hcf_so,Fields,savefunction=False)
Zeeman = Cr3_oct.Zeeman_y(L,S,Hcf_so,Fields,savefunction=False)

"""
############# Magnetic susceptibility ###############

X_z,T = Mn2_oct.Susceptibility_z(L,S,Hcf_so,SO)
X_x,T = Mn2_oct.Susceptibility_x(L,S,Hcf_so,SO)
X_y,T = Mn2_oct.Susceptibility_y(L,S,Hcf_so,SO)



plt.plot(T,X_z,'ro')
plt.plot(T,X_x,'bo')
plt.plot(T,X_y,'yo')
plt.xlabel('Temperature (K)')
plt.ylabel('Magnetic Susceptibility ($\mu_B$)')
plt.title('Magnetic Susceptibility')
#plt.legend('')
plt.grid(alpha=.4,linestyle='--')
plt.xlim(0,300)
plt.ylim(0,0.18)
plt.show()


exp_susceptibility = np.loadtxt(r'C:\Users\Student\OneDrive - University of Tennessee\Desktop\Research\MPS (M = Cr, Mn)\CrPS4\Digitized_susceptibility_c.dat')

"""
def Fit_B(B,X,E_exp):
        #Qeff = X[0]
        #rho = X[1]
        SO = B[8]
        #SOM_2 = rho*((2/(1-rho))**(3))*(Qeff/Z)*B[0:3]
        #SOM_4 = rho*((2/(1-rho))**(5))*(Qeff/Z)*B[3:8]
        #SOM = np.array([])
        #SOM = np.append(SOM,SOM_2)
        #SOM = np.append(SOM,SOM_4)
        SO_matrix = Cr3_oct.SO(ion,L,S,SO)
        H = B[0]*O20 + B[1]*O21 + B[2]*O22 + B[3]*O40 + B[4]*O41 + B[5]*O42 + B[6]*O43 + B[7]*O44 + SO_matrix
        E = np.linalg.eigh(H)
        E_val = E[0]
        #E_state = E[1].T
        E_excitation = E_val-E_val[0]
        E_calc = np.array([E_excitation[4], E_excitation[8], E_excitation[12], E_excitation[16], E_excitation[20], E_excitation[24]])
    
        return np.sum((E_calc/E_exp - 1)**2)

    ### rho is in between 0.1 - 0.3 while Qeff is in between 0 - 1

result = scipy.optimize.minimize(Fit_B, B, method = 'Powell', args = (X,E_exp),tol=1e-6)
B = result.x
print(result)
print("New B: " + str(B))

"""



