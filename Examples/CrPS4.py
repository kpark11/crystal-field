# -*- coding: utf-8 -*-
"""
Created on Tue Dec 14 13:57:19 2021

@author: Student
"""

import Crystal_Field_Calculations as cef
import numpy as np
import matplotlib.pyplot as plt
from scipy import optimize
###################################################################################
np.set_printoptions(suppress=True)

ion = 'Cr3'
L = 3
S = 3/2
Z = 2
Cr3_oct = cef.LS(L,S)
"""
Cr  =  np.array([5.231547,    5.482989,    6.099909])
S1  =  np.array([6.749071,    5.561359,    4.267496])
S2  =  np.array([4.114306,    3.628235,    4.928116])
S3  =  np.array([3.842206,    7.256470,    5.199562])
S4  =  np.array([3.714024,    5.561359,    7.932321])
S5  =  np.array([6.348789,    3.628235,    7.271701])
S6  =  np.array([6.620889,    7.256470,    7.000255])
"""



### HS Kim's cif file ####
Cr =  np.array([-0.201430,    3.695913,    6.136695])
S1 =  np.array([1.323463,    3.806899,    4.291391])
S2 =  np.array([-1.316974,    1.881688,    4.949858])
S3 =  np.array([-1.726324,    3.806899,    7.981999])
S4 =  np.array([-1.588674,    5.494905,    5.245647])
S5 =  np.array([0.914113,    1.881688,    7.323532])
S6 =  np.array([1.185813,    5.494905,    7.027743])




S1_d = S1 - Cr
S2_d = S2 - Cr
S3_d = S3 - Cr
S4_d = S4 - Cr
S5_d = S5 - Cr
S6_d = S6 - Cr

d = np.array([S1_d,S2_d,S3_d,S4_d,S5_d,S6_d])



########################################### Crystal Field Splitting ########################################


B = Cr3_oct.PC(ion,L,S,d,Z)


O20 = Cr3_oct.Olm(L,S,2,0)
O21 = Cr3_oct.Olm(L,S,2,1)
O22 = Cr3_oct.Olm(L,S,2,2)
O40 = Cr3_oct.Olm(L,S,4,0)
O41 = Cr3_oct.Olm(L,S,4,1)
O42 = Cr3_oct.Olm(L,S,4,2)
O43 = Cr3_oct.Olm(L,S,4,3)
O44 = Cr3_oct.Olm(L,S,4,4)




Hcf = B[0]*O20 + B[1]*O21 + B[2]*O22 + B[3]*O40 + B[4]*O41 + B[5]*O42 + B[6]*O43 + B[7]*O44
#Hcf = B[0]*O20 + B[2]*O22 + B[3]*O40 + B[5]*O42 + B[7]*O44





Ecf_val,Ecf_val_excitation,H_cf_vt = Cr3_oct.Diag(Hcf,printfunction=True)


####################################### Spin-Orbit Coupling #######################################

SO = 11.28255 #meV

SO_matrix = Cr3_oct.SO(ion,L,S)

Hcf_so = Hcf + SO_matrix

Ecf_so_val,Ecf_so_val_excitation,H_cf_so_vt = Cr3_oct.Diag(Hcf_so,printfunction=True)



#################################### Magnetic Moment ##########################################
print('--------- Magnetic Properties --------')




MagMoment_L_z,MagMoment_S_z,MagMoment_M_z = Cr3_oct.MagMoment_z(L,S,Hcf_so)

MagMoment_L_x,MagMoment_S_x,MagMoment_M_x = Cr3_oct.MagMoment_x(L,S,Hcf_so)

MagMoment_L_y,MagMoment_S_y,MagMoment_M_y = Cr3_oct.MagMoment_y(L,S,Hcf_so)



################## specific heat ###############

C_schottky,T = Cr3_oct.SpecificHeat(Ecf_val_excitation)
C_schottky,T = Cr3_oct.SpecificHeat(Ecf_so_val_excitation)


################## magnetic entropy ###############

entropy,T = Cr3_oct.Entropy(Ecf_val_excitation)
entropy,T = Cr3_oct.Entropy(Ecf_so_val_excitation)


############# Zeeman effect ###############
"""
Fields = 25 # T


Zeeman = Cr3_oct.Zeeman_z(L,S,Hcf_so,Fields,savefunction=False)
Zeeman = Cr3_oct.Zeeman_x(L,S,Hcf_so,Fields,savefunction=False)
Zeeman = Cr3_oct.Zeeman_y(L,S,Hcf_so,Fields,savefunction=False)

"""
############# Magnetic susceptibility ###############

X_z,T = Cr3_oct.Susceptibility_z(L,S,Hcf_so,SO)
X_x,T = Cr3_oct.Susceptibility_x(L,S,Hcf_so,SO)
X_y,T = Cr3_oct.Susceptibility_y(L,S,Hcf_so,SO)



plt.plot(T,X_z,label = r'z')
plt.plot(T,X_x,label = r'x')
plt.plot(T,X_y,label = r'y')
plt.xlabel('Temperature (K)')
plt.ylabel('Magnetic Susceptibility ($\mu_B$)')
plt.title('Magnetic Susceptibility')
plt.legend()
plt.grid(alpha=.4,linestyle='--')
plt.xlim(0,300)
plt.show()


#exp_susceptibility = np.loadtxt(r'C:\Users\Student\OneDrive - University of Tennessee\Desktop\Research\MPS (M = Cr, Mn)\CrPS4\Digitized_susceptibility_c.dat')





Exp = np.array([0.000000001,0.000000001,0.000000001,0.000000001,1.49,1.49,1.49,1.49,1.51,1.51,1.51,1.51,1.74,1.74,1.74,1.74,2,2,2,2,3,3,3,3,4,4,4,4])*1000
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
        H = B[0]*O20 + B[1]*O21 + B[2]*O22 + B[3]*O40 + B[4]*O41 + B[5]*O42 + B[6]*O43 + B[7]*O44 #+ SO_matrix
        E = np.linalg.eigh(H)
        E_val = E[0]
        #E_state = E[1].T
        E_excitation = E_val-E_val[0]
        
        return np.sum((E_excitation[4:16]/E_exp[4:16] - 1)**2)

    ### rho is in between 0.1 - 0.3 while Qeff is in between 0 - 1
t = Fit_B(B,Exp)



result = optimize.minimize(Fit_B, x0=B, args = (Exp,), method = 'Powell',tol=1e-20)
B = result.x
print(result)




Exp = np.array([0.000000001,0.000000001,0.000000001,0.000000001,1.2,1.2,1.2,1.2,1.48,1.48,1.48,1.48,1.51,1.51,1.51,1.51,1.75,1.75,1.75,1.75,2,2,2,2,2.76,2.76,2.76,2.76])*1000
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
        H = B[0]*O20 + B[1]*O21 + B[2]*O22 + B[3]*O40 + B[4]*O41 + B[5]*O42 + B[6]*O43 + B[7]*O44 #+ SO_matrix
        E = np.linalg.eigh(H)
        E_val = E[0]
        #E_state = E[1].T
        E_excitation = E_val-E_val[0]
        
        return np.sum((E_excitation[4:]/E_exp[4:] - 1)**2)

    ### rho is in between 0.1 - 0.3 while Qeff is in between 0 - 1
t = Fit_B(B,Exp)



result = optimize.minimize(Fit_B, x0=B, args = (Exp,), method = 'Powell',tol=1e-20)
B = result.x
print(result)



Hcf = B[0]*O20 + B[1]*O21 + B[2]*O22 + B[3]*O40 + B[4]*O41 + B[5]*O42 + B[6]*O43 + B[7]*O44
#Hcf = B[0]*O20 + B[2]*O22 + B[3]*O40 + B[5]*O42 + B[7]*O44





Ecf_val,Ecf_val_excitation,H_cf_vt = Cr3_oct.Diag(Hcf,printfunction=True)


####################################### Spin-Orbit Coupling #######################################

SO = 11.28255 #meV

SO_matrix = Cr3_oct.SO(ion,L,S)

Hcf_so = Hcf + SO_matrix

Ecf_so_val,Ecf_so_val_excitation,H_cf_so_vt = Cr3_oct.Diag(Hcf_so,printfunction=True)



