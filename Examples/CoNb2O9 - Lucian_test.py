# -*- coding: utf-8 -*-
"""
Created on Wed March 24 
@author: chemuser
"""
import Crystal_Field_Calcuations as cef
import numpy as np
import scipy.optimize


np.set_printoptions(suppress=True)

ion = 'Co2'
L = 3
S = 3/2
Z = 2
Co2_oct = cef.LS(L,S)

###################### Co1 ##############################
Co_oct = np.array([0, 0, 0])

O1 = np.array([1.3469, 1.227, 0.9227])

O2 = np.array([-1.3469, 1.227, -0.9227])

O3 = np.array([1.0653, -1.697, 0.6982])

O4 = np.array([-1.0653, -1.697, -0.6982])

O5 = np.array([-1.0653, -0.3798, 1.8242])

O6 = np.array([1.0653, -0.3798, -1.8242])

O1_d = O1 - Co_oct
O2_d = O2 - Co_oct
O3_d = O3 - Co_oct
O4_d = O4 - Co_oct
O5_d = O5 - Co_oct
O6_d = O6 - Co_oct


d = np.array([O1_d,O2_d,O3_d,O4_d,O5_d,O6_d])

B = Co2_oct.PC(ion,L,S,d,Z)

#################################### Crystal Field ##########################################


O20 = Co2_oct.Olm(L,S,2,0)
O21 = Co2_oct.Olm(L,S,2,1)
O22 = Co2_oct.Olm(L,S,2,2)
O40 = Co2_oct.Olm(L,S,4,0)
O41 = Co2_oct.Olm(L,S,4,1)
O42 = Co2_oct.Olm(L,S,4,2)
O43 = Co2_oct.Olm(L,S,4,3)
O44 = Co2_oct.Olm(L,S,4,4)


print('---------------- Crystal Field cubi ---------------')


"""
X = np.array([0.09991105, 0.21181041])
rho = X[1]
Qeff = X[0]
SOM = Co2_oct.SOM(ion,L,S,d,Z,rho,Qeff)
"""

Hcf = B[0]*O20 + B[1]*O21 + B[2]*O22 + B[3]*O40 + B[4]*O41 + B[5]*O42 + B[6]*O43 + B[7]*O44
Ecf_val,Ecf_val_excitation,H_cf_vt = Co2_oct.Diag(Hcf,printfunction=True)




#################################### Spin-Orbit coupling ##########################################
print('---------------- Crystal Field with Spin-Orbit coupling ---------------')


SO = -30 #meV

#SO = B[8]

SO_matrix = Co2_oct.SO(ion,L,S,SO)

Hcf_so = Hcf + SO_matrix

Ecf_so_val,Ecf_so_val_excitation,H_cf_so_vt = Co2_oct.Diag(Hcf_so,printfunction=True)




#################################### Magnetic Moment ##########################################
print('--------- Magnetic Properties --------')




MagMoment_L_z,MagMoment_S_z,MagMoment_M_z = Co2_oct.MagMoment_z(L,S,Hcf_so)

MagMoment_L_x,MagMoment_S_x,MagMoment_M_x = Co2_oct.MagMoment_x(L,S,Hcf_so)

MagMoment_L_y,MagMoment_S_y,MagMoment_M_y = Co2_oct.MagMoment_y(L,S,Hcf_so)


################### Testing ##################


E_exp = np.array([30.0,30.0,30.0,30.0, 50.0,50.0,50.0,50.0, 107.0,107.0,107.0,107.0, 136.0,136.0,136.0,136.0, 170.0,170.0,170.0,170.0, 737.0,737.0,737.0,737.0])

"""

def Fit_B(B,E_exp):
    Eexp = E_exp
    #Qeff = X[0]
    #rho = X[1]
    #SOM_2 = rho*((2/(1-rho))**(3))*(Qeff/Z)*B[0:3]
    #SOM_4 = rho*((2/(1-rho))**(5))*(Qeff/Z)*B[3:8]
    #SOM = np.array([])
    #SOM = np.append(SOM,SOM_2)
    #SOM = np.append(SOM,SOM_4)
    H = B[0]*O20 + B[1]*O21 + B[2]*O22 + B[3]*O40 + B[4]*O41 + B[5]*O42 + B[6]*O43 + B[7]*O44
    Ecf_val,Ecf_val_excitation,H_cf_vt = Co2_oct.Diag(H)
    E_calc = np.array([Ecf_so_val_excitation[4:]])
    
    return (np.sum(Eexp - E_calc))**2

### rho is in between 0.1 - 0.3 while Qeff is in between 0 - 1

result = scipy.optimize.minimize(Fit_B, B, method = 'Powell', args = (E_exp),tol=1e-16)
B = result.x
print(result)




def Fit_SOM(X,B,E_exp):
    Qeff,rho = X
    SOM_2 = rho*((2/(1-rho))**(3))*(Qeff/Z)*B[0:3]
    SOM_4 = rho*((2/(1-rho))**(5))*(Qeff/Z)*B[3:8]
    SOM = np.array([])
    SOM = np.append(SOM,SOM_2)
    SOM = np.append(SOM,SOM_4)
    SO_matrix = Co2_oct.SO(ion,L,S)
    H = SOM[0]*O20 + SOM[1]*O21 + SOM[2]*O22 + SOM[3]*O40 + SOM[4]*O41 + SOM[5]*O42 + SOM[6]*O43 + SOM[7]*O44 + SO_matrix
    Ecf_val,Ecf_val_excitation,H_cf_vt = Co2_oct.Diag(H)

    E_calc = np.array([Ecf_so_val_excitation[4:]])
    
    return np.sum(E_calc - E_exp)**2

### rho is in between 0.1 - 0.3 while Qeff is in between 0 - 1
    
bnds1 = (0,1)
bnds2 = (0,1)
bounds = [bnds1,bnds2]

result = scipy.optimize.minimize(Fit_SOM, X, method = 'Powell', args = (B,E_exp),bounds=bounds,tol=1e-16)
X = result.x
print(result)

"""