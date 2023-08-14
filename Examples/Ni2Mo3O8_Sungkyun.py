# -*- coding: utf-8 -*-
"""
Created on Fri Feb 17 11:26:31 2023

@author: brian
"""
import numpy as np
import Crystal_Field_Calculations as cef
import matplotlib.pyplot as plt
from scipy import optimize


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


B02_oct = B_oct[0]
B04_oct = B_oct[3]
B34_oct = B_oct[6]



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

SO_oct = -42


SO_matrix_oct = Ni2_oct.SO(ion,L,S,SO_oct)

Hcf_so_oct = Hcf_oct + SO_matrix_oct

Ecf_so_val_oct,Ecf_so_val_excitation_oct,H_cf_so_vt_oct = Ni2_oct.Diag(Hcf_so_oct,printfunction=True)

#################################### Spin-Spin interaction ##########################################
print('\n')
print('---------------- Crystal Field with Spin-Orbit coupling and SS interaction ---------------')


p = 1.1


Hss = Ni2_oct.SS(ion,L,S,p,)
Hcf_so_ss_oct = Hcf_so_oct + Hss
Ecf_so_ss_val_oct,Ecf_so_ss_val_excitation_oct,H_cf_so_ss_vt_oct = Ni2_oct.Diag(Hcf_so_ss_oct,printfunction=True)




#################################### Molecular Field ##########################################


print('---------------- Crystal Field with Spin-Orbit coupling and Molecular Field ---------------')

Hm=120

Hg_oct = Ni2_oct.Molecular_Field_Sz(ion,L,S,Hm)

Hcf_so_ss_g_oct = Hcf_so_ss_oct + Hg_oct

Ecf_so_ss_g_val_oct,Ecf_so_ss_g_val_excitation_oct,H_cf_so_ss_g_vt_oct = Ni2_oct.Diag(Hcf_so_ss_g_oct,printfunction=True)


#################################### Magnetic Moment ##########################################

print('\n')

print('---------- Magnetic moment----------')

MagMoment_L_z,MagMoment_S_z,MagMoment_M_z = Ni2_oct.MagMoment_z(L,S,Hcf_so_ss_oct)

MagMoment_L_x,MagMoment_S_x,MagMoment_M_x = Ni2_oct.MagMoment_x(L,S,Hcf_so_ss_oct)

MagMoment_L_y,MagMoment_S_y,MagMoment_M_y = Ni2_oct.MagMoment_y(L,S,Hcf_so_ss_oct)




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

B02_tet = B_tet[0]
#B02_tet  = B02
B04_tet = B_tet[3]
#B04_tet = B04
B34_tet = B_tet[6]
#B34_tet = B34

#################################### Crystal Field ##########################################

print('---------------- Crystal Field ---------------')




Hcf_tet = B02_tet*O20 + B04_tet*O40 + B34_tet*O43
Ecf_val_tet,Ecf_val_excitation_tet,H_cf_vt_tet = Ni2_tet.Diag(Hcf_tet,printfunction=True)





#################################### Spin-Orbit coupling ##########################################
print('---------------- Crystal Field with Spin-Orbit coupling ---------------')


SO_tet = -42

SO_matrix_tet = Ni2_tet.SO(ion,L,S,SO_tet)

Hcf_so_tet = Hcf_tet + SO_matrix_tet

Ecf_so_val_tet,Ecf_so_val_excitation_tet,H_cf_so_vt_tet = Ni2_tet.Diag(Hcf_so_tet,printfunction=True)

#################################### Spin-Spin interaction ##########################################
print('\n')
print('---------------- Crystal Field with Spin-Orbit coupling and SS interaction ---------------')

p = 1.1


Hss = Ni2_tet.SS(ion,L,S,p)
Hcf_so_ss_tet = Hcf_so_tet + Hss
Ecf_so_ss_val_tet,Ecf_so_ss_val_excitation_tet,H_cf_so_ss_vt_tet = Ni2_tet.Diag(Hcf_so_ss_tet,printfunction=True)


#################################### Molecular Field ##########################################


print('---------------- Crystal Field with Spin-Orbit coupling and Molecular Field ---------------')


Hm=120

Hg_tet = Ni2_tet.Molecular_Field_Sz(ion,L,S,Hm)

Hcf_so_ss_g_tet = Hcf_so_ss_tet + Hg_tet

Ecf_so_ss_g_val_tet,Ecf_so_ss_g_val_excitation_tet,H_cf_so_ss_g_vt_tet = Ni2_tet.Diag(Hcf_so_ss_g_tet,printfunction=True)




#################################### Magnetic Moment ##########################################

print('\n')

print('---------- Magnetic moment----------')

MagMoment_L_z,MagMoment_S_z,MagMoment_M_z = Ni2_tet.MagMoment_z(L,S,Hcf_so_ss_tet)

MagMoment_L_x,MagMoment_S_x,MagMoment_M_x = Ni2_tet.MagMoment_x(L,S,Hcf_so_ss_tet)

MagMoment_L_y,MagMoment_S_y,MagMoment_M_y = Ni2_tet.MagMoment_y(L,S,Hcf_so_ss_tet)



