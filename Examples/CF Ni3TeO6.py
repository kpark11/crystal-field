# -*- coding: utf-8 -*-
"""
Created on Thu Jun  3 15:14:24 2021

@author: Student
"""

import numpy as np
import Crystal_Field_Calculations as cef

Z = 2 # Oxygens
ion = "Ni2"
L = 3
S = 1

Ni = cef.LS(L,S)



Ni1 =  np.array([-2.554350,    4.424264,    6.967479])
O11 =  np.array([-1.848072,    2.980332,    8.335001])
O12 = np.array([-3.319889,    2.772391,    6.121267])
O13 =  np.array([-4.157971,    4.534575,    8.335001])
O14 =  np.array([-3.602144,    5.913176,    6.121267])
O15 =  np.array([-1.657007,    5.757884,    8.335001])
O16 =  np.array([-0.741017,    4.587224,    6.121267])

O11_d = O11 - Ni1
O12_d = O12 - Ni1
O13_d = O13 - Ni1
O14_d = O14 - Ni1
O15_d = O15 - Ni1
O16_d = O16 - Ni1

d1 = np.array([O11_d,O12_d,O13_d,O14_d,O15_d,O16_d])

B1 = Ni.PC(ion,L,S,d1,Z)
O20 = Ni.Olm(L,S,2,0)
O40 = Ni.Olm(L,S,4,0)
O43 = Ni.Olm(L,S,4,3)

Hcf1 = B1[0]*O20 + B1[3]*O40 + B1[6]*O43
Ecf_val,Ecf_val_excitation,H_cf_vt = Ni.Diag(Hcf1,printfunction=True)
SO_matrix = Ni.SO(ion,L,S)
Hcf_so1 = Hcf1 + SO_matrix
Ecf_so_val,Ecf_so_val_excitation,H_cf_so_vt = Ni.Diag(Hcf_so1,printfunction=True)


 
Ni2 =  np.array([-2.554350,    4.424264,    9.744282])
O21 =  np.array([-0.765539,    4.247146,   10.710267])
O22 =  np.array([-1.848072,    2.980332,    8.335001])
O23 =  np.array([-4.157971,    4.534575,    8.335001])
O24 =  np.array([-1.657007,    5.757884,    8.335001])
O25 =  np.array([-3.602144,    2.963667,   10.710267])
O26 =  np.array([-3.295367,    6.061979,   10.710267])


O21_d = O21 - Ni2
O22_d = O22 - Ni2
O23_d = O23 - Ni2
O24_d = O24 - Ni2
O25_d = O25 - Ni2
O26_d = O26 - Ni2


d2 = np.array([O21_d,O22_d,O23_d,O24_d,O25_d,O26_d])

B2 = Ni.PC(ion,L,S,d2,Z)
Hcf2 = B2[0]*O20 + B2[3]*O40 + B2[6]*O43
Ecf_val,Ecf_val_excitation,H_cf_vt = Ni.Diag(Hcf2,printfunction=True)
SO_matrix = Ni.SO(ion,L,S)
Hcf_so2 = Hcf2 + SO_matrix
Ecf_so_val,Ecf_so_val_excitation,H_cf_so_vt = Ni.Diag(Hcf_so2,printfunction=True)



Ni3  = np.array([-2.554350,    4.424264,    0.000000])
O31  = np.array([-4.402422,    4.455086,   -0.842999])
O32  = np.array([-1.657007,    2.808375,   -0.842999])
O33  = np.array([-3.295367,    3.112470,    1.532267])
O34  = np.array([-3.319889,    5.721900,    1.532267])
O35  = np.array([-1.047794,    4.438421,    1.532267])
O36  = np.array([-1.603621,    6.009330,   -0.842999])


O31_d = O31 - Ni3
O32_d = O32 - Ni3
O33_d = O33 - Ni3
O34_d = O34 - Ni3
O35_d = O35 - Ni3
O36_d = O36 - Ni3

d3 = np.array([O31_d,O32_d,O33_d,O34_d,O35_d,O36_d])


B3 = Ni.PC(ion,L,S,d3,Z)
Hcf3 = B3[0]*O20 + B3[3]*O40 + B3[6]*O43
Ecf_val,Ecf_val_excitation,H_cf_vt = Ni.Diag(Hcf3,printfunction=True)
SO_matrix = Ni.SO(ion,L,S)
Hcf_so3 = Hcf3 + SO_matrix
Ecf_so_val,Ecf_so_val_excitation,H_cf_so_vt = Ni.Diag(Hcf_so3,printfunction=True)
