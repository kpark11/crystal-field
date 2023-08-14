# -*- coding: utf-8 -*-
"""
Created on Tue Jun  8 16:01:38 2021

@author: Kiman Park
"""

import numpy as np
import crystalfield_test as cef
import matplotlib.pyplot as plt


######## Literature Yang et al, Physica B 370 (2005). #######
#  B20 = v - 2*sqrt(2)*v'
#  B40 = (4/3)*v + 2*sqrt(2)*v' - 14*Dq
#  B43 = -(7/sqrt(70))*((2/3)*v + sqrt(2)*v' + 20*Dq)



np.set_printoptions(suppress=True)

ion = 'Mn2'
L = 0
S = 5/2
Z = 2
Mn = cef.LS(L,S)

Te11  =  np.array([4.701347,    3.987663,    5.304571])
Te12  =  np.array([4.617244,    0.021975,    5.304571])
Te13  =  np.array([1.224909,    2.077654,    5.304571])
Te14  =  np.array([2.289591,    4.009639,    8.950430])
Mn1  =  np.array([3.514500,    2.029098,    7.137193])
Te15  =  np.array([2.411755,   -0.021975,    8.950430])
Te16  =  np.array([5.842153,    2.099629,    8.950430])

Te11d = Te11 - Mn1
Te12d = Te12 - Mn1
Te13d = Te13 - Mn1
Te14d = Te14 - Mn1
Te15d = Te15 - Mn1
Te16d = Te16 - Mn1

Mn1d = np.array([Te11d,Te12d,Te13d,Te14d,Te15d,Te16d])

B1 = Mn.PC(ion,L,S,Mn1d,Z)


Te1_11   =  np.array([-1.186847,    2.099629,    8.950430])
Te1_12   =  np.array([-1.102745,    6.065317,    8.950430])
Te1_13    =  np.array([1.224909,    2.077654,    5.304571])
Te1_14    =  np.array([2.289591,    4.009639,    8.950430])
Mn1_1   =  np.array([-0.000000,    4.058195,    7.117807])
Te1_15   =  np.array([-2.327653,    3.987663,    5.304571])
Te1_16   =  np.array([ 1.102745,    6.109268,    5.304571])



Te21    =  np.array([1.186847,    2.099629,   12.432071])
Te22    =  np.array([2.289591,    4.009639,    8.950430])
Mn2    =  np.array([3.514500,    2.029097,   10.691250])
Te23    =  np.array([4.617244,   -0.021975,   12.432071])
Te24    =  np.array([4.739409,    4.009639,   12.432071])
Te25    =  np.array([2.411755,   -0.021975,    8.950430])
Te26    =  np.array([5.842153,    2.099629,    8.950430])






