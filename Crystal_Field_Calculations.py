# -*- coding: utf-8 -*-
"""
Created on Sun Sep  6 15:57:45 2020

@author: brian
"""
import numpy as np
import matplotlib.pyplot as plt
#np.set_printoptions(suppress=True)

# 1st rule of Hund: The combination of electron spins, s,
# that gives the lowest energy is that with the highest value of (2s+1)

# 2nd rule of Hund: When the first rule has been satisfied there are 
# in general several possible l values having the same value of (2s+1), 
# that with the largest l will be the most stable.

# 3rd rule of Hund: The most stable configurations result from:
# j = | l-s |... less than half filled shells
# j = | l+s |... more than half filled shells
# j = s ... half filled shell (l = 0)



    
class LS():
    def __init__(self,L,S):
        self.L = L
        self.S = S
        if L == 0 or 1 or 2 or 3 or 4 or 5 or 6 or 7:
            self.L_array = np.arange(-L,L+1)
        else:
            L1 = int(L + 0.5)
            self.L_array = np.arange(-L1,L1) + 0.5
        if S == 0 or 1 or 2 or 3 or 4 or 5 or 6 or 7:
            self.S_array = np.arange(-S,S+1)
        else:
            S1 = int(S + 0.5)
            self.S_array = np.arange(-S1,S1) + 0.5
        self.L_degeneracy = int(2*L+1)
        self.S_degeneracy = int(2*S+1)       
        self.LS_degeneracy = self.L_degeneracy*self.S_degeneracy
        self.L_matrix = []
        self.S_matrix = []
        self.LS_matrix = []
        self.LS_matrix_1 = np.zeros((self.LS_degeneracy,self.LS_degeneracy))
        self.SO_matrix = np.zeros((self.LS_degeneracy,self.LS_degeneracy))
        
        for i in range(self.L_degeneracy):
            for k in range(self.S_degeneracy):
                self.L_matrix = np.append(self.L_matrix,self.L_array[i])
        for i in range(self.L_degeneracy):
            for k in range(self.S_degeneracy):
                self.S_matrix = np.append(self.S_matrix,self.S_array[k])
        self.L_tot = np.identity(self.LS_degeneracy)*L
        self.S_tot = np.identity(self.LS_degeneracy)*S
        self.i = np.identity(self.LS_degeneracy)

    
    def Lz(self,L,S):
        elements = LS(L,S)
        elements.LS_matrix = elements.L_matrix
        elements.LS_matrix_diag = np.diag(elements.LS_matrix)
        return elements.LS_matrix_diag
    
    def Lplus(self,L,S):
        elements = LS(L,S)
        Lplus = np.sqrt(elements.L*(elements.L+1)-elements.L_matrix*(elements.L_matrix+1))
        return Lplus    
    
    def Lp(self,L,S):
        elements = LS(L,S)
        Lplus = np.sqrt(elements.L*(elements.L+1)-elements.L_matrix*(elements.L_matrix+1))
        for i in range(elements.LS_degeneracy - elements.S_degeneracy):
            elements.LS_matrix_1[i][i+elements.S_degeneracy] = Lplus[i]
        return elements.LS_matrix_1
    
    def Lp2(self,L,S):
        elements = LS(L,S)
        Lp2 = np.dot(elements.Lp(L,S),elements.Lp(L,S))
        return Lp2
    
    def Lp3(self,L,S):
        elements = LS(L,S)
        Lp3 = np.dot(elements.Lp(L,S),elements.Lp2(L,S))
        return Lp3
    
    def Lp4(self,L,S):
        elements = LS(L,S)
        Lp4 = np.dot(elements.Lp(L,S),elements.Lp3(L,S))
        return Lp4
    
    def Lp5(self,L,S):
        elements = LS(L,S)
        Lp5 = np.dot(elements.Lp(L,S),elements.Lp4(L,S))
        return Lp5
    
    def Lp6(self,L,S):
        elements = LS(L,S)
        Lp6 = np.dot(elements.Lp(L,S),elements.Lp5(L,S))
        return Lp6
    
    def Lminus(self,L,S):
        elements = LS(L,S)
        Lminus = np.sqrt(elements.L*(elements.L+1)-elements.L_matrix*(elements.L_matrix-1))
        return Lminus
    
    def Lm(self,L,S):
        elements = LS(L,S)
        Lminus = np.sqrt(elements.L*(elements.L+1)-elements.L_matrix*(elements.L_matrix-1))
        for i in range(elements.LS_degeneracy - elements.S_degeneracy):
            elements.LS_matrix_1[i+elements.S_degeneracy][i] = Lminus[i + elements.S_degeneracy]
        return elements.LS_matrix_1
    
    def Lm2(self,L,S):
        elements = LS(L,S)
        Lm2 = np.dot(elements.Lm(L,S),elements.Lm(L,S))
        return Lm2
    
    def Lm3(self,L,S):
        elements = LS(L,S)
        Lm3 = np.dot(elements.Lm(L,S),elements.Lm2(L,S))
        return Lm3
    
    def Lm4(self,L,S):
        elements = LS(L,S)
        Lm4 = np.dot(elements.Lm(L,S),elements.Lm3(L,S))
        return Lm4
        
    def Lm5(self,L,S):
        elements = LS(L,S)
        Lm5 = np.dot(elements.Lm(L,S),elements.Lm4(L,S))
        return Lm5
        
    def Lm6(self,L,S):
        elements = LS(L,S)
        Lm6 = np.dot(elements.Lm(L,S),elements.Lm5(L,S))
        return Lm6

    def Sz(self,L,S):
        elements = LS(L,S)
        elements.LS_matrix = elements.S_matrix
        elements.LS_matrix = np.diag(elements.LS_matrix)
        return elements.LS_matrix
    
    def Splus(self,L,S):
        elements = LS(L,S)
        Splus = np.sqrt(elements.S*(elements.S+1)-elements.S_matrix*(elements.S_matrix+1))
        return Splus   
    
    def Sp(self,L,S):
        elements = LS(L,S)
        Splus = np.sqrt(elements.S*(elements.S+1)-elements.S_matrix*(elements.S_matrix+1))
        for i in range(elements.LS_degeneracy - 1):
            elements.LS_matrix_1[i][i+1] = Splus[i]
        return elements.LS_matrix_1
    
    def Sminus(self,L,S):
        elements = LS(L,S)
        Sminus = np.sqrt(elements.S*(elements.S+1)-elements.L_matrix*(elements.S_matrix-1))
        return Sminus
    
    def Sm(self,L,S):
        elements = LS(L,S)
        Sminus = np.sqrt(elements.S*(elements.S+1)-elements.S_matrix*(elements.S_matrix-1))
        for i in range(elements.LS_degeneracy - 1):
            elements.LS_matrix_1[i+1][i] = Sminus[i + 1]
        return elements.LS_matrix_1
    
    def Lx(self,L,S):
        elements = LS(L,S)
        Lx_matrix = 0.5*(elements.Lp(L,S) + elements.Lm(L,S))
        return Lx_matrix
    
    def Ly(self,L,S):
        elements = LS(L,S)
        Ly_matrix = 0.5j*(elements.Lp(L,S) - elements.Lm(L,S))
        return Ly_matrix
    
    def Sx(self,L,S):
        elements = LS(L,S)
        Sx_matrix = 0.5*(elements.Sp(L,S) + elements.Sm(L,S))
        return Sx_matrix
    
    def Sy(self,L,S):
        elements = LS(L,S)
        Sy_matrix = 0.5j*(elements.Sp(L,S) - elements.Sm(L,S))
        return Sy_matrix
    
    def Olm(self,L,S,l,m):
        
        element = LS(L,S)
        Lz = element.Lz(L,S)
        Lp = element.Lp(L,S)
        Lp2 = element.Lp2(L,S)
        Lp3 = element.Lp3(L,S)
        Lp4 = element.Lp4(L,S)
        Lp5 = element.Lp5(L,S)
        Lp6 = element.Lp6(L,S)
        Lm = element.Lm(L,S)
        Lm2 = element.Lm2(L,S)
        Lm3 = element.Lm3(L,S)
        Lm4 = element.Lm4(L,S)
        Lm5 = element.Lm5(L,S)
        Lm6 = element.Lm6(L,S)
        L_matrix = element.L_tot
        S_matrix = element.S_tot
        i = element.i
        X = L_matrix*(L_matrix+1*i)
        
        if [l,m] == [0,0]:
           return np.zeros((int(2*L_matrix+i), int(2*S_matrix+i)))
        elif [l,m] == [1,0]:
           Olm = Lz
        elif [l,m] == [1,1]:
           Olm = 0.5 *(Lp + Lm)
        elif [l,m] == [1,-1]:
           Olm = -0.5j *(Lp - Lm)
           
        elif [l,m] == [2,2]:
           Olm = 0.5 *(Lp2 + Lm2)
        elif [l,m] == [2,1]:
           Olm = 0.25*(np.dot(Lz,(Lp + Lm)) + np.dot((Lp + Lm),Lz))
        elif [l,m] == [2,0]:
           Olm = 3*Lz**2 - X*i
        elif [l,m] == [2,-1]:
           Olm = -0.25j*(np.dot(Lz,(Lp - Lm)) + np.dot((Lp - Lm),Lz))
        elif [l,m] == [2,-2]:
           Olm = -0.5j *(Lp2 - Lm2)
           
        elif [l,m] == [3,3]:
           Olm = 0.5 *(Lp3 + Lm3)
        elif [l,m] == [3,2]:
           Olm = 0.25 *(np.dot((Lp2 + Lm2),Lz) + np.dot(Lz,(Lp2 + Lm2)))
        elif [l,m] == [3,1]:
           Olm = 0.25*(np.dot((Lp + Lm),(5*Lz**2 - X - 0.5*i)) + np.dot((5*Lz**2 - X - 0.5*i),(Lp + Lm)))
        elif [l,m] == [3,0]:
           Olm = 5*Lz**3 - (3*X-1)*Lz
        elif [l,m] == [3,-1]:
           Olm = -0.25j*(np.dot((Lp - Lm),(5*Lz**2 - X - 0.5*i)) + np.dot((5*Lz**2 - X - 0.5*i),(Lp - Lm)))
        elif [l,m] == [3,-2]:
           Olm = -0.25j*(np.dot(Lz,(Lp2 - Lm2)) + np.dot((Lp2 - Lm2),Lz))
        elif [l,m] == [3,-3]:
           Olm = -0.5j *(Lp3 - Lm3)
    
        elif [l,m] == [4,4]:
           Olm = 0.5 *(Lp4 + Lm4)
        elif [l,m] == [4,3]:
           Olm = 0.25 *(np.dot((Lp3 + Lm3),Lz) + np.dot(Lz,(Lp3 + Lm3)))
        elif [l,m] == [4,2]:
           Olm = 0.25 *(np.dot((Lp2 + Lm2),(7*Lz**2 -X -5*i)) + np.dot((7*Lz**2 -X -5*i),(Lp2 + Lm2)))
        elif [l,m] == [4,1]:
           Olm = 0.25 *(np.dot((Lp + Lm),(7*Lz**3 -(3*X+1)*Lz)) + np.dot((7*Lz**3 -(3*X+1)*Lz),(Lp + Lm)))
        elif [l,m] == [4,0]:
           Olm = 35*Lz**4 - (30*X -25*i)*Lz**2 + 3*X**2 - 6*X
        elif [l,m] == [4,-4]:
           Olm = -0.5j *(Lp**4 - Lm**4)
        elif [l,m] == [4,-3]:
           Olm = -0.25j *(np.dot((Lp3 - Lm3),Lz) + np.dot(Lz,(Lp3 - Lm3)))
        elif [l,m] == [4,-2]:
           Olm = -0.25j *(np.dot((Lp2 - Lm2),(7*Lz**2 -X -5*i)) + np.dot((7*Lz**2 -X -5*i),(Lp2 - Lm2)))
        elif [l,m] == [4,-1]:
           Olm = -0.25j *(np.dot((Lp - Lm),(7*Lz**3 -(3*X+1*i)*Lz)) + np.dot((7*Lz**3 -(3*X+1*i)*Lz),(Lp - Lm)))
        
        elif [l,m] == [5,5]:
           Olm = (1/2)*(Lp5 + Lm5)
        elif [l,m] == [5,4]:
           Olm = (1/4)*(np.dot((Lp4 + Lm4),Lz) + np.dot(Lz,(Lp4 + Lm4)))
        elif [l,m] == [5,3]:
           Olm = (1/4)*(np.dot((Lp3 + Lm3),(9*Lz**2 - X - 33*i/2)) + np.dot((9*Lz**2 - X - 33*i/2),(Lp3 + Lm3)))
        elif [l,m] == [5,2]:
           Olm = (1/4)*(np.dot((Lp2 + Lm2),(3*Lz**3 - (X - 6*i)*Lz)) + np.dot((3*Lz**3 - (X - 6*i)*Lz),(Lp2 + Lm2)))
        elif [l,m] == [5,1]:
           Olm = (1/4)*(np.dot((Lp + Lm),(21*Lz**4 - 14*Lz**2*X + X**2 - X + 3*i/2)) + np.dot((21*Lz**4 - 14*Lz**2*X + X**2 - X + 3*i/2),(Lp + Lm)))
        elif [l,m] == [5,0]:
           Olm = (63*Lz**5 - (70*X - 105*i)*Lz**3 + (15*X**2 - 50*X + 12*i)*Lz)
        elif [l,m] == [5,-1]:
           Olm = (-1j/4)*(np.dot((Lp + Lm),(21*Lz**4 - 14*Lz**2*X + X**2 - X + 3*i/2)) + np.dot((21*Lz**4 - 14*Lz**2*X + X**2 - X + 3*i/2),(Lp + Lm)))
        elif [l,m] == [5,-2]:
           Olm = (-1j/4)*(np.dot((Lp2 + Lm2),(3*Lz**3 - (X - 6*i)*Lz)) + np.dot((3*Lz**3 - (X - 6*i)*Lz),(Lp2 + Lm2)))
        elif [l,m] == [5,-3]:
           Olm = (-1j/4)*(np.dot((Lp3 + Lm3),(9*Lz**2 - X - 33*i/2)) + np.dot((9*Lz**2 - X - 33*i/2),(Lp3 + Lm3))) 
        elif [l,m] == [5,-4]:
           Olm = (-1j/4)*(np.dot((Lp4 + Lm4),Lz) + np.dot(Lz,(Lp4 + Lm4))) 
        elif [l,m] == [5,-5]:
           Olm = (-1j/2)*(Lp5 + Lm5) 
        
        elif [l,m] == [6,6]:
           Olm = 0.5*(Lp6 + Lm6)
        elif [l,m] == [6,5]:
           Olm = 0.25*(np.dot((Lp5 + Lm5),Lz) + np.dot(Lz,(Lp5 + Lm5)))
        elif [l,m] == [6,4]:
           Olm = 0.25*(np.dot((Lp4 + Lm4),(11*Lz**2 -X -38*i)) + np.dot((11*Lz**2 -X -38*i),(Lp4 + Lm4)))
        elif [l,m] == [6,3]:
           Olm = 0.25*(np.dot((Lp3 + Lm3),(11*Lz**3 -(3*X+59*i)*Lz)) + np.dot((11*Lz**3 -(3*X+59*i)*Lz),(Lp3 + Lm3)))
        elif [l,m] == [6,2]:
           Olm = 0.25*(np.dot((Lp2 + Lm2),(33*Lz**4 -(18*X+123*i)*Lz**2 +X**2 +10*X +102*i)) + np.dot((33*Lz**4 -(18*X+123*i)*Lz**2 +X**2 +10*X +102*i),(Lp2 + Lm2)))
        elif [l,m] == [6,1]:
           Olm = 0.25*(np.dot((Lp +Lm),((33*Lz**5 -(30*X-15*i)*Lz**3) +(5*X**2 -10*X +12*i)*Lz)) + np.dot((33*Lz**5 -(30*X-15*i)*Lz**3 +(5*X**2 -10*X +12*i)*Lz),(Lp+ Lm)))
        elif [l,m] == [6,0]:
           Olm = 231*Lz**6 - (315*X-735*i)*Lz**4 + (105*X**2 -525*X +294*i)*Lz**2 - 5*X**3 + 40*X**2 - 60*X
        elif [l,m] == [6,-6]:
           Olm = -0.5j*(Lp6 - Lm6)
        elif [l,m] == [6,-5]:
           Olm = -0.25j*(np.dot((Lp5 - Lm5),Lz) + np.dot(Lz,(Lp5 - Lm5)))
        elif [l,m] == [6,-4]:
           Olm = -0.25j*(np.dot((Lp4 - Lm4),(11*Lz**2 -X -38*i)) + np.dot((11*Lz**2 -X -38*i),(Lp4 - Lm4)))
        elif [l,m] == [6,-3]:
           Olm = -0.25j*(np.dot((Lp**3 - Lm**3),(11*Lz**3 -(3*X+59*i)*Lz)) + np.dot((11*Lz**3 -(3*X+59*i)*Lz),(Lp**3 - Lm**3)))
        elif [l,m] == [6,-2]:
           Olm = -0.25j*(np.dot((Lp**2 - Lm**2),(33*Lz**4 -(18*X+123*i)*Lz**2 +X**2 +10*X +102*i)) + np.dot((33*Lz**4 -(18*X+123*i)*Lz**2 +X**2 +10*X +102*i),(Lp**2 - Lm**2)))
        elif [l,m] == [6,-1]:
           Olm = -0.25j*(np.dot((Lp - Lm),(33*Lz**5 -(30*X-15*i)*Lz**3 +(5*X**2 -10*X +12*i)*Lz) + np.dot((33*Lz**5 -(30*X-15*i)*Lz**3 +(5*X**2 -10*X +12*i)*Lz),(Lp - Lm))))
        
        return Olm

    def Spherical_position(self,xyz_position):
        r = np.sqrt(xyz_position[:,0]**2 + xyz_position[:,1]**2 + xyz_position[:,2]**2)
        theta = np.arctan2(xyz_position[:,1],xyz_position[:,0])
        phi = np.arccos(xyz_position[:,2]/r)
        
        print('r:[], theta:[], phi:[]')
        print(r,theta,phi)
        # r = [:,0]
        # phi = [:,1]
        # theta = [:,2]
        
        sp_position = np.array([r,theta,phi])
        
        return sp_position
    
    def SO(self,ion,L,S,SO=None):
        elements = LS(L,S)
        if SO is None:
            SO = elements.r_l(ion)
            SO_val = SO[3]/8.065548 # converting it from cm-1 to meV
            print('SO coupling value (meV): ' + str(SO_val))
            Lz = elements.Lz(L,S)
            Sz = elements.Sz(L,S)
            Lp = elements.Lp(L,S)
            Lm = elements.Lm(L,S)
            Sp = elements.Sp(L,S)
            Sm = elements.Sm(L,S)
            SO_matrix = Lz*Sz + (1/2)*(np.dot(Lp,Sm) + np.dot(Lm,Sp))
            H_SO = SO_val*SO_matrix
        else:
            print('SO coupling value (meV): ' + str(SO))
            Lz = elements.Lz(L,S)
            Sz = elements.Sz(L,S)
            Lp = elements.Lp(L,S)
            Lm = elements.Lm(L,S)
            Sp = elements.Sp(L,S)
            Sm = elements.Sm(L,S)
            SO_matrix = Lz*Sz + (1/2)*(np.dot(Lp,Sm) + np.dot(Lm,Sp))
            H_SO = SO*SO_matrix
        return H_SO
    
    
    
    def SS(self,ion,L,S,p_val=None):
        elements = LS(L,S)
        if p_val is None:
            p = elements.r_l(ion)
            p_val = p[5]/8.065548 # meV
            print('SS coupling value (meV): ' + str(p_val))
            Lz = elements.Lz(L,S)
            Sz = elements.Sz(L,S)
            Lp = elements.Lp(L,S)
            Lm = elements.Lm(L,S)
            Sp = elements.Sp(L,S)
            Sm = elements.Sm(L,S)
            L_tot = elements.L_tot
            S_tot = elements.S_tot
            SO_matrix = Lz*Sz + (1/2)*(np.dot(Lp,Sm) + np.dot(Lm,Sp))
            H_SS_matrix = np.dot(SO_matrix,SO_matrix) + (1/2)*(SO_matrix) - (1/3)*(L_tot*(L_tot + self.i)*S_tot*(S_tot+ self.i))
            H_SS = -p_val*H_SS_matrix
        else:
            print('SS coupling value (meV): ' + str(p_val))
            Lz = elements.Lz(L,S)
            Sz = elements.Sz(L,S)
            Lp = elements.Lp(L,S)
            Lm = elements.Lm(L,S)
            Sp = elements.Sp(L,S)
            Sm = elements.Sm(L,S)
            L_tot = elements.L_tot
            S_tot = elements.S_tot
            SO_matrix = Lz*Sz + (1/2)*(np.dot(Lp,Sm) + np.dot(Lm,Sp))
            H_SS_matrix = np.dot(SO_matrix,SO_matrix) + (1/2)*(SO_matrix) - (1/3)*(L_tot*(L_tot + self.i)*S_tot*(S_tot+ self.i))
            H_SS = -p_val*H_SS_matrix
        return H_SS
    
            
    def Molecular_Field_Sz(self,ion,L,S,Hm=None):
        elements = LS(L,S)
        if Hm is None:
            Hm = elements.r_l(ion)
            Hm_val = Hm[5]/8.065548
            print('Molecular Field value (meV): ' + str(Hm_val))
            elements = LS(L,S)
            Sz = elements.Sz(L,S)
            g = 2.002319
            mu_B = 5.7883818012*(10**-2) #meV/T
            #mu_B = 5.7883818012*(10**-5) * 8065.548 #1/(cmT)
            #print('g: ' + str(g))
            Hg_z = g*Hm_val*mu_B*Sz
        else:
            print('Molecular Field value (meV): ' + str(Hm))
            Sz = elements.Sz(L,S)
            g = 2.002319
            mu_B = 5.7883818012*(10**-2) #meV/T
            #mu_B = 5.7883818012*(10**-5) * 8065.548 #1/(cmT)
            #print('g: ' + str(g))
            Hg_z = g*Hm*mu_B*Sz
        return Hg_z
    
    def Molecular_Field_Sx(self,ion,L,S,Hm=None):
        elements = LS(L,S)
        if Hm is None:
            Hm = elements.r_l(ion)
            Hm_val = Hm[5]/8.065548
            elements = LS(L,S)
            Sx = elements.Sx(L,S)
            g = 2.002319
            mu_B = 5.7883818012*(10**-2) #meV/T
            #mu_B = 5.7883818012e-5 * 8065.548 #1/(cmT)
            #print('g: ' + str(g))
            Hg_x = g*Hm_val*mu_B*Sx
        else:
            Sx = elements.Sx(L,S)
            g = 2.002319
            mu_B = 5.7883818012*(10**-2) #meV/T
            #mu_B = 5.7883818012e-5 * 8065.548 #1/(cmT)
            #print('g: ' + str(g))
            Hg_x = g*Hm*mu_B*Sx
        return Hg_x
    
    def Molecular_Field_Sy(self,ion,L,S,Hm=None):
        elements = LS(L,S)
        if Hm is None:
            Hm = elements.r_l(ion)
            Hm_val = Hm[5]/8.065548
            elements = LS(L,S)
            Sy = elements.Sy(L,S)
            g = 2.002319
            mu_B = 5.7883818012*(10**-2) #meV/T
            #mu_B = 5.7883818012e-5 * 8065.548 #1/(cmT)
            #print('g: ' + str(g))
            Hg_y = g*Hm_val*mu_B*Sy
        else:
            Sy = elements.Sy(L,S)
            g = 2.002319
            mu_B = 5.7883818012*(10**-2) #meV/T
            #mu_B = 5.7883818012e-5 * 8065.548 #1/(cmT)
            #print('g: ' + str(g))
            Hm
            Hg_y = g*Hm*mu_B*Sy
        return Hg_y
        
    
    def PC(self,ion,L,S,d,Z):
        
        
        e = 1.60217646e-19 # C
        e_0 = 8.854e-22 # C/VAngstrom
        a0 = 0.52917721067    #Bohr radius in \AA
        
        ### B_lm = -|e|*p_ml*gamma_lm*<r_l>*theta_l ###
        ### theta_2 = alpha ###
        ### theta_4 = beta ###
        ###theta_6 = gamma ###
        elements = LS(L,S)
        
        x=[]
        for k in range(len(d)):
            x = np.append(x,d[k][0])
        #print(x)
        y=[]
        for k in range(len(d)):
            y = np.append(y,d[k][1])
        #print(y)
        z=[]
        for k in range(len(d)):
            z = np.append(z,d[k][2])
        #print(z)
        
        d_sp = elements.Spherical_position(d)
        r = d_sp[0,:]
        
        print('---------------- Ligand positions (Angstrom) ---------------')
        print('[x[i]    y[i]    z[i]]:')
        print(d)
        print("r[i]: " + str(r))

        r_l = elements.r_l(ion)
        r_2 = r_l[1]*a0**2
        r_4 = r_l[2]*a0**4
        print(r_l)
        
        theta = elements.theta(ion)
        alpha = theta[0]
        beta = theta[1]
        
        #r_2 = cef.r_l(L,S)[]
        
        ### http://www2.cpfs.mpg.de/~rotter/homepage_mcphase/manual/node131.html ####
        Zs22 = (1/4)*np.sqrt(15/np.pi)*(2*x*y/r**5)
        Zs21 = (1/2)*np.sqrt(15/np.pi)*(y*z/r**5)
        Z20 = (1/4)*np.sqrt(5/np.pi)*(3*z**2 - r**2)/r**5
        Zc21 = (1/2)*np.sqrt(15/np.pi)*(x*z/r**5)
        Zc22 = (1/4)*np.sqrt(15/np.pi)*((x**2 - y**2)/r**5)

        Zs44 = (3/16)*np.sqrt(35/np.pi)*(4*(x**3*y - x*y**3)/r**9)
        Zs43 = (3/8)*np.sqrt(70/np.pi)*((3*x**2*y - y**3)*z/r**9)
        Zs42 = (3/8)*np.sqrt(5/np.pi)*(2*x*y*(7*z**2 - r**2)/r**9)
        Zs41 = (3/4)*np.sqrt(5/(2*np.pi))*(y*z*(7*z**2 - 3*r**2)/r**9)
        Z40 =  (3/16)*np.sqrt(1/np.pi)*((35*z**4 - 30*z**2*r**2 + 3*r**4)/r**9)
        Zc41 = (3/4)*np.sqrt(5/(2*np.pi))*(x*z*(7*z**2 - 3*r**2)/r**9)
        Zc42 = (3/8)*np.sqrt(5/np.pi)*((x**2 - y**2)*(7*z**2 - r**2)/r**9)
        Zc43 = (3/8)*np.sqrt(70/np.pi)*((x**3 - 3*x*y**2)*z/r**9)
        Zc44 = (3/16)*np.sqrt(35/np.pi)*((x**4 - 6*x**2*y**2 + y**4)/r**9)

        Bs22 = (-Z/5)*(1/4)*np.sqrt(15/np.pi)*r_2*alpha*np.sum(Zs22)*e/e_0*1000
        Bs21 = (-Z/5)*(1/4)*np.sqrt(15/np.pi)*r_2*alpha*np.sum(Zs21)*e/e_0*1000
        B20 = (-Z/5)*(1/4)*np.sqrt(5/np.pi)*r_2*alpha*np.sum(Z20)*e/e_0*1000
        Bc21 = (-Z/5)*(1/2)*np.sqrt(15/np.pi)*r_2*alpha*np.sum(Zc21)*e/e_0*1000
        Bc22 = (-Z/5)*(1/4)*np.sqrt(15/np.pi)*r_2*alpha*np.sum(Zc22)*e/e_0*1000

        B22 = Bs22 + Bc22
        B21 = Bs21 + Bc21
        
        Bs44 = (-Z/9)*(3/16)*np.sqrt(35/np.pi)*r_4*beta*np.sum(Zs44)*e/e_0*1000
        Bs43 = (-Z/9)*(3/8)*np.sqrt(70/np.pi)*r_4*beta*np.sum(Zs43)*e/e_0*1000
        Bs42 = (-Z/9)*(3/8)*np.sqrt(5/np.pi)*r_4*beta*np.sum(Zs42)*e/e_0*1000
        Bs41 = (-Z/9)*(3/4)*np.sqrt(5/(2*np.pi))*r_4*beta*np.sum(Zs41)*e/e_0*1000
        B40 = (-Z/9)*(3/16)*np.sqrt(1/np.pi)*r_4*beta*np.sum(Z40)*e/e_0*1000
        Bc41 = (-Z/9)*(3/4)*np.sqrt(5/(2*np.pi))*r_4*beta*np.sum(Zc41)*e/e_0*1000
        Bc42 = (-Z/9)*(3/8)*np.sqrt(5/np.pi)*r_4*beta*np.sum(Zc42)*e/e_0*1000
        Bc43 = (-Z/9)*(3/8)*np.sqrt(70/np.pi)*r_4*beta*np.sum(Zc43)*e/e_0*1000
        Bc44 = (-Z/9)*(3/16)*np.sqrt(35/np.pi)*r_4*beta*np.sum(Zc44)*e/e_0*1000

        B44 = Bs44 + Bc44
        B43 = Bs43 + Bc43
        B42 = Bs42 + Bc42
        B41 = Bs41 + Bc41

        
        print("-------- Point Charge Modeling (meV) --------")
        print("B20: " + str(B20))
        print("B21: " + str(B21))
        print("B22: " + str(B22))
        print("B40: " + str(B40))
        print("B41: " + str(B41))
        print("B42: " + str(B42))
        print("B43: " + str(B43))
        print("B44: " + str(B44))

        
        B = np.array([B20, B21, B22, B40, B41, B42, B43, B44])
        
        return B
    

    
    def SOM(self,ion,L,S,d,Z,rho,Qeff):
        elements = LS(L,S)
        B = elements.PC(ion,L,S,d,Z)
        # rho should be in between 0.1 ~ 0.3 (Co in octahedral should be 0.13 to 0.16) 
        # and Qeff should be 0.6 ~ 1
        SOM_2 = rho*((2/(1-rho))**(3))*(Qeff/Z)*B[0:3]
        SOM_4 = rho*((2/(1-rho))**(5))*(Qeff/Z)*B[3:8]
        SOM = np.array([])
        SOM = np.append(SOM,SOM_2)
        SOM = np.append(SOM,SOM_4)
        
        print("-------- Using Simple Overlap Model (meV) -------")
        print("Qeff: " + str(Qeff))
        print("rho: " + str(rho))
        print("B20: " + str(SOM[0]))
        print("B21: " + str(SOM[1]))
        print("B22: " + str(SOM[2]))
        print("B40: " + str(SOM[3]))
        print("B41: " + str(SOM[4]))
        print("B42: " + str(SOM[5]))
        print("B43: " + str(SOM[6]))
        print("B44: " + str(SOM[7]))
        
        return SOM
    
    def r_l(self,ion):
        
        #From Abragam and Bleaney
        
        #   [r-3 (a.u.), r2 (a.u.), r4 (a.u.), lambda(exp) (cm-1), lambda(calc) (cm-1), p (cm-1)]
        
        #  2D    3d1
        if ion == "Sc2":
            r_l = np.array([0,0,0,86,79,0])
            return r_l
        elif ion == "Ti3":
            r_l = np.array([2.552,1.893,7.071,159,154,0])
        elif ion == "V4":
            r_l = np.array([3.684,1.377,3.593,255,248,0])
        
        #  3F    3d2
        elif ion == "Sc1":
            r_l = np.array([0,0,0,0,35,0])
        elif ion == "Ti2":
            r_l = np.array([2.133,2.447,13.17,61,60,0.16])
        elif ion == "V3":
            r_l = np.array([3.217,1.643,5.447,106,104,0.26])
        elif ion == "Cr4":
            r_l = np.array([4.484,1.227,2.906,163,164,0])
        
        #  4F    3d3
        elif ion == "Ti1":
            r_l = np.array([1.706,3.508,31.62,0,29,0])
        elif ion == "V2":
            r_l = np.array([2.748,2.070,9.605,57,55,0.11])
        elif ion == "Cr3":
            r_l = np.array([3.959,1.447,4.297,91,91,0.17])
        elif ion == "Mn4":
            r_l = np.array([5.361,1.104,2.389,135,134,0])
        
        
        #  5D    3d4
        elif ion == "V1":
            r_l = np.array([2.289,2.819,20.71,0,34,0])
        elif ion == "Cr2":
            r_l = np.array([3.451,1.781,7.211,59,58,0.12])
        elif ion == "Mn3":
            r_l = np.array([4.790,1.286,3.446,87,88,0.18])
        elif ion == "Fe4":
            r_l = np.array([6.332,1,1.986,125,129,0.25])
        
        #  6S    3d5
        elif ion == "Cr1":
            r_l = np.array([2.968,2.319,14.14,0,0,0])
        elif ion == "Mn2":
            r_l = np.array([4.250,1.548,5.513,0,0,0])
        elif ion == "Fe3":
            r_l = np.array([5.724,1.150,2.789,0,0,0])
        elif ion == "Co4":
            r_l = np.array([7.421,0.908,1.659,0,0,0])
        
        #  5D    3d6
        elif ion == "Mn1":
            r_l = np.array([3.683,2.026,10.87,-64,-64,0])
        elif ion == "Fe2":
            r_l = np.array([5.081,1.393,4.496,-114,-103,0.18])
        elif ion == "Co3":
            r_l = np.array([6.699,1.049,2.342,-145,0,0])
        elif ion == "Ni4":
            r_l = np.array([8.552,0.8371,1.423,-197,0,0])
        
        #  4F    3d7
        elif ion == "Fe1":
            r_l = np.array([0,1.774,8.385,-115,-119,0])
        elif ion == "Co2":
            r_l = np.array([6.035,1.251,3.655,-189,-178,0.24])
        elif ion == "Ni3":
            r_l = np.array([7.79,0.9582,1.971,-272,0,0])
        elif ion == "Cu4":
            r_l = np.array([9.814,0.7719,1.221,-320,0,0])
        
        #  3F    3d8
        elif ion == "Co1":
            r_l = np.array([5.388,1.576,6.637,-228,-228,0])
        #elif ion == "Ni2": # r^4 is taken from "NiO - from first principles" by Radwanski and Ropka
        #    r_l = np.array([7.094,1.130,10.5,-343,-324,0.53])
        elif ion == "Ni2":
            r_l = np.array([7.094,1.130,3.003,-343,-324,0.53])
        elif ion == "Cu3":
            r_l = np.array([9.018,0.8763,1.662,-438,0,0])
        
        #  2D    3d9
        elif ion == "Ni1":
            r_l = np.array([0,1.401,5.264,-605,0,0])
        elif ion == "Cu2":
            r_l = np.array([8.252,1.028,2.498,-830,-830,0])
        
        
        return r_l
    
    def theta(self,ion):
        if ion == "Sc2":
            theta = np.array([-2/21,2/63])
            return theta
        elif ion == "Ti3":
            theta = np.array([-2/21,2/63])  
        elif ion == "V4":
            theta = np.array([-2/21,2/63])
        
        #  3F    3d2
        elif ion == "Sc1":
            theta = np.array([-2/105,-2/315])
        elif ion == "Ti2":
            theta = np.array([-2/105,-2/315])
        elif ion == "V3":
            theta = np.array([-2/105,-2/315])
        elif ion == "Cr4":
            theta = np.array([-2/105,-2/315])
        
        #  4F    3d3
        elif ion == "Ti1":
            theta = np.array([2/105,2/315])
        elif ion == "V2":
            theta = np.array([2/105,2/315])
        elif ion == "Cr3":
            theta = np.array([2/105,2/315])
        elif ion == "Mn4":
            theta = np.array([2/105,2/315])
        
        
        #  5D    3d4
        elif ion == "V1":
            theta = np.array([2/21,-2/63])
        elif ion == "Cr2":
            theta = np.array([2/21,-2/63])
        elif ion == "Mn3":
            theta = np.array([2/21,-2/63])
        elif ion == "Fe4":
            theta = np.array([2/21,-2/63])
        
        #  6S    3d5
        elif ion == "Cr1":
            pass
        elif ion == "Mn2":
            pass
        elif ion == "Fe3":
            pass
        elif ion == "Co4":
            pass
        
        #  5D    3d6
        elif ion == "Mn1":
            theta = np.array([-2/21,2/63])
        elif ion == "Fe2":
            theta = np.array([-2/21,2/63])
        elif ion == "Co3":
            theta = np.array([-2/21,2/63])
        elif ion == "Ni4":
            theta = np.array([-2/21,2/63])
        
        #  4F    3d7
        elif ion == "Fe1":
            theta = np.array([-2/107,-2/315])
        elif ion == "Co2":
            theta = np.array([-2/107,-2/315])
        elif ion == "Ni3":
            theta = np.array([-2/107,-2/315])
        elif ion == "Cu4":
            theta = np.array([-2/107,-2/315])
        
        #  3F    3d8
        elif ion == "Co1":
            theta = np.array([2/105,2/315])
        elif ion == "Ni2":
            theta = np.array([2/105,2/315])
        elif ion == "Cu3":
            theta = np.array([2/105,2/315])
        
        #  2D    3d9
        elif ion == "Ni1":
            theta = np.array([2/21,-2/63])
        elif ion == "Cu2":
            theta = np.array([2/21,-2/63])

        return theta
    
    def Diag(self,H,printfunction=False,SaveFig=False):
        E = np.linalg.eigh(H)
        E_val = E[0]
        E_state = E[1].T

        E_excitation = E_val-E_val[0]
        
        if printfunction == False:
            pass
        else:
            print(E_val)
            print(E_excitation)
        
        
        if printfunction == False:
            pass
        else:
            fig = plt.figure()
            
            ax1 = fig.add_subplot(221)
            ax1.eventplot(E_val,orientation='vertical',linelength=0.05,linestyles='solid',colors='r')
            ax1.tick_params(axis='x',which='both',bottom=False,top=False)
            ax1.set_xlabel('Crystal Field Splitting')
            ax1.set_ylabel('Energy (meV)')
            ax1.set_xticklabels('')
            ax1.set_title('Point Charge Model')
            ax1.set_xlim([0.95,1.05])
            
            ax2 = fig.add_subplot(222)
            ax2.eventplot(E_excitation,orientation='horizontal',linelength=0.05,linestyles='solid',colors='b')
            ax2.tick_params(axis='y',which='both',bottom=False,top=False)
            ax2.set_xlabel('Energy (meV)')
            ax2.set_ylabel('Arb. Units')
            ax2.set_yticklabels('')
            ax2.set_yticks([])
            ax2.set_title('Crystal Field excitations')
            ax2.set_ylim([0.95,1.05])
        if SaveFig == False:
            pass
        else:
            plt.savefig("Crystal_Field.pdf")
        plt.show()
        
        
        
        return E_val,E_excitation,E_state
    
    def Zeeman_z(self,L,S,H,Fields,savefunction=False):
        print('--------- Applying Magnetic Field inducing Zeeman splitting in Z-axis --------')
        elements = LS(L,S)
        Fields = np.array(range(0,Fields+2,1))
        #mu_B = 9.27400999e-24 #J/T
        mu_B = (5.7883818012e-5)*1000 #meV
        Lz = elements.Lz(L,S)
        Sz = elements.Sz(L,S)
        g_e = 2.002324
        H_m = mu_B*(Lz + g_e*Sz)
        E_val = np.zeros((len(Fields),elements.LS_degeneracy))
        E_excitation = np.zeros((len(Fields),elements.LS_degeneracy))
        for k in range(len(Fields)):
            H_M = H + H_m*Fields[k]
            E_val[k],E_excitation[k],E_state = elements.Diag(H_M,printfunction=False)
        
        for k in range(Fields[k]):
            for k1 in range(elements.LS_degeneracy):
                plt.plot(E_excitation[k][k1],Fields[k],'|',markersize=8)
                
        plt.title(r'Zeeman Splitting (H $\Vert$ c)')
        plt.xlabel('Energy (meV)')
        plt.ylabel('Magnetic Field (T)')
        if savefunction == False:
            pass
        else:
            plt.savefig("Zeeman_z.pdf")
        plt.show()    
        
        return E_val
    
    def Zeeman_x(self,L,S,H,Fields,savefunction=False):
        print('--------- Applying Magnetic Field inducing Zeeman splitting in X-axis --------')
        elements = LS(L,S)
        Fields = np.array(range(0,Fields+2,1))
        #mu_B = 9.27400999e-24 #J/T
        mu_B = (5.7883818012e-5)*1000 #meV
        Lx = elements.Lx(L,S)
        Sx = elements.Sx(L,S)
        g_e = 2.002324
        H_m = mu_B*(Lx + g_e*Sx)
        E_val = np.zeros((len(Fields),elements.LS_degeneracy))
        E_excitation = np.zeros((len(Fields),elements.LS_degeneracy))
        for k in range(len(Fields)):
            H_M = H + H_m*Fields[k]
            E_val[k],E_excitation[k],E_state = elements.Diag(H_M,printfunction=False)
            
        for k in range(Fields[k]):
            for k1 in range(elements.LS_degeneracy):
                plt.plot(E_excitation[k][k1],Fields[k],'|',markersize=8)
                
        plt.title(r'Zeeman Splitting (H $\bot$ c)')
        plt.xlabel('Energy (meV)')
        plt.ylabel('Magnetic Field (T)')
        if savefunction == False:
            pass
        else:
            plt.savefig("Zeeman_x.pdf")
        plt.show() 
        
        return E_val
    
    def Zeeman_y(self,L,S,H,Fields,savefunction=False):
        print('--------- Applying Magnetic Field inducing Zeeman splitting in Y-axis --------')
        elements = LS(L,S)
        Fields = np.array(range(0,Fields+2,1))
        #mu_B = 9.27400999e-24 #J/T
        mu_B = (5.7883818012e-5)*1000 #meV
        Ly = elements.Ly(L,S)
        Sy = elements.Sy(L,S)
        g_e = 2.002324
        H_m = mu_B*(Ly + g_e*Sy)
        E_val = np.zeros((len(Fields),elements.LS_degeneracy))
        E_excitation = np.zeros((len(Fields),elements.LS_degeneracy))
        for k in range(len(Fields)):
            H_M = H + H_m*Fields[k]
            E_val[k],E_excitation[k],E_state = elements.Diag(H_M,printfunction=False)
        
        for k in range(Fields[k]):
            for k1 in range(elements.LS_degeneracy):
                plt.plot(E_excitation[k][k1],Fields[k],'|',markersize=8)
                
        plt.title(r'Zeeman Splitting (H $\bot$ c)')
        plt.xlabel('Energy (meV)')
        plt.ylabel('Magnetic Field (T)')
        if savefunction == False:
            pass
        else:
            plt.savefig("Zeeman_y.pdf")
        plt.show() 
        
        return E_val
    
    
    def MagMoment_z(self,L,S,H):
        print('\n')
        print('--------- Magnetic Properties (z-axis) --------')
        elements = LS(L,S)
        Lz = elements.Lz(L,S)
        mu_B = (5.7883818012e-5)*1000 #meV
        g_e = 2.0023193
        Sz = elements.Sz(L,S)
        B = 0.1
        H = H + mu_B*(Lz + g_e*Sz)*B
        H_e,H_v = np.linalg.eigh(H)
        E_final = H_e - H_e[0]
        H_vt = H_v.T
        print(E_final)
        print('----------------- Lz ---------------')
        Lz_m = np.zeros((elements.LS_degeneracy))
        for k in range(elements.LS_degeneracy):
            Lz_m[k] = np.dot(np.conjugate(H_vt[k]),np.dot(Lz,H_vt[k]))
        print(Lz_m)
        print('----------------- Sz ---------------')
        Sz = elements.Sz(L,S)
        Sz_m = np.zeros((elements.LS_degeneracy))
        for k in range(elements.LS_degeneracy):
                Sz_m[k] = np.dot(np.conjugate(H_vt[k]),np.dot(Sz,H_vt[k]))
        print(Sz_m)
        print('----------------- Mz ---------------')
        Mz = Lz_m + 2*Sz_m
        print(Mz)
        
        return Lz_m,Sz_m,Mz
        
    def MagMoment_x(self,L,S,H):
        elements = LS(L,S)
        print('\n')
        print('--------- Magnetic Properties (x-axis) --------')
        print("----------- Applying magnetic field (0.01 T) to break degeneracy -----------")
        mu_B = 5.7883818012e-2 # meV
        Lx = elements.Lx(L,S)
        g_e = 2.0023193
        Sx = elements.Sx(L,S)
        B = 0.1
        H = H + mu_B*(Lx + g_e*Sx)*B
        H_e,H_v = np.linalg.eigh(H)
        E_final = H_e - H_e[0]
        H_vt = H_v.T
        print(E_final)
        print('--------------- Lx ---------------')
        Lx_m = np.zeros((elements.LS_degeneracy))
        for k in range(elements.LS_degeneracy):
            Lx_m[k] = np.dot(np.conjugate(H_vt[k]),np.dot(Lx,H_vt[k]))
        print(Lx_m)
        print('--------------- Sx ---------------')
        Sx_m = np.zeros((elements.LS_degeneracy))
        for k in range(elements.LS_degeneracy):
            Sx_m[k] = np.dot(np.conjugate(H_vt[k]),np.dot(Sx,H_vt[k]))
        print(Sx_m)
        print('--------------- Mx ---------------')
        Mx = Lx_m + g_e*Sx_m
        print(Mx)
        
        return Lx_m,Sx_m,Mx
    
    def MagMoment_y(self,L,S,H):
        elements = LS(L,S)
        print('\n')
        print('--------- Magnetic Properties (y-axis) --------')
        print("----------- Applying magnetic field (0.01 T) to break degeneracy -----------")
        mu_B = 5.7883818012e-2 # meV
        Ly = elements.Ly(L,S)
        g_e = 2.0023193
        Sy = elements.Sy(L,S)
        B = 0.1
        H = H + mu_B*(Ly + g_e*Sy)*B
        H_e,H_v = np.linalg.eigh(H)
        E_final = H_e - H_e[0]
        H_vt = H_v.T
        print(E_final)
        print('--------------- Ly ---------------')
        Ly_m = np.zeros((elements.LS_degeneracy))
        for k in range(elements.LS_degeneracy):
            Ly_m[k] = np.dot(np.conjugate(H_vt[k]),np.dot(Ly,H_vt[k]))
        print(Ly_m)
        print('--------------- Sy ---------------')
        Sy_m = np.zeros((elements.LS_degeneracy))
        for k in range(elements.LS_degeneracy):
            Sy_m[k] = np.dot(np.conjugate(H_vt[k]),np.dot(Sy,H_vt[k]))
        print(Sy_m)
        print('--------------- My ---------------')
        My = Ly_m + g_e*Sy_m
        print(My)
        
        return Ly_m,Sy_m,My
    
    def SpecificHeat(self,E):
        print('\n')
        print('--------- Specific Heat (J/mol/K) --------')
        ################## specific heat ###############

        R = 8.314#*1.6021773e-22    # J/molK
        #kB_J = 1.380649e-23#*6.241506363094e+21  # J/K
        kB = 8.617333e-2 # meV/K
        #Na = R/kB_J
        T = np.array(range(1,301,1))  # K

        #a = dict(collections.Counter(E))
        #print(a)
        #degeneracies = np.array(list(a.items()))
        #E = E#*1.60217733e-22 # meV to J
        #E1 = degeneracies[:,0]*1.60217733e-22 # meV to J
        #g = degeneracies[:,1]
        C_schottky = np.zeros((len(T)))
        for k in range(len(T)):
            C_schottky[k] = R*np.sum(E*np.exp(-E/(kB*T[k]))*(E+E*np.sum(np.exp(-E/(kB*T[k])))-np.sum(E*np.exp(-E/(kB*T[k])))))/(T[k]**2*(1+np.sum(np.exp(-E/(kB*T[k]))))**2)
        
        plt.plot(T,C_schottky,'ro')
        plt.xlabel('Temperature (K)')
        plt.ylabel('specific heat (J/mol/K)')
        plt.title('Specific Heat')
        #plt.legend('')
        plt.grid(alpha=.4,linestyle='--')
        plt.xlim(0,300)
        plt.show()
        
        return C_schottky,T
    
    def Entropy(self,E):
        print('\n')
        print('--------- Magnetic entropy (J/mol/K) --------')
        ################## magnetic entropy ###############
        R = 8.314#*6.241506363094e+21     # J/molK
        #kB = 1.380649e-23#*6.241506363094e+21  # J/K
        kB = 8.617333e-2 # meV/K
        T = np.array(range(1,301,1))
        #a = dict(collections.Counter(E))
        #print(a)
        #degeneracies = np.array(list(a.items()))
        #E = E*1.60217733e-22 # meV to J
        #E1 = degeneracies[:,0]*1.60217733e-22 # meV to J
        #g = degeneracies[:,1]
        entropy = np.zeros((len(T)))
        for k in range(len(T)):
            entropy[k] = R*np.log(np.sum(np.exp(-E/(kB*T[k]))))
        
        plt.plot(T,entropy,'ro')
        plt.xlabel('Temperature (K)')
        plt.ylabel('Magnetic entropy (J/mol/K)')
        plt.title('Magnetic entropy')
        #plt.legend('')
        plt.grid(alpha=.4,linestyle='--')
        plt.xlim(0,300)
        plt.show()
        
        return entropy,T
    
    def Susceptibility_z(self,L,S,H,SO):
        elements = LS(L,S)
        T = np.array(range(1,301,1))
        R = 8.314#*6.241506363094e+21     # J/molK
        kB = 1.380649e-23#*6.241506363094e+21  # J/K
        Na = kB/R
        mu_B = 5.7883818012e-2 # meV
        g_e = 2.002324
        g_l = 1
        Lz = elements.Lz(L,S)
        Sz = elements.Sz(L,S)
        Jz = Lz + Sz
        B = 0.1
        H = H + mu_B*(Lz + g_e*Sz)*B
        H_e,H_v = np.linalg.eigh(H)
        E_val = H_e - H_e[0]
        E_state = H_v.T
        E_val = E_val/6.241506363094e+21 #meV to J
        X_z = np.zeros((len(T)))
        Z = np.zeros((len(T)))
        
        for k in range(len(T)):
            X1 = 0
            X2 = 0
            X3 = 0
            Z[k] = np.sum(np.exp(-E_val/(kB*T[k])))
            for k1 in range(len(E_val)):
                X3 += (-1/(kB*T[k]))*(np.sum(np.dot(np.conjugate(E_state[k1]),np.dot(Jz,E_state[k1])))**2*np.exp(-E_val[k1]/(kB*T[k])))**2
                for k2 in range(len(E_val)):
                    if k1 == k2:
                        X1 += np.sum(np.dot(np.conjugate(E_state[k2]),np.dot(Jz,E_state[k1])))**2/(kB*T[k])*np.exp(-E_val[k1]/(kB*T[k]))
                        X2 += 0
                        
                    else:
                        X1 += 0
                        X2 += 2*(np.sum(np.dot(np.conjugate(E_state[k2]),np.dot(Jz,E_state[k1])))**2/(E_val[k2]-E_val[k1]))*np.exp(-E_val[k1]/(kB*T[k]))
                        
            X_z[k] = (Na*SO**2*mu_B**2/Z[k])*(X1 + X2 - X3)
        
        
        
        plt.plot(T,X_z,'ro')
        plt.xlabel('Temperature (K)')
        plt.ylabel('Magnetic Susceptibility ($\mu_B$) (z-axis)')
        plt.title('Magnetic Susceptibility (z-axis)')
        #plt.legend('')
        plt.grid(alpha=.4,linestyle='--')
        plt.xlim(0,300)
        plt.show()
                    
        return X_z,T
    
    def Susceptibility_x(self,L,S,H,SO):
        elements = LS(L,S)
        T = np.array(range(1,301,1))
        R = 8.314#*6.241506363094e+21     # J/molK
        kB = 1.380649e-23#*6.241506363094e+21  # J/K
        Na = kB/R
        mu_B = 5.7883818012e-2 # meV
        g_e = 2.002324
        g_l = 1
        Lx = elements.Lx(L,S)
        Sx = elements.Sx(L,S)
        Jx = Lx + Sx
        B = 0.01
        H = H + mu_B*(Lx + g_e*Sx)*B
        H_e,H_v = np.linalg.eigh(H)
        E_val = H_e - H_e[0]
        E_state = H_v.T
        E_val = E_val/6.241506363094e+21 #meV to J
        X_x = np.zeros((len(T)))
        Z = np.zeros((len(T)))
        
        for k in range(len(T)):
            X1 = 0
            X2 = 0
            X3 = 0
            Z[k] = np.sum(np.exp(-E_val/(kB*T[k])))
            for k1 in range(len(E_val)):
                X3 += (-1/(kB*T[k]))*(np.sum(np.dot(np.conjugate(E_state[k1]),np.dot(Jx,E_state[k1])))**2*np.exp(-E_val[k1]/(kB*T[k])))**2
                for k2 in range(len(E_val)):
                    if k1 == k2:
                        X1 += np.sum(np.dot(np.conjugate(E_state[k2]),np.dot(Jx,E_state[k1])))**2/(kB*T[k])*np.exp(-E_val[k1]/(kB*T[k]))
                        X2 += 0
                        
                    else:
                        X1 += 0
                        X2 += 2*(np.sum(np.dot(np.conjugate(E_state[k2]),np.dot(Jx,E_state[k1])))**2/(E_val[k2]-E_val[k1]))*np.exp(-E_val[k1]/(kB*T[k]))
                        
            X_x[k] = (Na*SO**2*mu_B**2/Z[k])*(X1 + X2 - X3)
        
        
        plt.plot(T,X_x,'bo')
        plt.xlabel('Temperature (K)')
        plt.ylabel('Magnetic Susceptibility ($\mu_B$) (x-axis)')
        plt.title('Magnetic Susceptibility (x-axis)')
        #plt.legend('')
        plt.grid(alpha=.4,linestyle='--')
        plt.xlim(0,300)
        plt.show()
                    
        return X_x,T
    
    def Susceptibility_y(self,L,S,H,SO):
        elements = LS(L,S)
        T = np.array(range(1,301,1))
        R = 8.314#*6.241506363094e+21     # J/molK
        kB = 1.380649e-23#*6.241506363094e+21  # J/K
        Na = kB/R
        mu_B = 5.7883818012e-2 # meV
        g_e = 2.002324
        g_l = 1
        Ly = elements.Ly(L,S)
        Sy = elements.Sy(L,S)
        Jy = Ly + Sy
        B = 0.01
        H = H + mu_B*(Ly + g_e*Sy)*B
        H_e,H_v = np.linalg.eigh(H)
        E_val = H_e - H_e[0]
        E_state = H_v.T
        E_val = E_val/6.241506363094e+21 #meV to J
        X_y = np.zeros((len(T)))
        Z = np.zeros((len(T)))
        
        for k in range(len(T)):
            X1 = 0
            X2 = 0
            X3 = 0
            Z[k] = np.sum(np.exp(-E_val/(kB*T[k])))
            for k1 in range(len(E_val)):
                X3 += (-1/(kB*T[k]))*(np.sum(np.dot(np.conjugate(E_state[k1]),np.dot(Jy,E_state[k1])))**2*np.exp(-E_val[k1]/(kB*T[k])))**2
                for k2 in range(len(E_val)):
                    if k1 == k2:
                        X1 += np.sum(np.dot(np.conjugate(E_state[k2]),np.dot(Jy,E_state[k1])))**2/(kB*T[k])*np.exp(-E_val[k1]/(kB*T[k]))
                        X2 += 0
                        
                    else:
                        X1 += 0
                        X2 += 2*(np.sum(np.dot(np.conjugate(E_state[k2]),np.dot(Jy,E_state[k1])))**2/(E_val[k2]-E_val[k1]))*np.exp(-E_val[k1]/(kB*T[k]))
                        
            X_y[k] = (Na*SO**2*mu_B**2/Z[k])*(X1 + X2 - X3)
                    
        plt.plot(T,X_y,'yo')
        plt.xlabel('Temperature (K)')
        plt.ylabel('Magnetic Susceptibility ($\mu_B$) (y-axis)')
        plt.title('Magnetic Susceptibility (y-axis)')
        #plt.legend('')
        plt.grid(alpha=.4,linestyle='--')
        plt.xlim(0,300)
        plt.show()
                    
        return X_y,T


    """    
    E_exp = np.array([30.0, 50.0, 107.0, 136.0, 170.0, 737.0])
    X = np.array([Qeff,rho])



    def Fit_B(B,X,E_exp):
        #Qeff = X[0]
        #rho = X[1]
        SO = B[8]
        #SOM_2 = rho*((2/(1-rho))**(3))*(Qeff/Z)*B[0:3]
        #SOM_4 = rho*((2/(1-rho))**(5))*(Qeff/Z)*B[3:8]
        #SOM = np.array([])
        #SOM = np.append(SOM,SOM_2)
        #SOM = np.append(SOM,SOM_4)
        SO_matrix = Co2_oct.SO(ion,L,S,SO)
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
    



    X = np.array([Qeff,rho])
    def Fit_SOM(X,B,E_exp):
        Qeff,rho = X
        SOM_2 = rho*((2/(1-rho))**(3))*(Qeff/Z)*B[0:3]
        SOM_4 = rho*((2/(1-rho))**(5))*(Qeff/Z)*B[3:8]
        SOM = np.array([])
        SOM = np.append(SOM,SOM_2)
        SOM = np.append(SOM,SOM_4)
        SO_matrix = Co2_oct.SO(ion,L,S,SO)
        H = SOM[0]*O20 + SOM[1]*O21 + SOM[2]*O22 + SOM[3]*O40 + SOM[4]*O41 + SOM[5]*O42 + SOM[6]*O43 + SOM[7]*O44 + SO_matrix
        E = np.linalg.eigh(H)
        E_val = E[0]
        #E_state = E[1].T
        E_excitation = E_val-E_val[0]
        E_calc = np.array([E_excitation[4], E_excitation[8], E_excitation[12], E_excitation[16], E_excitation[20], E_excitation[24]])
        
    return np.sum((E_calc/E_exp - 1)**2)

    ### rho is in between 0.1 - 0.3 while Qeff is in between 0 - 1
    
    bnds1 = (0,1)
    bnds2 = (0,1)
    bounds = [bnds1,bnds2]

    result = scipy.optimize.minimize(Fit_SOM, X, method = 'Powell', args = (B,E_exp),bounds=bounds,tol=1615)
    X = result.x
    print(result)
    print("New X(Qeff,rho): " + str(X))
    """
        
                
'''
    
    

    def FormFactor(self,ion):
                #####https://www.ill.eu/sites/ccsl/ffacts/ffactnode6.html############
        ####  r0
        #Ion	A	a	B	b	C	c	D
        Sc0 = [0.2512,	90.0296,	0.3290,	39.4021, 0.4235, 14.3222,	-0.0043]
        Sc1 = [0.4889,	51.1603,	0.5203,	14.0764, -0.0286, 0.1792,	0.0185]
        Sc2	= [0.5048,	31.4035,	0.5186,	10.9897, -0.0241, 1.1831,	0.0000]
        Ti0	= [0.4657,	33.5898,	0.5490,	9.8791,	-0.0291, 0.3232,	0.0123]
        Ti1	= [0.5093,	36.7033,	0.5032,	10.3713, -0.0263, 0.3106,	0.0116]
        Ti2	= [0.5091,	24.9763,	0.5162,	8.7569,	-0.0281, 0.9160,	0.0015]
        Ti3	= [0.3571,	22.8413,	0.6688,	8.9306,	-0.0354, 0.4833,	0.0099]
        V0	= [0.4086,	28.8109,	0.6077,	8.5437,	-0.0295, 0.2768,	0.0123]
        V1	= [0.4444,	32.6479,	0.5683,	9.0971,	-0.2285, 0.0218,	0.2150]
        V2	= [0.4085,	23.8526,	0.6091,	8.2456,	-0.1676, 0.0415,	0.1496]
        V3	= [0.3598,	19.3364,	0.6632,	7.6172,	-0.3064, 0.0296,	0.2835]
        V4	= [0.3106,	16.8160,	0.7198,	7.0487,	-0.0521, 0.3020,	0.0221]
        Cr0	= [0.1135,	45.1990,	0.3481,	19.4931, 0.5477, 7.3542,	-0.0092]
        Cr1	= [-0.0977,	0.0470,	0.4544,	26.0054, 0.5579, 7.4892,	0.0831]
        Cr2	= [1.2024,	-0.0055,	0.4158,	20.5475, 0.6032, 6.9560,	-1.2218]
        Cr3	= [-0.3094,	0.0274,	0.3680,	17.0355, 0.6559, 6.5236,	0.2856]
        Cr4	= [-0.2320,	0.0433,	0.3101,	14.9518, 0.7182, 6.1726,	0.2042]
        Mn0	= [0.2438,	24.9629,	0.1472,	15.6728, 0.6189, 6.5403,	-0.0105]
        Mn1	= [-0.0138,	0.4213,	0.4231,	24.6680, 0.5905, 6.6545,	-0.0010]
        Mn2	= [0.4220,	17.6840,	0.5948,	6.0050,	0.0043,	-0.6090,	-0.0219]
        Mn3	= [0.4198,	14.2829,	0.6054,	5.4689,	0.9241,	-0.0088,	-0.9498]
        Mn4	= [0.3760,	12.5661,	0.6602,	5.1329,	-0.0372, 0.5630,	0.0011]
        Fe0	= [0.0706,	35.0085,	0.3589,	15.3583, 0.5819, 5.5606,	-0.0114]
        Fe1	= [0.1251,	34.9633,	0.3629,	15.5144, 0.5223, 5.5914,	-0.0105]
        Fe2	= [0.0263,	34.9597,	0.3668,	15.9435, 0.6188, 5.5935,	-0.0119]
        Fe3	= [0.3972,	13.2442,	0.6295,	4.9034,	-0.0314, 0.3496,	0.0044]
        Fe4	= [0.3782,	11.3800,	0.6556,	4.5920,	-0.0346, 0.4833,	0.0005]
        Co0	= [0.4139,	16.1616,	0.6013,	4.7805,	-0.1518, 0.0210,	0.1345]
        Co1	= [0.0990,	33.1252,	0.3645,	15.1768, 0.5470, 5.0081,	-0.0109]
        Co2	= [0.4332,	14.3553,	0.5857,	4.6077,	-0.0382, 0.1338,	0.0179]
        Co3	= [0.3902,	12.5078,	0.6324,	4.4574,	-0.1500, 0.0343,	0.1272]
        Co4	= [0.3515,	10.7785,	0.6778,	4.2343,	-0.0389, 0.2409,	0.0098]
        Ni0	= [-0.0172,	35.7392,	0.3174,	14.2689, 0.7136, 4.5661,	-0.0143]
        Ni1	= [0.0705,	35.8561,	0.3984,	13.8042, 0.5427, 4.3965,	-0.0118]
        Ni2	= [0.0163,	35.8826,	0.3916,	13.2233, 0.6052, 4.3388,	-0.0133]
        Ni3	= [-0.0134,	35.8677,	0.2678,	12.3326, 0.7614, 4.2369,	-0.0162]
        Ni4	= [-0.0090,	35.8614,	0.2776,	11.7904, 0.7474, 4.2011,	-0.0163]
        Cu0	= [0.0909,	34.9838,	0.4088,	11.4432, 0.5128, 3.8248,	-0.0124]
        Cu1	= [0.0749,	34.9656,	0.4147,	11.7642, 0.5238, 3.8497,	-0.0127]
        Cu2	= [0.0232,	34.9686,	0.4023,	11.5640, 0.5882, 3.8428,	-0.0137]
        Cu3	= [0.0031,	34.9074,	0.3582,	10.9138, 0.6531, 3.8279,	-0.0147]
        Cu4	= [-0.0132,	30.6817,	0.2801,	11.1626, 0.7490, 3.8172,	-0.0165]
        
        #### r2
        #Ion	A	a	B	b	C	c	D
        Sc0=	[10.8172,   54.3270,    4.7353,	14.8471,0.6071,	4.2180,	0.0011]
        Sc1=	[8.5021,	34.2851,	3.2116,	10.9940,0.4244,	3.6055,	0.0009]
        Sc2=	[4.3683,	28.6544,	3.7231,	10.8233,0.6074,	3.6678,	0.0014]
        Ti0=	[4.3583,	36.0556,	3.8230,	11.1328,0.6855,	3.4692,	0.0020]
        Ti1=	[6.1567,	27.2754,	2.6833,	8.9827,	0.4070,	3.0524,	0.0011]
        Ti2=	[4.3107,	18.3484,	2.0960,	6.7970,	0.2984,	2.5476,	0.0007]
        Ti3=	[3.3717,	14.4441,	1.8258,	5.7126,	0.2470,	2.2654,	0.0005]
        V0=	    [3.8099,	21.3471,	2.3295,	7.4089,	0.4333,	2.6324,	0.0015]
        V1=	    [4.7474,	23.3226,	2.3609,	7.8082,	0.4105,	2.7063,	0.0014]
        V2= 	[3.4386,	16.5303,	1.9638,	6.1415,	0.2997,	2.2669,	0.0009]
        V3=	    [2.3005,	14.6821,	2.0364,	6.1304,	0.4099,	2.3815,	0.0014]
        V4=	    [1.8377,	12.2668,	1.8247,	5.4578,	0.3979,	2.2483,	0.0012]
        V0=	    [3.8989,	20.4087,	2.2151,	7.0842,	0.3921,	2.5401,	0.0014]
        Cr0=	[3.4085,	20.1267,	2.1006,	6.8020,	0.4266,	2.3941,	0.0019]
        Cr1=	[3.7768,	20.3456,	2.1028,	6.8926,	0.4010,	2.4114,	0.0017]
        Cr2=	[2.6422,	16.0598,	1.9198,	6.2531,	0.4446,	2.3715,	0.0020]
        Cr3=	[1.6262,	15.0656,	2.0618,	6.2842,	0.5281,	2.3680,	0.0023]
        Cr4=	[1.0293,	13.9498,	1.9933,	6.0593,	0.5974,	2.3457,	0.0027]
        Mn0=	[2.6681,	16.0601,	1.7561,	5.6396,	0.3675,	2.0488,	0.0017]
        Mn1=	[3.2953,	18.6950,	1.8792,	6.2403,	0.3927,	2.2006,	0.0022]
        Mn2=	[2.0515,	15.5561,	1.8841,	6.0625,	0.4787,	2.2323,	0.0027]
        Mn3=	[1.2427,	14.9966,	1.9567,	6.1181,	0.5732,	2.2577,	0.0031]
        Mn4=	[0.7879,	13.8857,	1.8717,	5.7433,	0.5981,	2.1818,	0.0034]
        Fe0=	[1.9405,	18.4733,	1.9566,	6.3234,	0.5166,	2.1607,	0.0036]
        Fe1=	[2.6290,	18.6598,	1.8704,	6.3313,	0.4690,	2.1628,	0.0031]
        Fe2=	[1.6490,	16.5593,	1.9064,	6.1325,	0.5206,	2.1370,	0.0035]
        Fe3=	[1.3602,	11.9976,	1.5188,	5.0025,	0.4705,	1.9914,	0.0038]
        Fe4=	[1.5582,	8.2750,	    1.1863,	3.2794,	0.1366,	1.1068,	-0.0022]
        Co0=	[1.9678,	14.1699,	1.4911,	4.9475,	0.3844,	1.7973,	0.0027]
        Co1=	[2.4097,	16.1608,	1.5780,	5.4604,	0.4095,	1.9141,	0.0031]
        Co2=	[1.9049,	11.6444,	1.3159,	4.3574,	0.3146,	1.6453,	0.0017]
        Co3=	[1.7058,	8.8595,	    1.1409,	3.3086,	0.1474,	1.0899,	-0.0025]
        Co4=	[1.3110,	8.0252,	    1.1551,	3.1792,	0.1608,	1.1301,	-0.0011]
        Ni0=	[1.0302,	12.2521,	1.4669,	4.7453,	0.4521,	1.7437,	0.0036]
        Ni1=	[2.1040,	14.8655,	1.4302,	5.0714,	0.4031,	1.7784,	0.0034]
        Ni2=	[1.7080,	11.0160,	1.2147,	4.1031,	0.3150,	1.5334,	0.0018]
        Ni3=	[1.1612,	7.7000,	    1.0027,	3.2628,	0.2719,	1.3780,	0.0025]
        Ni4=	[1.1612,	7.7000,	    1.0027,	3.2628,	0.2719,	1.3780,	0.0025]
        Cu0=	[1.9182,	14.4904,	1.3329,	4.7301,	0.3842,	1.6394,	0.0035]
        Cu1=	[1.8814,	13.4333,	1.2809,	4.5446,	0.3646,	1.6022,	0.0033]
        Cu2=	[1.5189,	10.4779,	1.1512,	3.8132,	0.2918,	1.3979,	0.0017]
        Cu3=	[1.2797,	8.4502,	    1.0315,	3.2796,	0.2401,	1.2498,	0.0015]
        Cu4=	[0.9568,	7.4481,	    0.9099,	3.3964,	0.3729,	1.4936,	0.0049]

        
        ######  r4
        #Ion	A	a	B	b	C	c	D
        Sc0=	[1.3420,	10.2000,	0.3837,	3.0786,	0.0468,	0.1178,	-0.0328]
        Sc1=	[7.1167,	15.4872,	-6.6671,18.2692,0.4900,	2.9917,	0.0047]
        Sc2=	[-1.6684,	15.6475,	1.7742,	9.0624,	0.4075,	2.4116,	0.0042]
        Ti0=	[-2.1515,	11.2705,	2.5149,	8.8590,	0.3555,	2.1491,	0.0045]
        Ti1=	[-1.0383,	16.1899,	1.4699,	8.9239,	0.3631,	2.2834,	0.0044]
        Ti2=	[-1.3242,	15.3096,	1.2042,	7.8994,	0.3976,	2.1562,	0.0051]
        Ti3=	[-1.1117,	14.6349,	0.7689,	6.9267,	0.4385,	2.0886,	0.0060]
        V0=	    [-0.9633,	15.2729,	0.9274,	7.7315,	0.3891,	2.0530,	0.0063]
        V1=	    [-0.9606,	15.5451,	1.1278,	8.1182,	0.3653,	2.0973,	0.0056]
        V2=	    [-1.1729,	14.9732,	0.9092,	7.6131,	0.4105,	2.0391,	0.0067]
        V3=	    [-0.9417,	14.2045,	0.5284,	6.6071,	0.4411,	1.9672,	0.0076]
        V4=	    [-0.7654,	13.0970,	0.3071,	5.6739,	0.4476,	1.8707,	0.0081]
        Cr0=	[-0.6670,	19.6128,	0.5342,	6.4779,	0.3641,	1.9045,	0.0073]
        Cr1=	[-0.8309,	18.0428,	0.7252,	7.5313,	0.3828,	2.0032,	0.0073]
        Cr2=	[-0.8930,	15.6641,	0.5590,	7.0333,	0.4093,	1.9237,	0.0081]
        Cr3=	[-0.7327,	14.0727,	0.3268,	5.6741,	0.4114,	1.8101,	0.0085]
        Cr4=	[-0.6748,	12.9462,	0.1805,	6.7527,	0.4526,	1.7999,	0.0098]
        Mn0=	[-0.5452,	15.4713,	0.4406,	4.9024,	0.2884,	1.5430,	0.0059]
        Mn1=	[-0.7947,	17.8673,	0.6078,	7.7044,	0.3798,	1.9045,	0.0087]
        Mn2=	[-0.7416,	15.2555,	0.3831,	6.4693,	0.3935,	1.7997,	0.0093]
        Mn3=	[-0.6603,	13.6066,	0.2322,	6.2175,	0.4104,	1.7404,	0.0101]
        Mn4=	[-0.5127,	13.4613,	0.0313,	7.7631,	0.4282,	1.7006,	0.0113]
        Fe0=	[-0.5029,	19.6768,	0.2999,	3.7762,	0.2576,	1.4241,	0.0071]
        Fe1=	[-0.5109,	19.2501,	0.3896,	4.8913,	0.2810,	1.5265,	0.0069]
        Fe2=	[-0.5401,	17.2268,	0.2865,	3.7422,	0.2658,	1.4238,	0.0076]
        Fe3=	[-0.5507,	11.4929,	0.2153,	4.9063,	0.3468,	1.5230,	0.0095]
        Fe4=	[-0.5352,	9.5068,  	0.1783,	5.1750,	0.3584,	1.4689,	0.0097]
        Co0=	[-0.4221,	14.1952,	0.2900,	3.9786,	0.2469,	1.2859,	0.0063]
        Co1=	[-0.4115,	14.5615,	0.3580,	4.7170,	0.2644,	1.4183,	0.0074]
        Co2=	[-0.4759,	14.0462,	0.2747,	3.7306,	0.2458,	1.2504,	0.0057]
        Co3=	[-0.4466,	13.3912,	0.1419,	3.0110,	0.2773,	1.3351,	0.0093]
        Co4=	[-0.4091,	13.1937,	-0.0194, 3.4169,0.3534,	1.4214,	0.0112]
        Ni0=	[-0.4428,	14.4850,	0.0870,	3.2345,	0.2932,	1.3305,	0.0096]
        Ni1=	[-0.3836,	13.4246,	0.3116,	4.4619,	0.2471,	1.3088,	0.0079]
        Ni2=	[-0.3803,	10.4033,	0.2838,	3.3780,	0.2108,	1.1036,	0.0050]
        Ni3=	[-0.3715,	8.9516,  	0.1211,	2.9399,	0.2526,	1.1049,	0.0061]
        Ni4=	[-0.3509,	8.1572,  	0.2220,	2.1063,	0.1567,	0.9253,	0.0065]
        Cu0=	[-0.3204,	15.1324,	0.2335,	4.0205,	0.2312,	1.1957,	0.0068]
        Cu1=	[-0.3572,	15.1251,	0.2336,	3.9662,	0.2315,	1.1967,	0.0070]
        Cu2=	[-0.3914,	14.7400,	0.1275,	3.3840,	0.2548,	1.2552,	0.0103]
        Cu3=	[-0.3671,	14.0816,	-0.0078, 3.3149,0.3154,	1.3767,	0.0132]
        Cu4=	[-0.2915,	14.1243,	-0.1065, 4.2008,0.3247,	1.3516,	0.0148]





'''




class J():
    def __init__(self,Jvalue):
        self.J = Jvalue
        if J == 0 or 1 or 2 or 3 or 4 or 5 or 6 or 7 or 8 or 9:
            self.J_array = np.arange(-Jvalue,Jvalue+1)
        else:
            J1 = int(self.J + 0.5)
            self.J_array = np.arange(-J1,J1) + 0.5
        self.J_degeneracy = int(2*Jvalue+1)
        self.J_matrix = np.identity(self.J_degeneracy)*self.J_array
        self.J_tot = np.identity(self.J_degeneracy)*Jvalue
        self.i = np.identity(self.J_degeneracy)
        self.J_matrix_1 = np.zeros((self.J_degeneracy,self.J_degeneracy))
        self.J_tot = np.identity(self.J_degeneracy)*self.J

    def Jz(self,Jvalue):
        elements = J(Jvalue)
        elements.J_matrix = elements.J_matrix
        return elements.J_matrix
    
    def Jplus(self,Jvalue):
        elements = J(Jvalue)
        Jplus = np.sqrt(elements.J*(elements.J+1)-elements.J_array*(elements.J_array+1))
        return Jplus
        
    def Jp(self,Jvalue):
        elements = J(Jvalue)
        Jplus = elements.Jplus(Jvalue)
        for i in range(elements.J_degeneracy-1):
            elements.J_matrix_1[i][i+1] = Jplus[i]
        return elements.J_matrix_1
    
    def Jp2(self,Jvalue):
        elements = J(Jvalue)
        Jp2 = np.dot(elements.Jp(Jvalue),elements.Jp(Jvalue))
        return Jp2
    
    def Jp3(self,Jvalue):
        elements = J(Jvalue)
        Jp3 = np.dot(elements.Jp(Jvalue),elements.Jp2(Jvalue))
        return Jp3
    
    def Jp4(self,Jvalue):
        elements = J(Jvalue)
        Jp4 = np.dot(elements.Jp(Jvalue),elements.Jp3(Jvalue))
        return Jp4
    
    def Jp5(self,Jvalue):
        elements = J(Jvalue)
        Jp5 = np.dot(elements.Jp(Jvalue),elements.Jp4(Jvalue))
        return Jp5
    
    def Jp6(self,Jvalue):
        elements = J(Jvalue)
        Jp6 = np.dot(elements.Jp(Jvalue),elements.Jp5(Jvalue))
        return Jp6
    
    def Jminus(self,Jvalue):
        elements = J(Jvalue)
        Jminus = np.sqrt(elements.J*(elements.J+1)-elements.J_array*(elements.J_array-1))
        return Jminus
    
    def Jm(self,Jvalue):
        elements = J(Jvalue)
        Jminus = elements.Jminus(Jvalue)
        for i in range(elements.J_degeneracy-1):
            elements.J_matrix_1[i+1][i] = Jminus[i+1]
        return elements.J_matrix_1
    
    def Jm2(self,Jvalue):
        elements = J(Jvalue)
        Jm2 = np.dot(elements.Jm(Jvalue),elements.Jm(Jvalue))
        return Jm2
    
    def Jm3(self,Jvalue):
        elements = J(Jvalue)
        Jm3 = np.dot(elements.Jm(Jvalue),elements.Jm2(Jvalue))
        return Jm3
    
    def Jm4(self,Jvalue):
        elements = J(Jvalue)
        Jm4 = np.dot(elements.Jm(Jvalue),elements.Jm3(Jvalue))
        return Jm4
        
    def Jm5(self,Jvalue):
        elements = J(Jvalue)
        Jm5 = np.dot(elements.Jm(Jvalue),elements.Jm4(Jvalue))
        return Jm5
        
    def Jm6(self,Jvalue):
        elements = J(Jvalue)
        Jm6 = np.dot(elements.Jm(Jvalue),elements.Jm5(Jvalue))
        return Jm6
    
    def Ojm(self,Jvalue,j,m):
        
        element = J(Jvalue)
        Jz = element.Jz(Jvalue)
        Jp = element.Jp(Jvalue)
        Jp2 = element.Jp2(Jvalue)
        Jp3 = element.Jp3(Jvalue)
        Jp4 = element.Jp4(Jvalue)
        Jp5 = element.Jp5(Jvalue)
        Jp6 = element.Jp6(Jvalue)
        Jm = element.Jm(Jvalue)
        Jm2 = element.Jm2(Jvalue)
        Jm3 = element.Jm3(Jvalue)
        Jm4 = element.Jm4(Jvalue)
        Jm5 = element.Jm5(Jvalue)
        Jm6 = element.Jm6(Jvalue)
        J_matrix = element.J_tot
        i = element.i
        X = J_matrix*(J_matrix+1*i)
        
        if [j,m] == [0,0]:
           return np.zeros((int(2*J_matrix+i), int(2*J_matrix+i)))
        elif [j,m] == [1,0]:
           Ojm = Jz
        elif [j,m] == [1,1]:
           Ojm = 0.5 *(Jp + Jm)
        elif [j,m] == [1,-1]:
           Ojm = -0.5j *(Jp - Jm)
           
        elif [j,m] == [2,2]:
           Ojm = 0.5 *(Jp2 + Jm2)
        elif [j,m] == [2,1]:
           Ojm = 0.25*(np.dot(Jz,(Jp + Jm)) + np.dot((Jp + Jm),Jz))
        elif [j,m] == [2,0]:
           Ojm = 3*Jz**2 - X*i
        elif [j,m] == [2,-1]:
           Ojm = -0.25j*(np.dot(Jz,(Jp - Jm)) + np.dot((Jp - Jm),Jz))
        elif [j,m] == [2,-2]:
           Ojm = -0.5j *(Jp2 - Jm2)
           
        elif [j,m] == [3,3]:
           Ojm = 0.5 *(Jp3 + Jm3)
        elif [j,m] == [3,2]:
           Ojm = 0.25 *(np.dot((Jp2 + Jm2),Jz) + np.dot(Jz,(Jp2 + Jm2)))
        elif [j,m] == [3,1]:
           Ojm = 0.25*(np.dot((Jp + Jm),(5*Jz**2 - X - 0.5*i)) + np.dot((5*Jz**2 - X - 0.5*i),(Jp + Jm)))
        elif [j,m] == [3,0]:
           Ojm = 5*Jz**3 - (3*X-1)*Jz
        elif [j,m] == [3,-1]:
           Ojm = -0.25j*(np.dot((Jp - Jm),(5*Jz**2 - X - 0.5*i)) + np.dot((5*Jz**2 - X - 0.5*i),(Jp - Jm)))
        elif [j,m] == [3,-2]:
           Ojm = -0.25j*(np.dot(Jz,(Jp2 - Jm2)) + np.dot((Jp2 - Jm2),Jz))
        elif [j,m] == [3,-3]:
           Ojm = -0.5j *(Jp3 - Jm3)
    
        elif [j,m] == [4,4]:
           Ojm = 0.5 *(Jp4 + Jm4)
        elif [j,m] == [4,3]:
           Ojm = 0.25 *(np.dot((Jp3 + Jm3),Jz) + np.dot(Jz,(Jp3 + Jm3)))
        elif [j,m] == [4,2]:
           Ojm = 0.25 *(np.dot((Jp2 + Jm2),(7*Jz**2 -X -5*i)) + np.dot((7*Jz**2 -X -5*i),(Jp2 + Jm2)))
        elif [j,m] == [4,1]:
           Ojm = 0.25 *(np.dot((Jp + Jm),(7*Jz**3 -(3*X+1)*Jz)) + np.dot((7*Jz**3 -(3*X+1)*Jz),(Jp + Jm)))
        elif [j,m] == [4,0]:
           Ojm = 35*Jz**4 - (30*X -25*i)*Jz**2 + 3*X**2 - 6*X
        elif [j,m] == [4,-4]:
           Ojm = -0.5j *(Jp**4 - Jm**4)
        elif [j,m] == [4,-3]:
           Ojm = -0.25j *(np.dot((Jp3 - Jm3),Jz) + np.dot(Jz,(Jp3 - Jm3)))
        elif [j,m] == [4,-2]:
           Ojm = -0.25j *(np.dot((Jp2 - Jm2),(7*Jz**2 -X -5*i)) + np.dot((7*Jz**2 -X -5*i),(Jp2 - Jm2)))
        elif [j,m] == [4,-1]:
           Ojm = -0.25j *(np.dot((Jp - Jm),(7*Jz**3 -(3*X+1*i)*Jz)) + np.dot((7*Jz**3 -(3*X+1*i)*Jz),(Jp - Jm)))
        
        elif [j,m] == [5,5]:
           Ojm = (1/2)*(Jp5 + Jm5)
        elif [j,m] == [5,4]:
           Ojm = (1/4)*(np.dot((Jp4 + Jm4),Jz) + np.dot(Jz,(Jp4 + Jm4)))
        elif [j,m] == [5,3]:
           Ojm = (1/4)*(np.dot((Jp3 + Jm3),(9*Jz**2 - X - 33*i/2)) + np.dot((9*Jz**2 - X - 33*i/2),(Jp3 + Jm3)))
        elif [j,m] == [5,2]:
           Ojm = (1/4)*(np.dot((Jp2 + Jm2),(3*Jz**3 - (X - 6*i)*Jz)) + np.dot((3*Jz**3 - (X - 6*i)*Jz),(Jp2 + Jm2)))
        elif [j,m] == [5,1]:
           Ojm = (1/4)*(np.dot((Jp + Jm),(21*Jz**4 - 14*Jz**2*X + X**2 - X + 3*i/2)) + np.dot((21*Jz**4 - 14*Jz**2*X + X**2 - X + 3*i/2),(Jp + Jm)))
        elif [j,m] == [5,0]:
           Ojm = (63*Jz**5 - (70*X - 105*i)*Jz**3 + (15*X**2 - 50*X + 12*i)*Jz)
        elif [j,m] == [5,-1]:
           Ojm = (-1j/4)*(np.dot((Jp + Jm),(21*Jz**4 - 14*Jz**2*X + X**2 - X + 3*i/2)) + np.dot((21*Jz**4 - 14*Jz**2*X + X**2 - X + 3*i/2),(Jp + Jm)))
        elif [j,m] == [5,-2]:
           Ojm = (-1j/4)*(np.dot((Jp2 + Jm2),(3*Jz**3 - (X - 6*i)*Jz)) + np.dot((3*Jz**3 - (X - 6*i)*Jz),(Jp2 + Jm2)))
        elif [j,m] == [5,-3]:
           Ojm = (-1j/4)*(np.dot((Jp3 + Jm3),(9*Jz**2 - X - 33*i/2)) + np.dot((9*Jz**2 - X - 33*i/2),(Jp3 + Jm3))) 
        elif [j,m] == [5,-4]:
           Ojm = (-1j/4)*(np.dot((Jp4 + Jm4),Jz) + np.dot(Jz,(Jp4 + Jm4))) 
        elif [j,m] == [5,-5]:
           Ojm = (-1j/2)*(Jp5 + Jm5) 
        
        elif [j,m] == [6,6]:
           Ojm = 0.5*(Jp6 + Jm6)
        elif [j,m] == [6,5]:
           Ojm = 0.25*(np.dot((Jp5 + Jm5),Jz) + np.dot(Jz,(Jp5 + Jm5)))
        elif [j,m] == [6,4]:
           Ojm = 0.25*(np.dot((Jp4 + Jm4),(11*Jz**2 -X -38*i)) + np.dot((11*Jz**2 -X -38*i),(Jp4 + Jm4)))
        elif [j,m] == [6,3]:
           Ojm = 0.25*(np.dot((Jp3 + Jm3),(11*Jz**3 -(3*X+59*i)*Jz)) + np.dot((11*Jz**3 -(3*X+59*i)*Jz),(Jp3 + Jm3)))
        elif [j,m] == [6,2]:
           Ojm = 0.25*(np.dot((Jp2 + Jm2),(33*Jz**4 -(18*X+123*i)*Jz**2 +X**2 +10*X +102*i)) + np.dot((33*Jz**4 -(18*X+123*i)*Jz**2 +X**2 +10*X +102*i),(Jp2 + Jm2)))
        elif [j,m] == [6,1]:
           Ojm = 0.25*(np.dot((Jp +Jm),((33*Jz**5 -(30*X-15*i)*Jz**3) +(5*X**2 -10*X +12*i)*Jz)) + np.dot((33*Jz**5 -(30*X-15*i)*Jz**3 +(5*X**2 -10*X +12*i)*Jz),(Jp+ Jm)))
        elif [j,m] == [6,0]:
           Ojm = 231*Jz**6 - (315*X-735*i)*Jz**4 + (105*X**2 -525*X +294*i)*Jz**2 - 5*X**3 + 40*X**2 - 60*X
        elif [j,m] == [6,-6]:
           Ojm = -0.5j*(Jp6 - Jm6)
        elif [j,m] == [6,-5]:
           Ojm = -0.25j*(np.dot((Jp5 - Jm5),Jz) + np.dot(Jz,(Jp5 - Jm5)))
        elif [j,m] == [6,-4]:
           Ojm = -0.25j*(np.dot((Jp4 - Jm4),(11*Jz**2 -X -38*i)) + np.dot((11*Jz**2 -X -38*i),(Jp4 - Jm4)))
        elif [j,m] == [6,-3]:
           Ojm = -0.25j*(np.dot((Jp**3 - Jm**3),(11*Jz**3 -(3*X+59*i)*Jz)) + np.dot((11*Jz**3 -(3*X+59*i)*Jz),(Jp**3 - Jm**3)))
        elif [j,m] == [6,-2]:
           Ojm = -0.25j*(np.dot((Jp**2 - Jm**2),(33*Jz**4 -(18*X+123*i)*Jz**2 +X**2 +10*X +102*i)) + np.dot((33*Jz**4 -(18*X+123*i)*Jz**2 +X**2 +10*X +102*i),(Jp**2 - Jm**2)))
        elif [j,m] == [6,-1]:
           Ojm = -0.25j*(np.dot((Jp - Jm),(33*Jz**5 -(30*X-15*i)*Jz**3 +(5*X**2 -10*X +12*i)*Jz) + np.dot((33*Jz**5 -(30*X-15*i)*Jz**3 +(5*X**2 -10*X +12*i)*Jz),(Jp - Jm))))
        
        return Ojm
    
    
    def PC(self,ion,Jvalue,d,Z):
        
        
        e = 1.60217646e-19 # C
        e_0 = 8.854e-22 # C/VAngstrom
        a0 = 0.52917721067    #Bohr radius in \AA
        
        ### B_lm = -|e|*p_ml*gamma_lm*<r_l>*theta_l ###
        ### theta_2 = alpha ###
        ### theta_4 = beta ###
        ###theta_6 = gamma ###
        elements = J(Jvalue)
        x=[]
        for k in range(len(d)):
            x = np.append(x,d[k][0])
        #print(x)
        y=[]
        for k in range(len(d)):
            y = np.append(y,d[k][1])
        #print(y)
        z=[]
        for k in range(len(d)):
            z = np.append(z,d[k][2])
        #print(z)
        d_sp = elements.Spherical_position(d)
        r = d_sp[0,:]
        
        print('---------------- Ligand positions (Angstrom) ---------------')
        print('[x[i]    y[i]    z[i]]:')
        print(d)
        print("r[i]: " + str(r))
        
        r_l = elements.r_l(ion)
        r_2 = r_l[0]
        r_4 = r_l[1]
        r_6 = r_l[2]
        
        theta = elements.theta(ion)
        alpha = theta[4]
        beta = theta[5]
        gamma = theta[6]
        
        #r_2 = cef.r_l(L,S)[]
        
        ### http://www2.cpfs.mpg.de/~rotter/homepage_mcphase/manual/node131.html ####
        Zs22 = (1/4)*np.sqrt(15/np.pi)*(2*x*y/r**5)
        Zs21 = (1/2)*np.sqrt(15/np.pi)*(y*z/r**5)
        Z20 = (1/4)*np.sqrt(5/np.pi)*(3*z**2 - r**2)/r**5
        Zc21 = (1/2)*np.sqrt(15/np.pi)*(x*z/r**5)
        Zc22 = (1/4)*np.sqrt(15/np.pi)*((x**2 - y**2)/r**5)

        Zs44 = (3/16)*np.sqrt(35/np.pi)*(4*(x**3*y - x*y**3)/r**9)
        Zs43 = (3/8)*np.sqrt(70/np.pi)*((3*x**2*y - y**3)*z/r**9)
        Zs42 = (3/8)*np.sqrt(5/np.pi)*(2*x*y*(7*z**2 - r**2)/r**9)
        Zs41 = (3/4)*np.sqrt(5/(2*np.pi))*(y*z*(7*z**2 - 3*r**2)/r**9)
        Z40 =  (3/16)*np.sqrt(1/np.pi)*((35*z**4 - 30*z**2*r**2 + 3*r**4)/r**9)
        Zc41 = (3/4)*np.sqrt(5/(2*np.pi))*(x*z*(7*z**2 - 3*r**2)/r**9)
        Zc42 = (3/8)*np.sqrt(5/np.pi)*((x**2 - y**2)*(7*z**2 - r**2)/r**9)
        Zc43 = (3/8)*np.sqrt(70/np.pi)*((x**3 - 3*x*y**2)*z/r**9)
        Zc44 = (3/16)*np.sqrt(35/np.pi)*((x**4 - 6*x**2*y**2 + y**4)/r**9)

    
        Zs66= (231/64)*np.sqrt(26/(231*np.pi))*((6*x**5*y - 20*x**3*y**3 + 6*x*y**5)/r**13)
        Zs65 = np.sqrt(9009/(512*np.pi))*((5*x**4*y - 10*x**2*y**3 + y**5)*z/r**13)
        Zs64 = (21/32)*np.sqrt(13/(7*np.pi))*(4*(x**3*y - x*y**3)*(11*z**2 - r**2)/r**13)
        Zs63 = (1/32)*np.sqrt(2730/np.pi)*((3*x**2*y - y**3)*(11*z**3 - 3*x*r**2)/r**13)
        Zs62 = (1/64)*np.sqrt(2730/np.pi)*(2*x*y*(33*z**4 - 18*z**2*r**2 + r**4)/r**13)
        Zs61 = (1/8)*np.sqrt(273/(4*np.pi))*(y*z*(33*z**4 - 30*z**2*r**2 + 5*r**4)/r**13)
        Z60 = (1/32)*np.sqrt(13/np.pi)*((231*z**6 - 315*z**4*r**2 + 105*z**2*r**4 - 5*r**6)/r**13)
        Zc61 = (1/8)*np.sqrt(273/(4*np.pi))*(x*z*(33*z**4 - 30*z**2*r**2 + 5*r**4)/r**13)
        Zc62 = (1/64)*np.sqrt(2730/np.pi)*((x**2 - y**2)*(33*z**4 - 18*z**2*r**2 + r**4)/r**13)
        Zc63 = (1/32)*np.sqrt(2730/np.pi)*((x**3 - 3*x*y**2)*(11*z**3 - 3*z*r**2)/r**13)
        Zc64 = (21/32)*np.sqrt(13/(7*np.pi))*((x**4 - 6*x**2*y**2 + y**4)*(11*z**2 - r**2)/r**13)
        Zc65 = (np.sqrt(9009/(512*np.pi)))*((x**5 - 10*x**3*y**2 + 5*x*y**4)*z/r**13)
        Zc66 = (231/64)*np.sqrt(26/(231*np.pi))*((x**6 - 15*x**4*y**2 + 15*x**2*y**4 - y**6)/r**13)
        
        
        
        
        
        Bs22 = (-Z/5)*(1/4)*np.sqrt(15/np.pi)*r_2*alpha*np.sum(Zs22)*e/e_0*1000
        Bs21 = (-Z/5)*(1/4)*np.sqrt(15/np.pi)*r_2*alpha*np.sum(Zs21)*e/e_0*1000
        B20 = (-Z/5)*(1/4)*np.sqrt(5/np.pi)*r_2*alpha*np.sum(Z20)*e/e_0*1000
        Bc21 = (-Z/5)*(1/2)*np.sqrt(15/np.pi)*r_2*alpha*np.sum(Zc21)*e/e_0*1000
        Bc22 = (-Z/5)*(1/4)*np.sqrt(15/np.pi)*r_2*alpha*np.sum(Zc22)*e/e_0*1000

        B22 = Bs22 + Bc22
        B21 = Bs21 + Bc21
        
        Bs44 = (-Z/9)*(3/16)*np.sqrt(35/np.pi)*r_4*beta*np.sum(Zs44)*e/e_0*1000
        Bs43 = (-Z/9)*(3/8)*np.sqrt(70/np.pi)*r_4*beta*np.sum(Zs43)*e/e_0*1000
        Bs42 = (-Z/9)*(3/8)*np.sqrt(5/np.pi)*r_4*beta*np.sum(Zs42)*e/e_0*1000
        Bs41 = (-Z/9)*(3/4)*np.sqrt(5/(2*np.pi))*r_4*beta*np.sum(Zs41)*e/e_0*1000
        B40 = (-Z/9)*(3/16)*np.sqrt(1/np.pi)*r_4*beta*np.sum(Z40)*e/e_0*1000
        Bc41 = (-Z/9)*(3/4)*np.sqrt(5/(2*np.pi))*r_4*beta*np.sum(Zc41)*e/e_0*1000
        Bc42 = (-Z/9)*(3/8)*np.sqrt(5/np.pi)*r_4*beta*np.sum(Zc42)*e/e_0*1000
        Bc43 = (-Z/9)*(3/8)*np.sqrt(70/np.pi)*r_4*beta*np.sum(Zc43)*e/e_0*1000
        Bc44 = (-Z/9)*(3/16)*np.sqrt(35/np.pi)*r_4*beta*np.sum(Zc44)*e/e_0*1000

        B44 = Bs44 + Bc44
        B43 = Bs43 + Bc43
        B42 = Bs42 + Bc42
        B41 = Bs41 + Bc41
        
    
        Bs66 = (-Z/13)*(231/64)*np.sqrt(26/(231*np.pi))*r_6*gamma*np.sum(Zs66)*e/e_0*1000
        Bs65 = (-Z/13)*np.sqrt(9009/(512*np.pi))*r_6*gamma*np.sum(Zs65)*e/e_0*1000
        Bs64 = (-Z/13)*(21/32)*np.sqrt(13/(7*np.pi))*r_6*gamma*np.sum(Zs64)*e/e_0*1000
        Bs63 = (-Z/13)*(1/32)*np.sqrt(2730/np.pi)*r_6*gamma*np.sum(Zs63)*e/e_0*1000
        Bs62 = (-Z/13)*(1/64)*np.sqrt(2730/np.pi)*r_6*gamma*np.sum(Zs62)*e/e_0*1000
        Bs61 = (-Z/13)*(1/8)*np.sqrt(273/(4*np.pi))*r_6*gamma*np.sum(Zs61)*e/e_0*1000
        B60 = (-Z/13)*(1/32)*np.sqrt(13/np.pi)*r_6*gamma*np.sum(Z60)*e/e_0*1000
        Bc61 = (-Z/13)*(1/8)*np.sqrt(273/(4*np.pi))*r_6*gamma*np.sum(Zc61)*e/e_0*1000
        Bc62 = (-Z/13)*(1/64)*np.sqrt(2730/np.pi)*r_6*gamma*np.sum(Zc62)*e/e_0*1000
        Bc63 = (-Z/13)*(1/32)*np.sqrt(2730/np.pi)*r_6*gamma*np.sum(Zc63)*e/e_0*1000
        Bc64 = (-Z/13)*(21/32)*np.sqrt(13/(7*np.pi))*r_6*gamma*np.sum(Zc64)*e/e_0*1000
        Bc65 = (-Z/13)*(np.sqrt(9009/(512*np.pi)))*r_6*gamma*np.sum(Zc65)*e/e_0*1000
        Bc66 = (-Z/13)*(231/64)*np.sqrt(26/(231*np.pi))*r_6*gamma*np.sum(Zc66)*e/e_0*1000
        
        B66 = Bs66 + Bc66
        B65 = Bs65 + Bc65
        B64 = Bs64 + Bc64
        B63 = Bs63 + Bc63
        B62 = Bs62 + Bc62
        B61 = Bs61 + Bc61
        
        print("-------- Point Charge Modeling (meV) --------")
        print("B20: " + str(B20))
        print("B21: " + str(B21))
        print("B22: " + str(B22))
        print("B40: " + str(B40))
        print("B41: " + str(B41))
        print("B42: " + str(B42))
        print("B43: " + str(B43))
        print("B44: " + str(B44))
        print("B60: " + str(B60))
        print("B61: " + str(B61))
        print("B62: " + str(B62))
        print("B63: " + str(B63))
        print("B64: " + str(B64))
        print("B65: " + str(B65))
        print("B66: " + str(B66))
        
        B = np.array([B20, B21, B22, B40, B41, B42, B43, B44, B60, B61, B62, B63, B64, B65, B66])
        
        return B
    
    def Spherical_position(self,xyz_position):
        r = np.sqrt(xyz_position[:,0]**2 + xyz_position[:,1]**2 + xyz_position[:,2]**2)
        theta = np.arctan2(xyz_position[:,1],xyz_position[:,0])
        phi = np.arccos(xyz_position[:,2]/r)
        
        print('r:[], theta:[], phi:[]')
        print(r,theta,phi)
        # r = [:,0]
        # phi = [:,1]
        # theta = [:,2]
        
        sp_position = np.array([r,theta,phi])
    
        return sp_position
    
    def Diag(self,H,printfunction=True):
        E = np.linalg.eigh(H)
        E_val = E[0]
        E_state = E[1].T

        E_excitation = (E_val-E_val[0])
        
        print(E_val)
        print(E_excitation)
            
        fig = plt.figure()
            
        ax1 = fig.add_subplot(221)
        ax1.eventplot(E_val,orientation='vertical',linelength=0.05,linestyles='solid',colors='r')
        ax1.tick_params(axis='x',which='both',bottom=False,top=False)
        ax1.set_xlabel('Crystal Field Splitting')
        ax1.set_ylabel('Energy (meV)')
        ax1.set_xticklabels('')
        ax1.set_title('Point Charge Model')
        ax1.set_xlim([0.95,1.05])
            
        ax2 = fig.add_subplot(222)
        ax2.eventplot(E_excitation,orientation='horizontal',linelength=0.05,linestyles='solid',colors='b')
        ax2.tick_params(axis='y',which='both',bottom=False,top=False)
        ax2.set_xlabel('Energy (meV)')
        ax2.set_ylabel('Arb. Units')
        ax2.set_yticklabels('')
        ax2.set_yticks([])
        ax2.set_title('Crystal Field excitations')
        ax2.set_ylim([0.95,1.05])
            
        plt.show()
        
        return E_val,E_excitation,E_state
    
    def MagMoment_z(self,ion,Jvalue,H):
        print('\n')
        print('--------- Magnetic Properties (z-axis) --------')
        elements = J(Jvalue)
        Jz = elements.Jz(Jvalue)
        mu_B = (5.7883818012e-5)*1000 #meV
        g = elements.theta(ion)
        gj = g[3]
        B = 0.1
        H = H + mu_B*(Jz)*B
        H_e,H_v = np.linalg.eigh(H)
        E_final = H_e - H_e[0]
        H_vt = H_v.T
        print(E_final)
        print('----------------- Mz ---------------')
        Jz_m = np.zeros((elements.J_degeneracy))
        for k in range(elements.J_degeneracy):
            Jz_m[k] = gj*np.dot(np.conjugate(H_vt[k]),np.dot(Jz,H_vt[k]))
        print(Jz_m)
        
        return Jz_m
    
    def LandeFactor(self,Jvalue,L,S):
        g = 1 + (Jvalue*(Jvalue+1) + S*(S+1) - L*(L+1))/(2*Jvalue*(Jvalue+1))
        return g
        
    
        
    def r_l(self,ion):
        
        #From Abragam and Bleaney
        
        ##### ion2+ = [<r_2>, <r_4>, <r_6>]
        
        if ion == "Ce3":
            r_l = np.array([0.3666, 0.3108, 0.5119])
            return r_l
        elif ion == "Pr3":
            r_l = np.array([0.3350, 0.2614, 0.4030])
        elif ion == "Nd3":
            r_l = np.array([0.3120, 0.2282, 0.3300])
        elif ion == "Pm3":
            r_l = np.array([0.2899, 0.1991, 0.2755])
        elif ion == "Sm3":
            r_l = np.array([0.2728, 0.1772, 0.2317])
        elif ion == "Eu3":
            r_l = np.array([0.2569, 0.1584, 0.1985])
        elif ion == "Gd3":
            r_l = np.array([0.2428, 0.1427, 0.1720])
        elif ion == "Tb3":
            r_l = np.array([0.2302, 0.1295, 0.1505])
        elif ion == "Dy3":
            r_l = np.array([0.2188, 0.1180, 0.1328])
        elif ion == "Ho3":
            r_l = np.array([0.2085, 0.1081, 0.1181])
        elif ion == "Er3":
            r_l = np.array([0.1991, 0.0996, 0.1058])
        elif ion == "Tm3":
            r_l = np.array([0.1905, 0.0921, 0.0953])
        elif ion == "Yb3":
            r_l = np.array([0.1826, 0.0854, 0.0863])
        elif ion == "U4":
            r_l = np.array([0.5718, 0.5985, 1.0491])
        elif ion == "U3":
            r_l = np.array([0.6569, 0.8552, 1.9882])
        elif ion == "Np4":
            r_l = np.array([0.5276, 0.5100, 0.8300])
        elif ion == "Nd2":
            r_l = np.array([0.3898, 0.4191, 0.9980])
        elif ion == "Sm2":
            r_l = np.array([0.3352, 0.3028, 0.6271])
        elif ion == "Eu2":
            r_l = np.array([0.3075, 0.2641, 0.5178])
        elif ion == "Gd2":
            r_l = np.array([0.2879, 0.2333, 0.4359])
        elif ion == "Tb2":
            r_l = np.array([0.2711, 0.2082, 0.3729])
        elif ion == "Dy2":
            r_l = np.array([0.2557, 0.1875, 0.3235])
        elif ion == "Ho2":
            r_l = np.array([0.2425, 0.1701, 0.2837])
        elif ion == "Er2":
            r_l = np.array([0.2307, 0.1552, 0.2514])
        elif ion == "Tm2":
            r_l = np.array([0.2198, 0.1426, 0.2249])
        elif ion == "Yb2":
            r_l = np.array([0.2100, 0.1315, 0.2027])
        
        
        return r_l
 
    def theta(self,ion):
        
        #From Abragam and Bleaney
        
        ##### ion2+ = [l, s, j, gj, <alpha>, <beta>, <gamma>]
        
        if ion == "Ce3":
            theta = np.array([3, 1/2, 5/2, 6/7, -5.7143e-2, 63.4921e-4, 0.0000])
            return theta
        elif ion == "Pr3":
            theta = np.array([5, 1, 4, 4/5, -2.1010, -2.1010e-2, -7.3462e-4, 60.9940e-6])
        elif ion == "Ce2":
            theta = np.array([5, 1, 4, 4/5, -2.1010, -2.1010e-2, -7.3462e-4, 60.9940e-6])
        elif ion == "U4":
            theta = np.array([5, 1, 4, 4/5, -2.1010, -2.1010e-2, -7.3462e-4, 60.9940e-6])
        elif ion == "Nd3":
            theta = np.array([6, 3/2, 9/2, 8/11, -0.6428e-2,-2.9111e-4, -37.9880e-6])
        elif ion == "Pr2":
            theta = np.array([6, 3/2, 9/2, 8/11, -0.6428e-2,-2.9111e-4, -37.9880e-6])
        elif ion == "U3":
            theta = np.array([6, 3/2, 9/2, 8/11, -0.6428e-2,-2.9111e-4, -37.9880e-6])
        elif ion == "Np4":
            theta = np.array([6, 3/2, 9/2, 8/11, -0.6428e-2,-2.9111e-4, -37.9880e-6])
        elif ion == "Pm3":
            theta = np.array([6, 2, 4, 3/5, 0.7713, 0.7713e-2, 4.0755e-4, 60.7807e-6])
        elif ion == "Nd2":
            theta = np.array([6, 2, 4, 3/5, 0.7713, 0.7713e-2, 4.0755e-4, 60.7807e-6])
        elif ion == "Sm3":
            theta = np.array([5, 5/2, 5/2, 2/7, 4.1270e-2, 25.0120e-4, 0.0000])
        elif ion == "Pm2":
            theta = np.array([5, 5/2, 5/2, 2/7, 4.1270e-2, 25.0120e-4, 0.0000])
        elif ion == "Eu3":
            theta = np.array([3, 3, 0, 0, 0.0000, 0.0000, 0.0000])
        elif ion == "Sm2":
            theta = np.array([3, 3, 0, 0, 0.0000, 0.0000, 0.0000])
        elif ion == "Gd3":
            theta = np.array([0, 7/2, 7/2, 2, 0.0000, 0.0000, 0.0000])
        elif ion == "Eu2":
            theta = np.array([0, 7/2, 7/2, 2, 0.0000, 0.0000, 0.0000])
        elif ion == "Tb3":
            theta = np.array([3, 3, 6, 3/2, -1.0101e-2, 1.2244e-4, -1.1212e-6])
        elif ion == "Gd2":
            theta = np.array([3, 3, 6, 3/2, -1.0101e-2, 1.2244e-4, -1.1212e-6])
        elif ion == "Dy3":
            theta = np.array([5, 5/2, 15/2, 4/3, -0.6349e-2, -0.5920e-4, 1.0350e-6])
        elif ion == "Tb2":
            theta = np.array([5, 5/2, 15/2, 4/3, -0.6349e-2, -0.5920e-4, 1.0350e-6])
        elif ion == "Ho3":
            theta = np.array([6, 2, 8, 5/4, -0.2222e-2, -0.3330e-4, -1.2937e-6])
        elif ion == "Dy2":
            theta = np.array([6, 2, 8, 5/4, -0.2222e-2, -0.3330e-4, -1.2937e-6])
        elif ion == "Er3":
            theta = np.array([6, 3/2, 15/2, 6/5, 0.2540e-2, 0.4440e-4, 2.0699e-6])
        elif ion == "Ho2":
            theta = np.array([6, 3/2, 15/2, 6/5, 0.2540e-2, 0.4440e-4, 2.0699e-6])
        elif ion == "Tm3":
            theta = np.array([5, 1, 6, 7/6, 1.0101e-2, 1.6325e-4, -5.6061e-6])
        elif ion == "Er2":
            theta = np.array([5, 1, 6, 7/6, 1.0101e-2, 1.6325e-4, -5.6061e-6])
        elif ion == "Yb3":
            theta = np.array([3, 1/2, 7/2, 8/7, 3.1746e-2, -17.3160e-4, 148.0001e-6])
        elif ion == "Tm2":
            theta = np.array([3, 1/2, 7/2, 8/7, 3.1746e-2, -17.3160e-4, 148.0001e-6])
        
        return theta
"""

Dq = 302.8 # cm-1
A2 = -33.3 # cm-1
A4 = -399.4 # cm-1

B02 = A2/3
print(B02)
B04 = A4/12 + Dq/18
print(B04)
B34 = -20*Dq/(9*np.sqrt(2))
print(B34)



#This is the correct units

ahc = 1.43996e4  #Constant to get the energy in units of meV = alpha*hbar*c
a0 = 0.52917721067    #Bohr radius in \AA
r_2 = 0.1826
r_4 = 0.0854
r_6 = 0.0863
alpha = 3.1746e-2
beta = -17.3160e-4
gamma = 148e-6
e = 1.60217646e-19
e_0 = 8.854e-22
r_test = np.array([4,4])

Z20 = (-Z/r**3)*(1/4)*np.sqrt(5/np.pi)*((3*z**2 - r**2)/r**2)

gamma02 = -(0.2/5)*(1/4)*np.sqrt(5/np.pi)*(3*r_test**2 - r_test**2)/r_test**5

B02 = -(1/4)*np.sqrt(5/np.pi)*r_2*alpha*np.sum(gamma02)*e/e_0*1000

print(B02)

gamma04 = -(0.2/9)*2*(3/16)*np.sqrt(1/np.pi)*(35*4**4 - 30*4**2*4**2 + 3*4**4)/4**9

B04 = -(3/16)*np.sqrt(1/np.pi)*r_4*beta*gamma04*e/e_0*1000

print(B04)

gamma06 = -(0.2/13)*2*(1/32)*np.sqrt(13/np.pi)*(231*4**6 -315*4**4*4**2 + 105*4**2*4**4 - 5*4**6)/4**13

B06 = - (1/32)*np.sqrt(13/np.pi)*r_6*gamma*gamma06*e/e_0*1000

print(B06)

"""


