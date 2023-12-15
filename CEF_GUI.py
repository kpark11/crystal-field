# -*- coding: utf-8 -*-
"""
Created on Wed Dec 29 17:18:53 2021

@author: brian
"""
import sys
cef_dir = r'C:\Users\kit\OneDrive - University of Tennessee\Desktop\Research\Python program\CEF calculations\CEF'
sys.path.append(cef_dir)
import Crystal_Field_Calculations as cef
from tkinter import *
from tkinter.ttk import Combobox
from tkinter.filedialog import asksaveasfilename, askopenfilename, askdirectory
from tkinter import scrolledtext
import numpy as np
import matplotlib.pyplot as plt

def r_l():
    f = open('r_l.txt','r')
    lines = f.readlines()
    f.close()
    i = []
    for line in lines:
        if "#" in line:
            pass
        elif len(line) == 1:
            pass
        else:
            i.append(line.split(sep=','))
    return i
rl = r_l()
ions = []
for i in rl:
    ions.append(i[0])
    
def theta():
    f = open('theta.txt')
    lines = f.readlines()
    f.close()
    i = []
    for line in lines:
        if '#' in line:
            pass
        elif len(line) == 1 or len(line) == 9 or len(line) == 10:
            pass
        else:
            i.append(line.split(sep=','))
    return i
theta_ion = theta()


def pos(positions):
    positions = positions.split('\n')
    x = np.array([])
    y = np.array([])
    z = np.array([])
    for position_xyz in positions:
        print(position_xyz)
        pp = position_xyz.split()
        print(pp)
        
        i = pp[1].split(':')
        x.append(float(i[1]))

        i = pp[2].split(':')
        y.append(float(i[1]))

        i = pp[3].split(':')
        z.append(float(i[1]))

    print(p)
    return p


class MyWindow:
    def __init__(self, win):
        self.ion_l = Label(win, text = 'Choose an ion with correct oxidation number')
        self.ion_l.place(x=10, y=10)
        
        self.ion = Combobox(win)
        self.ion['values'] = ions
        self.ion.place(x=10, y=30)
        
        self.L_l = Label(win, text = 'L: ')
        self.L_l.place(x=10, y=55)
        
        self.L = Entry(win,width=3)
        self.L.place(x=10, y=75)
        
        self.S_l = Label(win, text = 'S: ')
        self.S_l.place(x=50, y=55)
        
        self.S = Entry(win,width=3)
        self.S.place(x=50, y=75)
        
        self.Z_l = Label(win, text = 'Z for ligands (e.g. 2 for oxygen)')
        self.Z_l.place(x=10, y=95)
        
        self.Z = Entry(win,width=3)
        self.Z.place(x=10, y=115)
    
        self.num_lig_label = Label(win,text='Number of ligands (e.g. 4 for tetrahedral, 6 for octahedral')
        self.num_lig_label.place(x=10,y=140)
        self.num_lig = Spinbox(window, from_=0, to=10, width=5)
        self.num_lig.place(x=10,y=160)
        
        def insert_lig():
            self.lig_pos.delete('1.0', END)
            ligands = int(self.num_lig.get())
            self.lig_pos.insert(END, 'ion:    X:    Y:    Z:    \n')
            for i in range(ligands):
                self.lig_pos.insert(END, 'ligand' + str(i) + ':  ' + '  X:    Y:    Z:    \n')
            
        self.insert_lig = Button(win,text='insert',command=insert_lig)
        self.insert_lig.place(x=10,y=180)
        
        self.lig_pos = scrolledtext.ScrolledText(window,width=60,height=5)
        self.lig_pos.place(x=10,y=220)
        
        
        def initialize():
            try: 
                L_value = int(self.L.get())
                S_value = int(self.S.get())
                
                ligands = int(self.num_lig.get())
                lig_positions = self.lig_pos.get('1.0',END)
                lig_positions = pos(lig_positions)
                CEF = cef.LS(L_value,S_value)
                
                
                
                
                
                
                self.results1.insert(END,"-------------------------------\n")
                self.results1.insert(END, 'L: ' + self.L.get() + '\n')
                self.results1.insert(END, 'S: ' + self.S.get() + '\n')
                self.results1.insert(END, 'ligand numbers: ' + self.num_lig.get() + '\n')
                self.results1.insert(END, 'ligand positions:\n')
                self.results1.insert(END, self.lig_pos.get('1.0',END))
                
                
            except:
                self.results1.insert(END,"-------------------------------\n")
                self.results1.insert(END,'something is wrong\n')
                print('something is wrong')
                
        self.btn = Button(win,text='initialize',command=initialize)
        self.btn.place(x=10,y=310)
        
        
        self.results = Label(win,text='results:')
        self.results.place(x=10,y=550)
        self.results1 = scrolledtext.ScrolledText(window,width=60,height=5)
        self.results1.place(x=10,y=570)
        
        
        def clear_results():
            self.results1.delete('1.0', END)
        self.btn = Button(win,text='clear',command=clear_results)
        self.btn.place(x=10,y=660)
        
        
        
    '''
    def add(self):
        self.t3.delete(0, 'end')
        num1=int(self.t1.get())
        num2=int(self.t2.get())
        result=num1+num2
        self.t3.insert(END, str(result))
    def sub(self, event):
        self.t3.delete(0, 'end')
        num1=int(self.t1.get())
        num2=int(self.t2.get())
        result=num1-num2
        self.t3.insert(END, str(result))
    '''
    

window=Tk()
mywin=MyWindow(window)
window.title('Crystal Field calculation')
window.geometry("800x700+10+10")
window.mainloop()



