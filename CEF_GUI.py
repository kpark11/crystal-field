# -*- coding: utf-8 -*-
"""
Created on Wed Dec 29 17:18:53 2021

@author: brian
"""
import sys
import os
cef_dir = r'C:\Users\brian\OneDrive - University of Tennessee\Desktop\Research\Python program\CEF calculations\CEF'
sys.path.append(cef_dir)
import Crystal_Field_Calculations as cef
import tkinter as tk
from tkinter.ttk import *
from tkinter.ttk import Combobox, Progressbar, Menubutton, Style, Separator
from tkinter.filedialog import asksaveasfilename, askopenfilename, askdirectory
from tkinter import scrolledtext, Menu, Canvas
from ttkthemes import ThemedTk,ThemedStyle

from PIL import ImageTk, Image, ImageDraw
import cv2

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.figure import Figure 
from matplotlib.backends.backend_tkagg import (FigureCanvasTkAgg,  
NavigationToolbar2Tk) 




def convert_to_float(frac_str):
    try:
        return float(frac_str)
    except ValueError:
        try:
            num, denom = frac_str.split('/')
        except ValueError:
            return None
        try:
            leading, num = num.split(' ')
        except ValueError:
            return float(num) / float(denom)        
        if float(leading) < 0:
            sign_mult = -1
        else:
            sign_mult = 1
        return float(leading) + sign_mult * (float(num) / float(denom))

def remove_values_from_list(the_list, val):
   return [value for value in the_list if value != val]


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

class StartUp:
    def __init__(self,t):

        self.win = ThemedTk(theme=t)
        color = Style().lookup("TFrame", "background", default="black")
        self.win.configure(bg=color)
        
        self.win.title("Initializing Crystal Field Calculation Program")
        self.win.eval('tk::PlaceWindow . center')
        
        self.win.attributes('-disabled', True)
        
        self.frame = Frame(self.win)
        self.frame.grid(row=0,column=0,padx= 200,pady=100,stick="nsew")
        
        #canvas = Canvas(self.win, height=HEIGHT, width=WIDTH)
        #canvas.pack()
        
        #logo = tk.PhotoImage(file='SphericalP_logo.png')
        #logolabel = tk.Label(frame, image=logo, bg="black")
        #logolabel.place(relx=0, rely=0)
        progStatus = tk.StringVar()
        progStatus.set("Initializing...")
        
        label_1 = Label(self.frame, textvariable=progStatus)
        label_1.grid(row=0,column=0)

        progress = Progressbar(self.frame, orient = tk.HORIZONTAL, length = 100, mode = 'determinate')
        progress.grid(row=1,column=0)

        progress['value'] = 0
        self.win.update()
        self.win.after(1000, progStatus.set("Loading Calibration"))
        progress['value'] = 40
        self.win.update()
        self.win.after(1000, progStatus.set("Ready"))
        progress['value'] = 80
        self.win.update()
        self.win.after(500, self.win.destroy)
        progress['value'] = 100
        self.win.update()
    


def pos(positions):
    positions = positions.split('\n')
    
    positions = remove_values_from_list(positions,'')
    positions = remove_values_from_list(positions,' ')
    Pos =[]
    print(positions)
    for position_xyz in positions:
        pp = position_xyz.split()
        xx = pp[1].split(':')
        x = float(xx[1])
        yy = pp[2].split(':')
        y = float(yy[1])
        zz = pp[3].split(':')
        z = float(zz[1])
        Pos.append([x,y,z])
    Pos = np.array(Pos)
    for i in range(len(Pos)-1):
        Pos[i+1] = Pos[i+1] - Pos[0]
    Pos = Pos[1:]

    return Pos

def typeB(B):
    return "-------- Point Charge Modeling (meV) --------\n" +\
    "B20: " + str(B[0]) + ' \n' +\
    "B21: " + str(B[1]) + ' \n' +\
    "B22: " + str(B[2]) + ' \n' +\
    "B40: " + str(B[3]) + ' \n' +\
    "B41: " + str(B[4]) + ' \n' +\
    "B42: " + str(B[5]) + ' \n' +\
    "B43: " + str(B[6]) + ' \n' +\
    "B44: " + str(B[7]) + ' \n'
    

def bar(): 
    import time 
    progress['value'] = 20
    root.update_idletasks() 
    time.sleep(1) 
  
    progress['value'] = 40
    root.update_idletasks() 
    time.sleep(1) 
  
    progress['value'] = 50
    root.update_idletasks() 
    time.sleep(1) 
  
    progress['value'] = 60
    root.update_idletasks() 
    time.sleep(1) 
  
    progress['value'] = 80
    root.update_idletasks() 
    time.sleep(1) 
    progress['value'] = 100
    


class MyWindow:
    def __init__(self, win):
        win.after(500)
        
        def hello():
            x = 0
        
        menubar = Menu(win)
        # create a pulldown menu, and add it to the menu bar
        filemenu = Menu(menubar, tearoff=0)
        
        filemenu.config(bg=color,fg='black',activebackground='black',activeforeground='white',relief=tk.FLAT)

        
        filemenu.add_command(label="Open Existing Work Space", command=hello)
        filemenu.add_command(label="Save Work Space", command=hello)
        filemenu.add_separator()
        filemenu.add_command(label="Create New Work Space", command=hello)
        menubar.add_cascade(label="Work Space", menu=filemenu)
        
        # create more pulldown menus ()
        editmenu = Menu(menubar, tearoff=0)
        #editmenu.add_command(label="Calibrate", command=lambda: menuAction(0))
        editmenu.add_command(label="Simulation Mode", command=hello)
        #editmenu.add_command(label="Global Variables", command=lambda: menuAction(1))
        editmenu.add_command(label="Log", command=hello)
        menubar.add_cascade(label="Diagnostics", menu=editmenu)

        helpmenu = Menu(menubar, tearoff=0)
        helpmenu.add_command(label="Operation Manual", command=hello)
        #helpmenu.add_command(label="About", command=snpAboutL)
        menubar.add_cascade(label="Help", menu=helpmenu)
        
        filemenu.add_command(label="Open Existing Work Space")
        # display the menu
        win.config(menu=menubar)
        

        
        
        
        ####################################################################################
        
        frame = Frame(win)
        frame.grid(row=0,column=0)
        
        ####################################################################################
        
        frame1 = Frame(frame)
        frame1.grid(row=0,column=0,padx=5)

        
       



        self.ion_l = Label(frame1, text = 'Choose an ion with correct oxidation number')
        #self.ion_l.place(x=10, y=10)
        #self.ion_l.pack(side = LEFT,padx=10)
        self.ion_l.grid(row=0,column=0,pady=5)
        
        self.ion = Combobox(frame1)
        self.ion['values'] = ions
        #self.ion.place(x=10, y=30)
        #self.ion.pack(side=LEFT,padx=10)
        self.ion.grid(row=0,column=1,pady=5)
        
        self.L_l = Label(frame1, text = 'L: ')
        #self.L_l.place(x=10, y=55)
        #self.L_l.pack(side=LEFT,padx=10)
        self.L_l.grid(row=1,column=0,pady=5)
        
        self.L = Entry(frame1,width=3)
        #self.L.place(x=10, y=75)
        #self.L.pack(side=LEFT,padx=10)
        self.L.grid(row=1,column=1,pady=5)
        
        self.S_l = Label(frame1, text = 'S: ')
        self.S_l.place(x=50, y=55)
        #self.S_l.pack(side=LEFT,padx=10)
        self.S_l.grid(row=2,column=0,pady=5)
        
        
        self.S = Entry(frame1,width=3)
        self.S.place(x=50, y=75)
        #self.S.pack(side=LEFT,padx=10)
        self.S.grid(row=2,column=1,pady=5)
        
        
        self.Z_l = Label(frame1, text = 'Z for ligands (e.g. 2 for oxygen)')
        #self.Z_l.place(x=10, y=95)
        self.Z_l.grid(row=3,column=0,pady=5)
        
        self.Z = Entry(frame1,width=3)
        self.Z.place(x=10, y=115)
        self.Z.grid(row=3,column=1,pady=5)
    
        self.num_lig_label = Label(frame1,text='Number of ligands (e.g. 4 for tetrahedral, 6 for octahedral')
        #self.num_lig_label.place(x=10,y=140)
        self.num_lig_label.grid(row=0,column=2,pady=5)
        
        self.num_lig = Spinbox(frame1, from_=0, to=10, width=5)
        #self.num_lig.place(x=10,y=160)
        self.num_lig.grid(row=1,column=2,pady=5)
        
        def insert_lig():
            self.lig_pos.delete('1.0', tk.END)
            ligands = int(self.num_lig.get())
            self.lig_pos.insert(tk.END, 'ion:    X:    Y:    Z:    \n')
            for i in range(ligands):
                self.lig_pos.insert(tk.END, 'ligand' + str(i) + ':  ' + '  X:    Y:    Z:    \n')
                
        self.insert_lig = Button(frame1,text='insert',command=insert_lig)
        #self.insert_lig.place(x=10,y=180)
        self.insert_lig.grid(row=2,column=2,pady=5)
        
        self.lig_pos = scrolledtext.ScrolledText(frame1,width=50,height=5)
        #self.lig_pos.place(x=10,y=220)
        self.lig_pos.grid(row=3,column=2,pady=10)
            

        

        
        

        def initialize():
            self.progress["value"] = 0
            win.update()
            try: 
                self.progress['value'] = 5
                win.update()
                
                L_value = int(self.L.get())
                S_value = convert_to_float(self.S.get())

                CEF = cef.LS(L_value,S_value)
                ligands = int(self.num_lig.get())
                lig_positions = self.lig_pos.get('1.0',tk.END)
                lig_positions = pos(lig_positions)
                
                ion_value = self.ion.get().replace('"', '')
                ion_value = ion_value.replace('"','')
                Z_value = int(self.Z.get())
                self.progress['value'] = 20
                win.update()
                
                ######################################################################
                B = CEF.PC(ion_value,L_value,S_value,lig_positions,Z_value)
                
                O20 = CEF.Olm(L_value,S_value,2,0)
                O21 = CEF.Olm(L_value,S_value,2,1)
                O22 = CEF.Olm(L_value,S_value,2,2)
                O40 = CEF.Olm(L_value,S_value,4,0)
                O41 = CEF.Olm(L_value,S_value,4,1)
                O42 = CEF.Olm(L_value,S_value,4,2)
                O43 = CEF.Olm(L_value,S_value,4,3)
                O44 = CEF.Olm(L_value,S_value,4,4)
                O2m1 = CEF.Olm(L_value,S_value,2,-1)
                O2m2 = CEF.Olm(L_value,S_value,2,-2)
                O4m1 = CEF.Olm(L_value,S_value,4,-1)
                O4m2 = CEF.Olm(L_value,S_value,4,-2)
                O4m3 = CEF.Olm(L_value,S_value,4,-3)
                O4m4 = CEF.Olm(L_value,S_value,4,-4)
                ######################################################################
                Hcf = B[0]*O20 + B[1]*O21 + B[2]*O22 + B[3]*O40 + B[4]*O41 + B[5]*O42 + B[6]*O43 + B[7]*O44
                Ecf_val,Ecf_val_excitation,H_cf_vt = CEF.Diag(Hcf)
                ######################################################################
                SO_matrix,SO_val = CEF.SO(ion_value,L_value,S_value)
                Hcf_so = Hcf + SO_matrix
                Ecf_so_val,Ecf_so_val_excitation,H_cf_so_vt = CEF.Diag(Hcf_so)
                ######################################################################
                Hss,SS_val = CEF.SS(ion_value,L_value,S_value)
                Hcf_so_ss = Hcf_so + Hss
                Ecf_so_ss_val,Ecf_so_ss_val_excitation,H_cf_so_ss_vt = CEF.Diag(Hcf_so_ss)
                ######################################################################
                Hg,g_val = CEF.Molecular_Field_Sz(ion_value,L_value,S_value)
                Hcf_so_ss_g = Hcf_so_ss + Hg
                Ecf_so_ss_g_val,Ecf_so_ss_g_val_excitation,H_cf_so_ss_g_vt = CEF.Diag(Hcf_so_ss_g)
                ######################################################################
                def PlotEnergies(E_val,E_excitation,canvas,title,savePlot=False):
                    global output, fig

                    fig = Figure()
                    
                    ax1 = fig.add_subplot(221)
                    ax1.eventplot(E_val,orientation='vertical',linelength=0.05,linestyles='solid',colors='r')
                    ax1.tick_params(axis='x',which='both',bottom=False,top=False)
                    ax1.set_xlabel(title)
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
                    fig.tight_layout()
                    
                    output = FigureCanvasTkAgg(fig, master=canvas)
                    
                    output.draw()
                    output.get_tk_widget().pack()

                
                
                def PopUp():
                    self.popup = tk.Tk()
                    popup_style = ThemedStyle(self.popup)
                    popup_style.theme_use(t)
                    self.popup.configure(bg=color)
                    self.popup.wm_title("Plots")
                    
                    self.frame_popup = Frame(self.popup,width=1200,height=700)
                    self.frame_popup.grid(row=0,column=0)
                    
                    self.cv = Canvas(self.frame_popup)
                    self.cv.grid(row=0,column=0)
                    self.cv1 = Canvas(self.frame_popup)
                    self.cv1.grid(row=0,column=1)
                    self.cv2 = Canvas(self.frame_popup)
                    self.cv2.grid(row=1,column=0)
                    self.cv3 = Canvas(self.frame_popup)
                    self.cv3.grid(row=1,column=1)
                    
                    self.frame_popup1 = Frame(self.popup)
                    self.frame_popup1.grid(row=0,column=1)
                    self.popup_btn = Button(self.frame_popup1,text='Save',command=SavePlot)
                    self.popup_btn.grid(row=0,column=1)
                    self.popup_btn1 = Button(self.frame_popup1,text='Clear',command=ClearCanvas)
                    self.popup_btn1.grid(row=1,column=1)
                    self.popup_btn2 = Button(self.frame_popup1,text='RePlot',command=RePlot)
                    self.popup_btn2.grid(row=2,column=1)
                
                def ClearCanvas():
                    self.cv.destroy()
                    self.cv1.destroy()
                    self.cv2.destroy()
                    self.cv3.destroy()
                    self.cv = Canvas(self.frame_popup)
                    self.cv.grid(row=0,column=0)
                    self.cv1 = Canvas(self.frame_popup)
                    self.cv1.grid(row=0,column=1)
                    self.cv2 = Canvas(self.frame_popup)
                    self.cv2.grid(row=1,column=0)
                    self.cv3 = Canvas(self.frame_popup)
                    self.cv3.grid(row=1,column=1)
                    
                    
                def RePlot():
                    ClearCanvas()
                    PlotEnergies(Ecf_val,Ecf_val_excitation,self.cv,r'Crystal Field')
                    PlotEnergies(Ecf_so_val,Ecf_so_val_excitation,self.cv1,r'Spin-Orbit coupling added')
                    PlotEnergies(Ecf_so_ss_val,Ecf_so_ss_val_excitation,self.cv2,r'Spin-Spin correlation added')
                    PlotEnergies(Ecf_so_ss_g_val,Ecf_so_ss_g_val_excitation,self.cv3,r'Molecular field added')
                    
                def SavePlot():
                    ClearCanvas()
                    current_directory = askdirectory()
                    file_name = "test"
                    file_path = os.path.join(current_directory,file_name)
                    print(file_path)
                    PlotEnergies(Ecf_val,Ecf_val_excitation,self.cv,r'Crystal Field')
                    fig.savefig(file_path+'1.png',dpi=600)
                    PlotEnergies(Ecf_so_val,Ecf_so_val_excitation,self.cv1,r'Spin-Orbit coupling added')
                    fig.savefig(file_path+'2.png',dpi=600)
                    PlotEnergies(Ecf_so_ss_val,Ecf_so_ss_val_excitation,self.cv2,r'Spin-Spin correlation added')
                    fig.savefig(file_path+'3.png',dpi=600)
                    PlotEnergies(Ecf_so_ss_g_val,Ecf_so_ss_g_val_excitation,self.cv3,r'Molecular field added')
                    fig.savefig(file_path+'4.png',dpi=600)
                    

                    
                    # asksaveasfilename, askopenfilename, askdirectory
                    
                    
                PopUp()

                PlotEnergies(Ecf_val,Ecf_val_excitation,self.cv,r'Crystal Field')
                PlotEnergies(Ecf_so_val,Ecf_so_val_excitation,self.cv1,r'Spin-Orbit coupling added')
                PlotEnergies(Ecf_so_ss_val,Ecf_so_ss_val_excitation,self.cv2,r'Spin-Spin correlation added')
                PlotEnergies(Ecf_so_ss_g_val,Ecf_so_ss_g_val_excitation,self.cv3,r'Molecular field added')
                
                
                self.progress['value'] = 50
                win.update()
                
                
                self.results1.insert(tk.END,"----------- Initializing ---------\n")
                self.results1.insert(tk.END, 'L: ' + self.L.get() + '\n')
                self.results1.insert(tk.END, 'S: ' + self.S.get() + '\n')
                self.results1.insert(tk.END, 'ligand numbers: ' + self.num_lig.get() + '\n')
                self.results1.insert(tk.END, 'ligand positions:\n')
                self.results1.insert(tk.END, self.lig_pos.get('1.0',tk.END))
                self.results1.insert(tk.END, typeB(B))
                self.results1.insert(tk.END,'SO value (meV): ' + str(SO_val) + '\n')
                self.results1.insert(tk.END,'SS value (meV): ' + str(SS_val) + '\n')
                self.results1.insert(tk.END,'molecular field (meV): ' + str(g_val) + '\n')
                self.progress['value'] = 100
                win.update()
                
            except:
                self.results1.insert(tk.END,"-------------------------------\n")
                self.results1.insert(tk.END,'something is wrong\n')
                print('something is wrong')
                self.progress['value'] = 100
                win.update()
                
                
        ####################################################################################
        
        
        frame2 = Frame(frame)
        frame2.grid(row=1,column=0,padx=5)
        
                
        self.btn = Button(frame2,text='initialize',command=initialize)
        #self.btn.place(x=10,y=310)
        self.btn.grid(row=0,column=0,pady=5)
                    
        self.progress = Progressbar(frame2, orient = tk.HORIZONTAL, maximum = 100, length = 100, mode = 'determinate')
        #self.progress.place(x=70,y=313,width=200)
        self.progress.grid(row=0,column=1,padx=30)
        
               
        
        ####################################################################################
        frame3 = Frame(win)
        frame3.grid(row=1,column=0,padx=5)
        
        
        self.results = Label(frame3,text='results:')
        #self.results.place(x=10,y=550)
        self.results.grid(row=0,column=0,pady=5)
        
        self.results1 = scrolledtext.ScrolledText(frame3,width=80,height=15)
        self.results1.place(x=10,y=570)
        self.results1.grid(row=1,column=0,pady=5)
        
        def clear_results():
            self.results1.delete('1.0', tk.END)
            
        self.btn = Button(frame3,text='clear',command=clear_results)
        #self.btn.place(x=10,y=660)
        self.btn.grid(row=1,column=1,padx=5)

        ####################################################################################

        
        

#styles = ThemedStyle().tk.call('ttk::themes')
#t = styles[12] # 'arc'
#t = styles[17] # 'yaru'
#t = styles[14] # 'radiance'
#t = styles[9] # 'equilux'
t = 'arc'
start=StartUp(t)


window = tk.Tk()

style = ThemedStyle(window)
'''
def styleChange():
    style.theme_use(t)
    color = Style().lookup("TFrame", "background", default="black")
    window.configure(bg=color)
'''
style.theme_use(t)
color = Style().lookup("TFrame", "background", default="black")
window.configure(bg=color)


mywin=MyWindow(window)
window.title('Crystal Field calculation')


#screen_width = window.winfo_screenwidth()
#screen_height = window.winfo_screenheight()
#window.geometry('%dx%d' % (screen_width-20, screen_height-100))
window.eval('tk::PlaceWindow . center')

#window.attributes('-disabled', True)
#window.resizable(0,0)


window.mainloop()



