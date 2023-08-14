# -*- coding: utf-8 -*-
"""
Created on Wed Dec 29 17:18:53 2021

@author: brian
"""
import Crystal_Field_Calcuations as cef
from PyQt5.QtWidgets import QApplication
from PyQt5.QtWidgets import QMainWindow
from PyQt5.QtWidgets import QLineEdit
from PyQt5.QtWidgets import QLabel
from PyQt5.QtWidgets import QPushButton
from PyQt5.QtWidgets import QMessageBox
from PyQt5.QtWidgets import QTableWidget
from PyQt5.QtWidgets import QTableWidgetItem
from PyQt5.QtWidgets import QHeaderView
from PyQt5.QtWidgets import QVBoxLayout


from PyQt5.QtCore import Qt
from PyQt5.QtCore import pyqtSlot
from PyQt5.QtGui import QPalette
from PyQt5 import QtCore
import sys
import numpy as np






class App(QMainWindow):

    def __init__(self):
        super().__init__()
        self.title = 'Crystal Field Calculation'
        self.left = 50
        self.top = 50
        self.width = 1800
        self.height = 900
        self.initUI()
    
    def initUI(self):
        self.setWindowTitle(self.title)
        self.setGeometry(self.left, self.top, self.width, self.height)
        
        
        # ion label
        self.ion_label = QLabel('Ion:',self)
        self.ion_label.move(20,20)
        self.ion_label.resize(90,40)
        
        
        # ion textbox
        self.ion = QLineEdit(self)
        self.ion.move(80, 20)
        self.ion.resize(40,40)
        
        
        # L label
        self.L_label = QLabel('L:',self)
        self.L_label.move(20,60)
        self.L_label.resize(80,40)
        
        
        # L textbox
        self.L = QLineEdit(self)
        self.L.move(80, 60)
        self.L.resize(40,40)
        
        # S label
        self.S_label = QLabel('S:',self)
        self.S_label.move(20,100)
        self.S_label.resize(80,40)
        
        
        # S textbox
        self.S = QLineEdit(self)
        self.S.move(80, 100)
        self.S.resize(40,40)
        
        
        # Z label
        self.Z_label = QLabel('Z:',self)
        self.Z_label.move(20,140)
        self.Z_label.resize(80,40)
        
        
        # Z textbox
        self.Z = QLineEdit(self)
        self.Z.move(80, 140)
        self.Z.resize(40,40)
        
        
        # SO label
        self.SO_label = QLabel('Spin-Orbit (meV):',self)
        self.SO_label .move(20,180)
        self.SO_label .resize(160,40)
        
        
        # SO textbox
        self.SO = QLineEdit(self)
        self.SO.move(160, 180)
        self.SO.resize(40,40)
        
        """
        # Ion_position label
        self.Ion_position = QLabel('Ion_position (x,y,z):',self)
        self.Ion_position .move(20,240)
        self.Ion_position .resize(110,40)
        
        
        # Ion_position textbox for x
        self.Ion_position_x = QLineEdit(self)
        self.Ion_position_x.move(140, 240)
        self.Ion_position_x.resize(40,40)
        # Ion_position textbox for y
        self.Ion_position_y = QLineEdit(self)
        self.Ion_position_y.move(180, 240)
        self.Ion_position_y.resize(40,40)
        # Ion_position textbox for z
        self.Ion_position_z = QLineEdit(self)
        self.Ion_position_z.move(220, 240)
        self.Ion_position_z.resize(40,40)
        
        
        # Ligand 1 label
        self.Li1 = QLabel('Ligand 1 (x,y,z):',self)
        self.Li1 .move(20,280)
        self.Li1 .resize(110,40)
        
        
        # Ligand 1 textbox for x
        self.Li1_x = QLineEdit(self)
        self.Li1_x.move(140, 280)
        self.Li1_x.resize(40,40)
        # Ligand 1 textbox for y
        self.Li1_y = QLineEdit(self)
        self.Li1_y.move(180, 280)
        self.Li1_y.resize(40,40)
        # Ligand 1 textbox for z
        self.Li1_z = QLineEdit(self)
        self.Li1_z.move(220, 280)
        self.Li1_z.resize(40,40)
        
        
        # Ligand 2 label
        self.Li2 = QLabel('Ligand 2 (x,y,z):',self)
        self.Li2 .move(20,320)
        self.Li2 .resize(110,40)
        
        
        # Ligand 2 textbox for x
        self.Li2_x = QLineEdit(self)
        self.Li2_x.move(140, 320)
        self.Li2_x.resize(40,40)
        # Ligand 2 textbox for y
        self.Li2_y = QLineEdit(self)
        self.Li2_y.move(180, 320)
        self.Li2_y.resize(40,40)
        # Ligand 2 textbox for z
        self.Li2_z = QLineEdit(self)
        self.Li2_z.move(220, 320)
        self.Li2_z.resize(40,40)
        
        
        # Ligand 3 label
        self.Li3 = QLabel('Ligand 3 (x,y,z):',self)
        self.Li3 .move(20,360)
        self.Li3 .resize(110,40)
        
        
        # Ligand 3 textbox for x
        self.Li3_x = QLineEdit(self)
        self.Li3_x.move(140, 360)
        self.Li3_x.resize(40,40)
        # Ligand 3 textbox for y
        self.Li3_y = QLineEdit(self)
        self.Li3_y.move(180, 360)
        self.Li3_y.resize(40,40)
        # Ligand 3 textbox for z
        self.Li3_z = QLineEdit(self)
        self.Li3_z.move(220, 360)
        self.Li3_z.resize(40,40)
        
        
        
        # Ligand 4 label
        self.Li4 = QLabel('Ligand 4 (x,y,z):',self)
        self.Li4 .move(20,400)
        self.Li4 .resize(110,40)
        
        
        # Ligand 4 textbox for x
        self.Li4_x = QLineEdit(self)
        self.Li4_x.move(140, 400)
        self.Li4_x.resize(40,40)
        # Ligand 4 textbox for y
        self.Li4_y = QLineEdit(self)
        self.Li4_y.move(180, 400)
        self.Li4_y.resize(40,40)
        # Ligand 4 textbox for z
        self.Li4_z = QLineEdit(self)
        self.Li4_z.move(220, 400)
        self.Li4_z.resize(40,40)
        
        
        # Ligand 5 label
        self.Li5 = QLabel('Ligand 5 (x,y,z):',self)
        self.Li5 .move(20,440)
        self.Li5 .resize(110,40)
        
        
        # Ligand 5 textbox for x
        self.Li5_x = QLineEdit(self)
        self.Li5_x.move(140, 440)
        self.Li5_x.resize(40,40)
        # Ligand 5 textbox for y
        self.Li5_y = QLineEdit(self)
        self.Li5_y.move(180, 440)
        self.Li5_y.resize(40,40)
        # Ligand 5 textbox for z
        self.Li5_z = QLineEdit(self)
        self.Li5_z.move(220, 440)
        self.Li5_z.resize(40,40)
        
        
        # Ligand 6 label
        self.Li6 = QLabel('Ligand 6 (x,y,z):',self)
        self.Li6 .move(20,480)
        self.Li6 .resize(110,40)
        
        
        # Ligand 6 textbox for x
        self.Li6_x = QLineEdit(self)
        self.Li6_x.move(140, 480)
        self.Li6_x.resize(40,40)
        # Ligand 6 textbox for y
        self.Li6_y = QLineEdit(self)
        self.Li6_y.move(180, 480)
        self.Li6_y.resize(40,40)
        # Ligand 6 textbox for z
        self.Li6_z = QLineEdit(self)
        self.Li6_z.move(220, 480)
        self.Li6_z.resize(40,40)
        
        
        
        # Ligand 7 label
        self.Li7 = QLabel('Ligand 7 (x,y,z):',self)
        self.Li7 .move(20,520)
        self.Li7 .resize(110,40)
        
        
        # Ligand 7 textbox for x
        self.Li7_x = QLineEdit(self)
        self.Li7_x.move(140, 520)
        self.Li7_x.resize(40,40)
        # Ligand 7 textbox for y
        self.Li7_y = QLineEdit(self)
        self.Li7_y.move(180, 520)
        self.Li7_y.resize(40,40)
        # Ligand 7 textbox for z
        self.Li7_z = QLineEdit(self)
        self.Li7_z.move(220, 520)
        self.Li7_z.resize(40,40)
        
        
        # Ligand 8 label
        self.Li8 = QLabel('Ligand 8 (x,y,z):',self)
        self.Li8 .move(20,560)
        self.Li8 .resize(110,40)
        
        
        # Ligand 8 textbox for x
        self.Li8_x = QLineEdit(self)
        self.Li8_x.move(140, 560)
        self.Li8_x.resize(40,40)
        # Ligand 8 textbox for y
        self.Li8_y = QLineEdit(self)
        self.Li8_y.move(180, 560)
        self.Li8_y.resize(40,40)
        # Ligand 8 textbox for z
        self.Li8_z = QLineEdit(self)
        self.Li8_z.move(220, 560)
        self.Li8_z.resize(40,40)
        

        """
        
        # Positions of ions and ligands
        self.Ligands = QLabel('Ion_position (x,y,z):',self)
        self.Ligands.move(20,220)
        self.Ligands.resize(110,40)
        
        
        # Ion_position textbox for x
        self.Ligands_t = QLineEdit(self)
        self.Ligands_t.move(20, 260)
        self.Ligands_t.resize(240,340)
        
        
        
        # B parameters label
        self.B_parameters_label = QLabel('B parameters',self)
        self.B_parameters_label .move(300,20)
        self.B_parameters_label .resize(200,40)
        # Create a table for B parameters
        self.B_parameters = QTableWidget(self)
        self.B_parameters.move(300 , 50)
        self.B_parameters.resize(300,400)
        self.B_parameters.setRowCount(8)
        self.B_parameters.setColumnCount(2)
        self.B_parameters.setItem(0,0, QTableWidgetItem("B20 (meV)"))
        self.B_parameters.setItem(1,0, QTableWidgetItem("B21 (meV)"))
        self.B_parameters.setItem(2,0, QTableWidgetItem("B22 (meV)"))
        self.B_parameters.setItem(3,0, QTableWidgetItem("B40 (meV)"))
        self.B_parameters.setItem(4,0, QTableWidgetItem("B41 (meV)"))
        self.B_parameters.setItem(5,0, QTableWidgetItem("B42 (meV)"))
        self.B_parameters.setItem(6,0, QTableWidgetItem("B43 (meV)"))
        self.B_parameters.setItem(7,0, QTableWidgetItem("B44 (meV)"))
        
        self.B_parameters.horizontalHeader().setStretchLastSection(True)
        self.B_parameters.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)
        self.B_parameters.verticalHeader().setStretchLastSection(True)
        self.B_parameters.verticalHeader().setSectionResizeMode(QHeaderView.Stretch)




        
        
        # Energies label
        self.Energies_label = QLabel('Energies (meV)',self)
        self.Energies_label .move(650,20)
        self.Energies_label .resize(100,40)
        # Create a table for Energies
        self.Energies = QTableWidget(self)
        self.Energies.move(650 , 50)
        self.Energies.resize(200,400)
        self.Energies.setRowCount(1)
        self.Energies.setColumnCount(1)
        self.Energies.setUpdatesEnabled(True)
        
        self.Energies.horizontalHeader().setStretchLastSection(True)
        self.Energies.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)
        self.Energies.verticalHeader().setStretchLastSection(True)
        self.Energies.verticalHeader().setSectionResizeMode(QHeaderView.Stretch)

        
        
        # Energies + SO label
        self.Energies_SO_label = QLabel('Energies + SO (meV)',self)
        self.Energies_SO_label .move(900,20)
        self.Energies_SO_label .resize(150,40)
        # Create a table for Energies
        self.Energies_SO = QTableWidget(self)
        self.Energies_SO.move(900 , 50)
        self.Energies_SO.resize(200,400)
        self.Energies_SO.setRowCount(1)
        self.Energies_SO.setColumnCount(1)
        self.Energies_SO.setUpdatesEnabled(True)

        
        self.Energies_SO.horizontalHeader().setStretchLastSection(True)
        self.Energies_SO.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)
        self.Energies_SO.verticalHeader().setStretchLastSection(True)
        self.Energies_SO.verticalHeader().setSectionResizeMode(QHeaderView.Stretch)
        
        
        #self.graphWidget = pg.PlotWidget()
        #self.setCentralWidget(self.graphWidget)
        
        
        
        # Just some button to start the calculation
        self.button = QPushButton('Calculate', self)
        # adding action to the button
        self.button.clicked.connect(self.on_click)
        self.button.move(120,620)
        self.button.resize(150,50)
        
        
        
        self.show()
        

    @pyqtSlot()
    def on_click(self):
        ion_text = self.ion.text()
        L_text = self.L.text()
        S_text = self.S.text()
        Z_text = self.Z.text()
        SO_text = self.SO.text()
        QMessageBox.question(self, 'Crystal Field', "You typed: ion = " + ion_text +\
        ", L = " + L_text + ", S = " + S_text + ", Z = " + Z_text + ", SO = " +\
        SO_text, QMessageBox.Ok, QMessageBox.Ok)
        
        ion = ion_text
        L = int(L_text)
        S = int(S_text)
        Z = int(Z_text)
        calc = cef.LS(L,S)
        Degeneracy = (2*L+1)*(2*S+1)
        
        O20 = calc.Olm(L,S,2,0)
        O21 = calc.Olm(L,S,2,1)
        O22 = calc.Olm(L,S,2,2)
        O40 = calc.Olm(L,S,4,0)
        O41 = calc.Olm(L,S,4,1)
        O42 = calc.Olm(L,S,4,2)
        O43 = calc.Olm(L,S,4,3)
        O44 = calc.Olm(L,S,4,4)
        
        d = np.zeros((9,3))
        positions = [[0,0,0]]
        
        Ligands_positions = np.array([self.Ligands_t.split()])
        """
        d[0] = np.array([float(self.Ligands_positions[1]),    float(self.Ligands_positions[2]),    float(self.Ligands_positions[3])])
        d[1] = np.array([float(self.Ligands_positions[5]),    float(self.Ligands_positions[6]),    float(self.Ligands_positions[7])])
        d[2] = np.array([float(self.Ligands_positions[9]),    float(self.Ligands_positions[10]),    float(self.Ligands_positions[11])])
        d[3] = np.array([float(self.Ligands_positions[13]),    float(self.Ligands_positions[14]),    float(self.Ligands_positions[15])])
        d[4] = np.array([float(self.Ligands_positions[17]),    float(self.Ligands_positions[2]),    float(self.Ligands_positions[3])])
        d[5] = np.array([float(self.Ligands_positions[21]),    float(self.Ligands_positions[2]),    float(self.Ligands_positions[3])])
        d[6] = np.array([float(self.Ligands_positions[25]),    float(self.Ligands_positions[2]),    float(self.Ligands_positions[3])])
        d[7] = np.array([float(self.Ligands_positions[1]),    float(self.Ligands_positions[2]),    float(self.Ligands_positions[3])])
        d[8] = np.array([float(self.Ligands_positions[1]),    float(self.Ligands_positions[2]),    float(self.Ligands_positions[3])])
        for k in range(8):
            if d[k+1][0] == 0 and d[k+1][1] == 0 and d[k+1][2] == 0:
                pass
            else:
                d[k+1] = d[k+1] - d[0]
                positions = np.vstack((positions,np.array([d[k+1]])))
        positions = np.delete(positions,0,0)
        """
        
        B = calc.PC(ion,L,S,positions,Z)
        
        
        QMessageBox.question(self, 'Crystal Field', 'Results are down below' +\
        '\n B20: ' + str(B[0]) +\
        '\n B21: ' + str(B[1]) +\
        '\n B22: ' + str(B[2]) +\
        '\n B40: ' + str(B[3]) +\
        '\n B41: ' + str(B[4]) +\
        '\n B42: ' + str(B[5]) +\
        '\n B43: ' + str(B[6]) +\
        '\n B44: ' + str(B[7]), QMessageBox.Ok, QMessageBox.Ok)

        self.B_parameters.setItem(0,1, QTableWidgetItem(str(B[0])))
        self.B_parameters.setItem(1,1, QTableWidgetItem(str(B[1])))
        self.B_parameters.setItem(2,1, QTableWidgetItem(str(B[2])))
        self.B_parameters.setItem(3,1, QTableWidgetItem(str(B[3])))
        self.B_parameters.setItem(4,1, QTableWidgetItem(str(B[4])))
        self.B_parameters.setItem(5,1, QTableWidgetItem(str(B[5])))
        self.B_parameters.setItem(6,1, QTableWidgetItem(str(B[6])))
        self.B_parameters.setItem(7,1, QTableWidgetItem(str(B[7])))

        
        Hcf = B[0]*O20 + B[1]*O21 + B[2]*O22 + B[3]*O40 + B[4]*O41 + B[5]*O42 + B[6]*O43 + B[7]*O44
        self.Ecf_val,self.Ecf_val_excitation,self.H_cf_vt = calc.Diag(Hcf,printfunction=True)
        if SO_text == '' or int(SO_text) == 0:
            SO_matrix = calc.SO(ion,L,S)
        else:
            SO_matrix = calc.SO(ion,L,S,float(SO_text))
        Hcf_so = Hcf + SO_matrix
        self.Ecf_so_val,self.Ecf_so_val_excitation,H_cf_so_vt = calc.Diag(Hcf_so,printfunction=True)
        
        
        for k in range(Degeneracy):
            self.Energies.setItem(k,0, QTableWidgetItem(str(self.Ecf_val_excitation[k])))
            self.Energies.insertRow(self.Energies.rowCount())

            self.Energies_SO.setItem(k,0, QTableWidgetItem(str(self.Ecf_so_val_excitation[k])))
            self.Energies_SO.insertRow(self.Energies_SO.rowCount())
        
        self.Energies.removeRow(Degeneracy)
        self.Energies_SO.removeRow(Degeneracy)
        
        
        
        self.fig = QApplication([])
        
        
        
        
            
        ax1 = self.fig.add_subplot(221)
        ax1.eventplot(self.Ecf_val_excitation,orientation='horizontal',linelength=0.05,linestyles='solid',colors='r')
        ax1.tick_params(axis='x',which='both',bottom=False,top=False)
        ax1.set_xlabel('Crystal Field Splitting')
        ax1.set_ylabel('Energy (meV)')
        ax1.set_xticklabels('')
        ax1.set_title('Point Charge Model')
        ax1.set_xlim([0.95,1.05])
        
        ax2 = self.fig.add_subplot(222)
        ax2.eventplot(self.Ecf_so_val_excitation,orientation='horizontal',linelength=0.05,linestyles='solid',colors='b')
        ax2.tick_params(axis='y',which='both',bottom=False,top=False)
        ax2.set_xlabel('Energy (meV)')
        ax2.set_ylabel('Arb. Units')
        ax2.set_yticklabels('')
        ax2.set_yticks([])
        ax2.set_title('Crystal Field excitations')
        ax2.set_ylim([0.95,1.05])
        









        
if __name__ == '__main__':
    app = QApplication(sys.argv)
    
    # setting the colors
    palette = QPalette()
    palette.setColor(QPalette.Window, Qt.lightGray)
    palette.setColor(QPalette.ButtonText, Qt.red)
    palette.setColor(QPalette.WindowText, Qt.darkBlue)
    app.setPalette(palette)
    ex = App()

    sys.exit(app.exec_())


