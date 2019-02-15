#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 28 14:04:17 2018

@author: boss
"""
import numpy as np
import matplotlib
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2TkAgg
from matplotlib.figure import Figure
from matplotlib import style

import tkinter as tk
from tkinter import ttk
from tkinter import messagebox

LARGE_FONT= ("Verdana", 12)  #type of font and size
EXTRALG_FONT= ("Verdana", 18) 
EXTRALG_FONT2= ("Verdana", 14) 
style.use("ggplot") #ggplot, dark_background,

#species
S=['H2', 'H', 'O', 'O2', 'OH', 'H2O', 'HO2', 'H2O2', 'C', 'CH', 'CH2', 'CH2(S)',
   'CH3', 'CH4', 'CO', 'CO2', 'HCO', 'CH2O', 'CH2OH', 'CH3O', 'CH3OH', 'C2H',
   'C2H2', 'C2H3', 'C2H4', 'C2H5', 'C2H6', 'HCCO', 'CH2CO', 'HCCOH', 'N', 'NH',
   'NH2', 'NH3', 'NNH', 'NO', 'NO2', 'N2O', 'HNO', 'CN', 'HCN', 'H2CN', 'HCNN',
   'HCNO', 'HOCN', 'HNCO', 'NCO', 'N2', 'AR', 'C3H7', 'C3H8', 'CH2CHO', 'CH3CHO']
#specificspecies
SS = SpecificSpecies
#reactions
R=(np.linspace(1, 325,325))
RR = ['2 O + M <=> O2 + M','O + H + M <=> OH + M','O + H2 <=> H + OH','O + HO2 <=> OH + O2',
      'O + H2O2 <=> OH + HO2','O + CH <=> H + CO','O + CH2 <=> H + HCO','O + CH2(S) <=> H2 + CO',
      'O + CH2(S) <=> H + HCO','O + CH3 <=> H + CH2O','O + CH4 <=> OH + CH3','O + CO (+ M) <=> CO2 (+ M)',
      'O + HCO <=> OH + CO','O + HCO <=> H + CO2','O + CH2O <=> OH + HCO','O + CH2OH <=> OH + CH2O',
      'O + CH3O <=> OH + CH2O','O + CH3OH <=> OH + CH2OH','O + CH3OH <=> OH + CH3O','O + C2H <=> CH + CO',
      'O + C2H2 <=> H + HCCO','O + C2H2 <=> OH + C2H','O + C2H2 <=> CO + CH2','O + C2H3 <=> H + CH2CO',
      'O + C2H4 <=> CH3 + HCO','O + C2H5 <=> CH3 + CH2O','O + C2H6 <=> OH + C2H5','O + HCCO <=> H + 2 CO',
      'O + CH2CO <=> OH + HCCO','O + CH2CO <=> CH2 + CO2','O2 + CO <=> O + CO2','O2 + CH2O <=> HO2 + HCO',
      'H + O2 + M <=> HO2 + M','H + 2 O2 <=> HO2 + O2','H + O2 + H2O <=> HO2 + H2O','H + O2 + N2 <=> HO2 + N2',
      'H + O2 + AR <=> HO2 + AR','H + O2 <=> O + OH','2 H + M <=> H2 + M','2 H + H2 <=> 2 H2',
      '2 H + H2O <=> H2 + H2O','2 H + CO2 <=> H2 + CO2','H + OH + M <=> H2O + M','H + HO2 <=> O + H2O',
      'H + HO2 <=> O2 + H2','H + HO2 <=> 2 OH','H + H2O2 <=> HO2 + H2','H + H2O2 <=> OH + H2O',
      'H + CH <=> C + H2','H + CH2 (+ M) <=> CH3 (+ M)','H + CH2(S) <=> CH + H2','H + CH3 (+ M) <=> CH4 (+ M)',
      'H + CH4 <=> CH3 + H2','H + HCO (+ M) <=> CH2O (+ M)','H + HCO <=> H2 + CO','H + CH2O (+ M) <=> CH2OH (+ M)',
      'H + CH2O (+ M) <=> CH3O (+ M)','H + CH2O <=> HCO + H2','H + CH2OH (+ M) <=> CH3OH (+ M)','H + CH2OH <=> H2 + CH2O',
      'H + CH2OH <=> OH + CH3','H + CH2OH <=> CH2(S) + H2O','H + CH3O (+ M) <=> CH3OH (+ M)','H + CH3O <=> H + CH2OH',
      'H + CH3O <=> H2 + CH2O','H + CH3O <=> OH + CH3','H + CH3O <=> CH2(S) + H2O','H + CH3OH <=> CH2OH + H2',
      'H + CH3OH <=> CH3O + H2','H + C2H (+ M) <=> C2H2 (+ M)','H + C2H2 (+ M) <=> C2H3 (+ M)','H + C2H3 (+ M) <=> C2H4 (+ M)',
      'H + C2H3 <=> H2 + C2H2','H + C2H4 (+ M) <=> C2H5 (+ M)','H + C2H4 <=> C2H3 + H2','H + C2H5 (+ M) <=> C2H6 (+ M)',
      'H + C2H5 <=> H2 + C2H4','H + C2H6 <=> C2H5 + H2','H + HCCO <=> CH2(S) + CO','H + CH2CO <=> HCCO + H2',
      'H + CH2CO <=> CH3 + CO','H + HCCOH <=> H + CH2CO','H2 + CO (+ M) <=> CH2O (+ M)','OH + H2 <=> H + H2O',
      '2 OH (+ M) <=> H2O2 (+ M)','2 OH <=> O + H2O','OH + HO2 <=> O2 + H2O','OH + H2O2 <=> HO2 + H2O',
      'OH + H2O2 <=> HO2 + H2O','OH + C <=> H + CO','OH + CH <=> H + HCO','OH + CH2 <=> H + CH2O',
      'OH + CH2 <=> CH + H2O','OH + CH2(S) <=> H + CH2O','OH + CH3 (+ M) <=> CH3OH (+ M)','OH + CH3 <=> CH2 + H2O',
      'OH + CH3 <=> CH2(S) + H2O','OH + CH4 <=> CH3 + H2O','OH + CO <=> H + CO2','OH + HCO <=> H2O + CO',
      'OH + CH2O <=> HCO + H2O','OH + CH2OH <=> H2O + CH2O','OH + CH3O <=> H2O + CH2O','OH + CH3OH <=> CH2OH + H2O',
      'OH + CH3OH <=> CH3O + H2O','OH + C2H <=> H + HCCO','OH + C2H2 <=> H + CH2CO','OH + C2H2 <=> H + HCCOH',
      'OH + C2H2 <=> C2H + H2O','OH + C2H2 <=> CH3 + CO','OH + C2H3 <=> H2O + C2H2','OH + C2H4 <=> C2H3 + H2O',
      'OH + C2H6 <=> C2H5 + H2O','OH + CH2CO <=> HCCO + H2O','2 HO2 <=> O2 + H2O2','2 HO2 <=> O2 + H2O2',
      'HO2 + CH2 <=> OH + CH2O','HO2 + CH3 <=> O2 + CH4','HO2 + CH3 <=> OH + CH3O','HO2 + CO <=> OH + CO2',
      'HO2 + CH2O <=> HCO + H2O2','C + O2 <=> O + CO','C + CH2 <=> H + C2H','C + CH3 <=> H + C2H2',
      'CH + O2 <=> O + HCO','CH + H2 <=> H + CH2','CH + H2O <=> H + CH2O','CH + CH2 <=> H + C2H2',
      'CH + CH3 <=> H + C2H3','CH + CH4 <=> H + C2H4','CH + CO (+ M) <=> HCCO (+ M)','CH + CO2 <=> HCO + CO',
      'CH + CH2O <=> H + CH2CO','CH + HCCO <=> CO + C2H2','CH2 + O2 => OH + H + CO','CH2 + H2 <=> H + CH3',
      '2 CH2 <=> H2 + C2H2','CH2 + CH3 <=> H + C2H4','CH2 + CH4 <=> 2 CH3','CH2 + CO (+ M) <=> CH2CO (+ M)',
     'CH2 + HCCO <=> C2H3 + CO','CH2(S) + N2 <=> CH2 + N2','CH2(S) + AR <=> CH2 + AR','CH2(S) + O2 <=> H + OH + CO',
     'CH2(S) + O2 <=> CO + H2O','CH2(S) + H2 <=> CH3 + H','CH2(S) + H2O (+ M) <=> CH3OH (+ M)','CH2(S) + H2O <=> CH2 + H2O',
     'CH2(S) + CH3 <=> H + C2H4','CH2(S) + CH4 <=> 2 CH3','CH2(S) + CO <=> CH2 + CO','CH2(S) + CO2 <=> CH2 + CO2',
     'CH2(S) + CO2 <=> CO + CH2O','CH2(S) + C2H6 <=> CH3 + C2H5','CH3 + O2 <=> O + CH3O','CH3 + O2 <=> OH + CH2O',
     'CH3 + H2O2 <=> HO2 + CH4','2 CH3 (+ M) <=> C2H6 (+ M)','2 CH3 <=> H + C2H5','CH3 + HCO <=> CH4 + CO',
     'CH3 + CH2O <=> HCO + CH4','CH3 + CH3OH <=> CH2OH + CH4','CH3 + CH3OH <=> CH3O + CH4','CH3 + C2H4 <=> C2H3 + CH4',
     'CH3 + C2H6 <=> C2H5 + CH4','HCO + H2O <=> H + CO + H2O','HCO + M <=> H + CO + M','HCO + O2 <=> HO2 + CO',
     'CH2OH + O2 <=> HO2 + CH2O','CH3O + O2 <=> HO2 + CH2O','C2H + O2 <=> HCO + CO','C2H + H2 <=> H + C2H2',
     'C2H3 + O2 <=> HCO + CH2O','C2H4 (+ M) <=> H2 + C2H2 (+ M)','C2H5 + O2 <=> HO2 + C2H4','HCCO + O2 <=> OH + 2 CO',
     '2 HCCO <=> 2 CO + C2H2','N + NO <=> N2 + O','N + O2 <=> NO + O','N + OH <=> NO + H',
     'N2O + O <=> N2 + O2','N2O + O <=> 2 NO','N2O + H <=> N2 + OH','N2O + OH <=> N2 + HO2',
     'N2O (+ M) <=> N2 + O (+ M)','HO2 + NO <=> NO2 + OH','NO + O + M <=> NO2 + M','NO2 + O <=> NO + O2',
     'NO2 + H <=> NO + OH','NH + O <=> NO + H','NH + H <=> N + H2','NH + OH <=> HNO + H',
     'NH + OH <=> N + H2O','NH + O2 <=> HNO + O','NH + O2 <=> NO + OH','NH + N <=> N2 + H',
     'NH + H2O <=> HNO + H2','NH + NO <=> N2 + OH','NH + NO <=> N2O + H','NH2 + O <=> OH + NH',
     'NH2 + O <=> H + HNO','NH2 + H <=> NH + H2','NH2 + OH <=> NH + H2O','NNH <=> N2 + H',
     'NNH + M <=> N2 + H + M','NNH + O2 <=> HO2 + N2','NNH + O <=> OH + N2','NNH + O <=> NH + NO',
     'NNH + H <=> H2 + N2','NNH + OH <=> H2O + N2','NNH + CH3 <=> CH4 + N2','H + NO + M <=> HNO + M',
     'HNO + O <=> NO + OH','HNO + H <=> H2 + NO','HNO + OH <=> NO + H2O','HNO + O2 <=> HO2 + NO',
     'CN + O <=> CO + N','CN + OH <=> NCO + H','CN + H2O <=> HCN + OH','CN + O2 <=> NCO + O',
     'CN + H2 <=> HCN + H','NCO + O <=> NO + CO','NCO + H <=> NH + CO','NCO + OH <=> NO + H + CO',
     'NCO + N <=> N2 + CO','NCO + O2 <=> NO + CO2','NCO + M <=> N + CO + M','NCO + NO <=> N2O + CO',
     'NCO + NO <=> N2 + CO2','HCN + M <=> H + CN + M','HCN + O <=> NCO + H','HCN + O <=> NH + CO',
     'HCN + O <=> CN + OH','HCN + OH <=> HOCN + H','HCN + OH <=> HNCO + H','HCN + OH <=> NH2 + CO',
     'H + HCN (+ M) <=> H2CN (+ M)','H2CN + N <=> N2 + CH2','C + N2 <=> CN + N','CH + N2 <=> HCN + N',
     'CH + N2 (+ M) <=> HCNN (+ M)','CH2 + N2 <=> HCN + NH','CH2(S) + N2 <=> NH + HCN','C + NO <=> CN + O',
     'C + NO <=> CO + N','CH + NO <=> HCN + O','CH + NO <=> H + NCO','CH + NO <=> N + HCO',
     'CH2 + NO <=> H + HNCO','CH2 + NO <=> OH + HCN','CH2 + NO <=> H + HCNO','CH2(S) + NO <=> H + HNCO',
     'CH2(S) + NO <=> OH + HCN','CH2(S) + NO <=> H + HCNO','CH3 + NO <=> HCN + H2O','CH3 + NO <=> H2CN + OH',
     'HCNN + O <=> CO + H + N2','HCNN + O <=> HCN + NO','HCNN + O2 <=> O + HCO + N2','HCNN + OH <=> H + HCO + N2',
     'HCNN + H <=> CH2 + N2','HNCO + O <=> NH + CO2','HNCO + O <=> HNO + CO','HNCO + O <=> NCO + OH',
     'HNCO + H <=> NH2 + CO','HNCO + H <=> H2 + NCO','HNCO + OH <=> NCO + H2O','HNCO + OH <=> NH2 + CO2',
     'HNCO + M <=> NH + CO + M','HCNO + H <=> H + HNCO','HCNO + H <=> OH + HCN','HCNO + H <=> NH2 + CO',
     'HOCN + H <=> H + HNCO','HCCO + NO <=> HCNO + CO','CH3 + N <=> H2CN + H','CH3 + N <=> HCN + H2',
     'NH3 + H <=> NH2 + H2','NH3 + OH <=> NH2 + H2O','NH3 + O <=> NH2 + OH','NH + CO2 <=> HNO + CO',
     'CN + NO2 <=> NCO + NO','NCO + NO2 <=> N2O + CO2','N + CO2 <=> NO + CO','O + CH3 => H + H2 + CO',
     'O + C2H4 <=> H + CH2CHO','O + C2H5 <=> H + CH3CHO','OH + HO2 <=> O2 + H2O','OH + CH3 => H2 + CH2O',
     'CH + H2 (+ M) <=> CH3 (+ M)','CH2 + O2 => 2 H + CO2','CH2 + O2 <=> O + CH2O','CH2 + CH2 => 2 H + C2H2',
     'CH2(S) + H2O => H2 + CH2O','C2H3 + O2 <=> O + CH2CHO','C2H3 + O2 <=> HO2 + C2H2','O + CH3CHO <=> OH + CH2CHO',
     'O + CH3CHO => OH + CH3 + CO','O2 + CH3CHO => HO2 + CH3 + CO','H + CH3CHO <=> CH2CHO + H2','H + CH3CHO => CH3 + H2 + CO',
     'OH + CH3CHO => CH3 + H2O + CO','HO2 + CH3CHO => CH3 + H2O2 + CO','CH3 + CH3CHO => CH3 + CH4 + CO','H + CH2CO (+ M) <=> CH2CHO (+ M)',
     'O + CH2CHO => H + CH2 + CO2','O2 + CH2CHO => OH + CO + CH2O','O2 + CH2CHO => OH + 2 HCO','H + CH2CHO <=> CH3 + HCO',
     'H + CH2CHO <=> CH2CO + H2','OH + CH2CHO <=> H2O + CH2CO','OH + CH2CHO <=> HCO + CH2OH','CH3 + C2H5 (+ M) <=> C3H8 (+ M)',
     'O + C3H8 <=> OH + C3H7','H + C3H8 <=> C3H7 + H2','OH + C3H8 <=> C3H7 + H2O','C3H7 + H2O2 <=> HO2 + C3H8',
     'CH3 + C3H8 <=> C3H7 + CH4','CH3 + C2H4 (+ M) <=> C3H7 (+ M)','O + C3H7 <=> C2H5 + CH2O','H + C3H7 (+ M) <=> C3H8 (+ M)',
     'H + C3H7 <=> CH3 + C2H5','OH + C3H7 <=> C2H5 + CH2OH','HO2 + C3H7 <=> O2 + C3H8','HO2 + C3H7 => OH + C2H5 + CH2O',
     'CH3 + C3H7 <=> 2 C2H5']



root = tk.Tk() 



#labels,Buttons,frames, optionmenus
topFrame = tk.Frame(root)
topFrame.pack()
botFrame = tk.Frame(root)
botFrame.pack()

        
label1 = tk.Label(topFrame, text="Temperature [K]", font=LARGE_FONT)
label1.grid(row=1,column=0)
var1=tk.StringVar()
set1 = tk.OptionMenu(topFrame,var1, T[0],T[1])#, T[2], T[3],T[4], T[5],
                            #T[6],T[7], T[8], T[9])
set1.configure(font=LARGE_FONT)
set1.grid(row=1,column=1)
       
label2 = tk.Label(topFrame, text="Pressure [atm]", font=LARGE_FONT)
label2.grid(row=1,column=3)
var2=tk.StringVar()
set2 = tk.OptionMenu(topFrame,var2, P[0], P[1])#, P[2], P[3], P[4], P[5],
                             #P[6], P[7], P[8], P[9])
set2.configure(font=LARGE_FONT)
set2.grid(row=1,column=4)
       
       
label3 = tk.Label(topFrame, text="Equivalence Ratio", font=LARGE_FONT)
label3.grid(row=2,column=0)
var3=tk.StringVar()
set3 = tk.OptionMenu(topFrame,var3, Phi[0], Phi[1])
set3.configure(font=LARGE_FONT)
set3.grid(row=2,column=1)

label4 = tk.Label(topFrame, text="Fuel", font=LARGE_FONT)
label4.grid(row=2,column=3)
var4=tk.StringVar()
set4 = tk.OptionMenu(topFrame,var4, CH4[0], CH4[1])
set4.configure(font=LARGE_FONT)
set4.grid(row=2,column=4)
     
label5 = tk.Label(topFrame, text="Species", font=LARGE_FONT)
label5.grid(row=3,column=0)

var5=tk.StringVar()
set5 = tk.OptionMenu(topFrame,var5, SS[0], SS[1], SS[2], SS[3], SS[4])
set5.configure(font=LARGE_FONT)
set5.grid(row=3,column=1)

var6=tk.StringVar()
set6 = tk.OptionMenu(topFrame,var6, S[0], S[1], S[2], S[3], S[4], S[5], S[6],
                     S[7], S[8], S[9], S[10], S[11], S[12], S[13], S[14], S[15],
                     S[16], S[17], S[18], S[19], S[20], S[21], S[22], S[23],
                     S[24], S[25], S[26], S[27], S[28], S[29], S[30], S[31], 
                     S[32], S[33], S[34], S[35], S[36], S[37], S[38], S[39],
                     S[40], S[41], S[42], S[43], S[44], S[45], S[46], S[47],
                     S[48], S[49], S[50], S[51], S[52])
set6.configure(font=LARGE_FONT)
set6.grid(row=3,column=2)

var7=tk.StringVar()
set7 = tk.OptionMenu(topFrame,var7, S[0], S[1], S[2], S[3], S[4], S[5], S[6],
                     S[7], S[8], S[9], S[10], S[11], S[12], S[13], S[14], S[15],
                     S[16], S[17], S[18], S[19], S[20], S[21], S[22], S[23],
                     S[24], S[25], S[26], S[27], S[28], S[29], S[30], S[31], 
                     S[32], S[33], S[34], S[35], S[36], S[37], S[38], S[39],
                     S[40], S[41], S[42], S[43], S[44], S[45], S[46], S[47],
                     S[48], S[49], S[50], S[51], S[52])
set7.configure(font=LARGE_FONT)
set7.grid(row=3,column=3)

var8=tk.StringVar()
set8 = tk.OptionMenu(topFrame,var8, S[0], S[1], S[2], S[3], S[4], S[5], S[6],
                     S[7], S[8], S[9], S[10], S[11], S[12], S[13], S[14], S[15],
                     S[16], S[17], S[18], S[19], S[20], S[21], S[22], S[23],
                     S[24], S[25], S[26], S[27], S[28], S[29], S[30], S[31], 
                     S[32], S[33], S[34], S[35], S[36], S[37], S[38], S[39],
                     S[40], S[41], S[42], S[43], S[44], S[45], S[46], S[47],
                     S[48], S[49], S[50], S[51], S[52])
set8.configure(font=LARGE_FONT)
set8.grid(row=3,column=4)

var9=tk.StringVar()
set9 = tk.OptionMenu(topFrame,var9, S[0], S[1], S[2], S[3], S[4], S[5], S[6],
                     S[7], S[8], S[9], S[10], S[11], S[12], S[13], S[14], S[15],
                     S[16], S[17], S[18], S[19], S[20], S[21], S[22], S[23],
                     S[24], S[25], S[26], S[27], S[28], S[29], S[30], S[31], 
                     S[32], S[33], S[34], S[35], S[36], S[37], S[38], S[39],
                     S[40], S[41], S[42], S[43], S[44], S[45], S[46], S[47],
                     S[48], S[49], S[50], S[51], S[52])
set9.configure(font=LARGE_FONT)
set9.grid(row=3,column=5)
        
label6 = tk.Label(topFrame, text="Rxn # (1 - 325)", font=LARGE_FONT)
label6.grid(row=4,column=0)

entry1=tk.Entry(topFrame,width=5)
entry1.grid(row=4,column=1)

entry2=tk.Entry(topFrame,width=5)
entry2.grid(row=4,column=2)

label7 = tk.Label(topFrame, text="Time [s], (1 - 5)", font=LARGE_FONT)
label7.grid(row=5,column=0)

entry3=tk.Entry(topFrame,width=5)
entry3.grid(row=5,column=1)

#functions
def do():
    Temp = float(var1.get())
    Press = float(var2.get())
    phi = float(var3.get())
    Fuel = float(var4.get())
    specie1 = var5.get() 
    specie2 = var6.get()
    specie3 = var7.get()
    specie4 = var8.get()
    specie5 = var9.get()
    species =[specie1,specie2,specie3,specie4,specie5]
    oxygen = 2/phi*Fuel
    nitrogen = 1 - oxygen - Fuel
    mixture={'CH4':Fuel, 'O2':oxygen, 'N2': nitrogen} 
    rxncount = int(entry1.get())-1
    entry2=tk.Entry(topFrame,width=5)
    entry2.grid(row=4,column=2)
    entryText = tk.StringVar()
    entryText.set(RR[rxncount])
    entry2 = tk.Entry(topFrame, textvariable = entryText)
    entry2.grid(row=4,column=2)
    
    for i in range(0,len(MoleFraction)):
        x=MoleFraction[i]
        x1=x[0]
        x2=x[1]
        x3=x[2]
        if mixture == x1 and Temp == x2 and Press == x3:
            mixcount = i
    speciecount=[0,0,0,0,0]
    for i in range(0, len(S)):
        if specie1 == S[i]:
            speciecount[0] = i
        if specie2 == S[i]:
            speciecount[1] = i
        if specie3 == S[i]:
            speciecount[2] = i
        if specie4 == S[i]:
            speciecount[3] = i
        if specie5 == S[i]:
            speciecount[4] = i
#    print(speciecount)
    print(mixcount)
#    plota(mixcount, speciecount,rxncount,specie1)
    plotb(mixcount, speciecount,rxncount,species)
    sensitive_rxn(mixcount,species,speciecount)
   
button1 = tk.Button(topFrame, text="Plot it!", command=do)
button1.grid(row=7, column=1)    

def sensitive_rxn(mixcount, species,speciecount):
    err = 100
    loc = 0
    time = float(entry3.get())
    y=All_tTP_AMo[mixcount]
#    y=All_tTP_SMo_AMo_SMa_AMa[mixcount]
    for i in range(0,len(y)-1):
        y1=y[i] #iterates through time steps of chosen simulation
        y11=y1[0] #grabs first item in above list (should be time)
        if np.absolute(y11 - time) <= err:
            loc = i
            err = np.absolute(y11 - time)
    print(loc) # prints location of least error
    

    specnum = 0
    specspecies=var5.get()
    for i in range (0, len(SS)):
        if specspecies == SS[i]:
            specnum = i
    print(specnum)
    print(loc+specnum*len(y))
    x=All_Ranks[mixcount]
    x1=x[loc*len(SS)+specnum]
    rxn_1=int(x1.index(1))
    rxn_2=int(x1.index(2))
    rxn_3=int(x1.index(3))
    rxn_4=int(x1.index(4))
    rxn_5=int(x1.index(5))
    rxn_=[rxn_1,rxn_2,rxn_3,rxn_4,rxn_5]
    rxn__ = [x+1 for x in rxn_] #use for labels of plot lines
    R=(rxn_,rxn__)
    print(rxn_)
    print(rxn__)
    plota(mixcount, species, speciecount, rxn_, rxn__)

    
def plotb(mixcount, speciecount,rxncount,species):
    time=[]
    molfrac1=[]
    molfrac2=[]
    molfrac3=[]
    molfrac4=[]
    molfrac5=[]
    y=All_tTP_AMo[mixcount]
#    y=All_tTP_SMo_AMo_SMa_AMa[mixcount]
    for i in range(0,len(y)-1):
        y1=y[i]
        y21=y1[0] #time
        y22_=np.absolute(y1[3])  #AMo
        y22_1=y22_[int(speciecount[0])]
        y22_2=y22_[int(speciecount[1])]
        y22_3=y22_[int(speciecount[2])]
        y22_4=y22_[int(speciecount[3])]
        y22_5=y22_[int(speciecount[4])]
        time.append(y21)
        molfrac1.append(y22_1)
        molfrac2.append(y22_2)
        molfrac3.append(y22_3)
        molfrac4.append(y22_4)
        molfrac5.append(y22_5)
    
    maxtime=max(time)
    maxX=-99999999999999999
    
    if max(molfrac1)> maxX:
        maxX=max(molfrac1)
    if max(molfrac2)> maxX:
        maxX=max(molfrac2)
    if max(molfrac3)> maxX:
        maxX=max(molfrac3)
    if max(molfrac4)> maxX:
        maxX=max(molfrac4)
    if max(molfrac5)> maxX:
        maxX=max(molfrac5)
        
    line1_1.set_xdata(time)
    line1_1.set_ydata(molfrac1)
    line1_1.set_label(species[0])
    line1_1
    
    line1_2.set_xdata(time)
    line1_2.set_ydata(molfrac2)
    line1_2.set_label(species[1])
    line1_2
    
    line1_3.set_xdata(time)
    line1_3.set_ydata(molfrac3)
    line1_3.set_label(species[2])
    line1_3
    
    line1_4.set_xdata(time)
    line1_4.set_ydata(molfrac4)
    line1_4.set_label(species[3])
    line1_4
    
    line1_5.set_xdata(time)
    line1_5.set_ydata(molfrac5)
    line1_5.set_label(species[4])
    line1_5
    
    b.legend()
    b.axis((0-0.1*maxtime, maxtime+0.1*maxtime, 0-0.1*maxX, maxX+0.1*maxX))
    canvas.draw()
    

def plota(mixcount, species, speciecount, rxn_, rxn__):
    time=[]
    sens=[]
    sens1=[]
    sens2=[]
    sens3=[]
    sens4=[]
    y=All_time_Sens[mixcount]
    for i in range(0,len(y)):
        y1=y[i]
        y21=y1[0] #time
        y22=y1[1] #sensitivity
        y22_=np.absolute(y22[int(speciecount[0]), int(rxn_[0])])
        y22_1=np.absolute(y22[int(speciecount[0]),int(rxn_[1])])
        y22_2=np.absolute(y22[int(speciecount[0]),int(rxn_[2])])
        y22_3=np.absolute(y22[int(speciecount[0]),int(rxn_[3])])
        y22_4=np.absolute(y22[int(speciecount[0]),int(rxn_[4])])
        time.append(y21)
        sens.append(y22_)
        sens1.append(y22_1)
        sens2.append(y22_2)
        sens3.append(y22_3)
        sens4.append(y22_4)
        
    maxtime=max(time)
    maxX=-99999999999999999
    
    if max(sens)> maxX:
        maxX=max(sens)
    if max(sens1)> maxX:
        maxX=max(sens1)
    if max(sens1)> maxX:
        maxX=max(sens2)
    if max(sens3)> maxX:
        maxX=max(sens3)
    if max(sens4)> maxX:
        maxX=max(sense4)
        
    line.set_xdata( time)
    line.set_ydata( sens )
    line.set_label(rxn__[0])
    line
    
    line1.set_xdata( time)
    line1.set_ydata( sens1 )
    line1.set_label(rxn__[1])
    line1
    
    line2.set_xdata( time)
    line2.set_ydata( sens2 )
    line2.set_label(rxn__[2])
    line2
    
    line3.set_xdata( time)
    line3.set_ydata( sens3 )
    line3.set_label(rxn__[3])
    line3
    
    line4.set_xdata( time)
    line4.set_ydata( sens4 )
    line4.set_label(rxn__[4])
    line4
    
    atitle= "Sensitivity Vs. t of " + species[0]
    a.set_title(atitle)
    a.legend()
    a.axis((0-0.1*max((time)), (max((time))+0.1*max((time))), 0-0.1*max((sens)), (max((sens))+0.1*max((sens)))))
    canvas.draw()

    


#figures
f = Figure(figsize=(11,4), dpi=100)
#first plot
a = f.add_subplot(121)
atitle= "Sensitivity Vs. t"
a.set_title(atitle)
a.set_xlabel('Time [s]')
a.set_ylabel('Sensitivity')
#second plot
b = f.add_subplot(122)
btitle= "X_i Vs. t"
b.set_title(btitle)
b.set_xlabel('Time [s]')
b.set_ylabel('X - Mole Fraction')

canvas = FigureCanvasTkAgg(f, botFrame)
canvas.get_tk_widget().pack(side=tk.BOTTOM, fill=tk.BOTH, expand = True)
#plota
line, = a.plot([], [], "r-")
line1, = a.plot([], [], "b--")
line2, = a.plot([], [], "k-.")
line3, = a.plot([], [], "g:")
line4, = a.plot([], [], "m|")

#plotb
line1_1, = b.plot([], [], "r-")
line1_2, = b.plot([], [], "b--")
line1_3, = b.plot([], [], "k-.")
line1_4, = b.plot([], [], "g:")
line1_5, = b.plot([], [], "m|")

        
toolbar = NavigationToolbar2TkAgg(canvas, botFrame)
toolbar.update()
canvas._tkcanvas.pack(side=tk.TOP, fill=tk.BOTH, expand = True)



root.mainloop()