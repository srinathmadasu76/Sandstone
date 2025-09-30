# -*- coding: utf-8 -*-
"""
Created on Fri Sep 26 17:02:01 2025

@author: srina
"""

# -*- coding: utf-8 -*-
"""
Created on Fri Sep  5 11:05:09 2025

@author: srina
"""

# -*- coding: utf-8 -*-
"""
Created on Mon Aug 18 10:00:36 2025

@author: srina
"""

# -*- coding: utf-8 -*-
"""
Created on Tue Feb 14 12:57:17 2023

@author: srina
"""
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
import numpy as np
import streamlit as st
import pandas as pd



kperm = 20 #mD
g = 3
time = 1 #day

st.header('Sandstone Acidizing Simulation - Linqx')   
      
por = st.number_input('Enter initial porosity: ', value=0.2)
kperm = st.number_input('Enter inital permeability mD: ', value=100)
#concsand = st.number_input('Enter inital sand concentration lb-mol/ft3: ')
#conckfelspar = st.number_input('Enter inital potassium felspar concentration lb-mol/ft3: ')
#concNafelspar = st.number_input('Enter inital sodium felspar tration lb-mol/ft3: ')
#concclay = st.number_input('Enter inital clay concentration lb-mol/ft3: ')

conchfwt = st.number_input('Enter HF wt %: ',min_value=0.,max_value=100.,step=1e-6,format="%.5f", value=3.)

conchf = conchfwt*1.17*10/20

conchclwt = st.number_input('Enter HCL wt %: ',min_value=0.,max_value=100.,step=1e-6,format="%.5f", value=10.)

conchcl = conchclwt*1.184*10/36.46

time = st.number_input('Enter time seconds: ', value=360)
Temperature = st.number_input('Enter Temperature F: ', value=100)

g = st.number_input('Enter permeability model exponent: ', value=3)
RockVolume = st.number_input('Enter rock volume cm3: ', value=1000.)

CalciteRockVolume = st.number_input('Enter calcite rock fraction volume : ', value=0.3)
QuartzRockVolume = st.number_input('Enter quatrz rock fraction volume : ', value=0.2)
KfelsparRockVolume = st.number_input('Enter k-felspar rock fraction volume : ', value=0.2)
ClayRockVolume = st.number_input('Enter clay rock fraction volume : ', value=0.3)

DolomiteRockVolume = st.number_input('Enter dolomite rock fraction volume : ', value=0.01)
SideriteRockVolume = st.number_input('Enter siderite rock fraction volume : ', value=0.01)
AnkeriteRockVolume = st.number_input('Enter ankerite rock fraction volume : ', value=0.01)
ZeoliteRockVolume = st.number_input('Enter ankerite rock fraction volume : ', value=0.01)
NafelsparRockVolume = st.number_input('Enter Na-felspar rock fraction volume : ', value=0.01)
KaoliniteRockVolume = st.number_input('Enter kaolinite rock fraction volume : ', value=0.01)
ChloriteRockVolume = st.number_input('Enter chlorite rock fraction volume : ', value=0.01)
MuscoviteRockVolume = st.number_input('Enter muscovite rock fraction volume : ', value=0.01)
MontmorilloniteRockVolume = st.number_input('Enter muscovite rock fraction volume : ', value=0.01)

R = 8.3144 #J/molK
T = 273+(Temperature-32)*5/9 #K

MolwtCalcite = 100.1 #g/mol
DensityCalcite = 2.71 #g/cm3
SurfaceareaCalcite = 0.38 #(m2/g)
kcoCalcite = 7.314*np.power(10,7) #mol/m2s
EaCalcite = 62.8 #KJ/mol


MolwtQuartz =  60.1 #g/mol
DensityQuartz = 2.65 #g/cm3
SurfaceareaQuartz = 0.05 #(m2/g)
kcoQuartz = 2.32*np.float_power(10,-8) #mol/m2s
print("kco=",kcoQuartz)
#kcoQuartz = 1*np.float_power(10,-4)/10000 #mol/m2s
EaQuartz = 9.56 #KJ/mol
#EaQuartz = 1.8932 #KJ/mol

Molwtkfelspar =  278.4 #g/mol
Densitykfelspar = 2.5 #g/cm3
Surfaceareakfelspar = 1.79 #(m2/g)
kcokfelspar = 0.127 #mol/m2s
Eakfelspar = 38.9 #KJ/mol

Molwtclay = 398 #g/mol
Densityclay = 2.75 #g/cm3
Surfaceareaclay = 37 #(m2/g)
kcoclay = 0.0275 #mol/m2s
Eaclay = 54.4 #KJ/mol

Molwthcl = 36.5 #g/mol
Molwthf = 20 #g/mol
Densityhcl = 1.07 #g/cm3
Densityhf = 1.07 #g/cm3


Molwth2sif6 = 144.09 #g/mol
Densityh2sif6 = 1.22 #g/cm3

Molwtsilica = 70 #g/mol
Densitysilica = 2.1 #g/cm3

Molwtsilicagel = 192 #g/mol
Densitysilicagel = 2.2 #g/cm3


Molwthalf2 = 64.97 #g/mol
Densityalf2 = 3.1 #g/cm3

Molwtcaf2= 78
Densitycaf2= 3.18

Molwtdolomite = 184.4 #g/mol
Densitydolomite = 2.85 #g/cm3

Molwtsiderite = 115.9 #g/mol
Densitysiderite = 3.85 #g/cm3

Molwtankerite = 200.2 #g/mol
Densityankerite = 3.35 #g/cm3

Molwtzeolite = 220.2 #g/mol
Densityzeolite = 2.3 #g/cm3

Molwtnafelspar =  262.2 #g/mol
Densitynafelspar = 2.62 #g/cm3
Surfaceareanafelspar = 1.79 #(m2/g)
kconafelspar = 0.127 #mol/m2s
Eanafelspar = 38.9 #KJ/mol

Molwtkaolinite =  258.1 #g/mol
Densitykaolinite = 2.59 #g/cm3

Molwtchlorite =  598.3 #g/mol
Densitychlorite = 2.88 #g/cm3

Molwtmuscovite =  381.3 #g/mol
Densitymuscovite = 2.6 #g/cm3

Molwtmontmorillonite =  429.6 #g/mol
Densitymontmorillonite = 2.6 #g/cm3

Molwtnasif6 =  188.06 #g/mol
Densitynasif6 = 2.68 #g/cm3


Molwtksif6 =  220.27 #g/mol
Densityksif6 = 2.7 #g/cm3

kcalcite_hcl = SurfaceareaCalcite*kcoCalcite*np.exp(-EaCalcite*1000/(R*T))
kcalcite_hf = SurfaceareaCalcite*0.001*np.exp(-40*1000/(R*T))
kquartz_hf = SurfaceareaQuartz*kcoQuartz*np.exp(-EaQuartz*1000/(R*T))
kfelspar_hf = Surfaceareakfelspar*kcokfelspar*np.exp(-Eakfelspar*1000/(R*T))*(1+ 0.0566*np.exp(956/T))
kclay_hf = Surfaceareaclay*kcoclay*np.exp(-Eaclay*1000/(R*T))
kfelsparh2sif6 = 122775*np.exp(-14*4184/(R*T))
kclayh2sif6=7.82*np.power(10,22)*np.exp(-37.8*4184/(R*T))
kalf2_hcl = 5.525*np.exp(-8.12*4184/(8.314*T))
kdolomite_hcl = 1.72*np.power(10,11)*np.exp(-15*4184/(R*T))
kdolomite_hf = 1.72*np.power(10,9)*np.exp(-15*4184/(R*T))
ksiderite_hcl = 1.72*np.power(10,10)*np.exp(-15*4184/(R*T))
kankerite_hcl = 1.72*np.power(10,11)*np.exp(-15*4184/(R*T))
kankerite_hf = 1.72*np.power(10,9)*np.exp(-15*4184/(R*T))
ksilica_hf = 3.12*np.power(10,6)*np.exp(-9.8*4184/(R*T))
kzeolite_hf = 9.4*np.power(10,4)*np.exp(-7.8*4184/(R*T))
kzeolite_h2sif6 = 4.14*np.power(10,5)*np.exp(-7.8*4184/(R*T))
nafelspar_hf = Surfaceareanafelspar*kconafelspar*np.exp(-Eanafelspar*1000/(R*T))*(1+ 0.0566*np.exp(956/T))
nafelsparh2sif6 = 122775*np.exp(-14*4184/(R*T))
kkaolinite_hf = 4.7*np.power(10,4)*np.exp(-7.8*4184/(R*T))
kkaolinite_h2sif6 = 2.07*np.power(10,5)*np.exp(-7.8*4184/(R*T))
kchlorite_hf = 9.4*np.power(10,4)*np.exp(-7.8*4184/(R*T))
kchlorite_h2sif6 = 4.14*np.power(10,5)*np.exp(-7.8*4184/(R*T))
kmuscovite_hf = 9.4*np.power(10,4)*np.exp(-7.8*4184/(R*T))
kmuscovite_h2sif6 = 4.14*np.power(10,5)*np.exp(-7.8*4184/(R*T))
kmontmorillonite_hf = 1.88*np.power(10,5)*np.exp(-7.8*4184/(R*T))
kmontmorillonite_h2sif6 = 8.28*np.power(10,5)*np.exp(-7.8*4184/(R*T))

if conchf==0.:
    multhf = 0.
else:
    multhf = 1
if conchcl==0.:
    multhcl = 0.
else:
    multhcl = 1 
    
if CalciteRockVolume==0.:
    multcalcite = 0.
else:
    multcalcite = 1
    
if QuartzRockVolume==0.:
    multquartz = 0.
else:
    multquartz = 1

if KfelsparRockVolume==0.:
    multkfelspar = 0.
else:
    multkfelspar = 1
    
if ClayRockVolume==0.:
    multclay = 0.
else:
    multclay = 1
    
if DolomiteRockVolume==0.:
    multdolomite = 0.
else:
    multdolomite = 1
    
if SideriteRockVolume==0.:
    multisiderite = 0.
else:
    multisiderite = 1

if AnkeriteRockVolume==0.:
    multankerite = 0.
else:
    multankerite = 1
    
if ZeoliteRockVolume==0.:
    multzeolite = 0.
else:
    multzeolite = 1   
    
if NafelsparRockVolume==0.:
    multnafelspar = 0.
else:
    multnafelspar = 1

if KaoliniteRockVolume==0.:
    multnakaolinite = 0.
else:
    multkaolinite = 1
    
if ChloriteRockVolume==0.:
    multnachlorite = 0.
else:
    multchlorite = 1

if MuscoviteRockVolume==0.:
    multnmuscovite = 0.
else:
    multmuscovite = 1
    
if MontmorilloniteRockVolume==0.:
    multmontmorillonite = 0.
else:
    multmontmorillonite = 1
    
#print(rquartz)
#print("hf=",hftest)
KspNa = 19189*np.exp(-45001.5/(R*T))
KspK = 197000*np.exp(-62726/(R*T))

if st.button('Calculate'):
    
    def rhs(s, v): 
              
          #rcalcite = -SurfaceareaCalcite*kcoCalcite*np.exp(-EaCalcite*1000/(R*T))*np.power(v[2],0.63)*MolwtCalcite/DensityCalcite
          rhcl = -kcalcite_hcl*multcalcite*np.power(v[0],0.63)*2-kfelspar_hf*multkfelspar*(1+np.power(v[0],0.4))*np.power(v[1],1.2)*3-18*kfelsparh2sif6*multkfelspar*np.power(v[6],0.6)*(v[4])*np.power(v[4]*Molwtkfelspar/(Molwtkfelspar*v[4]+v[7]*Molwtsilica),0.1)-6*kclayh2sif6*multclay*np.power(v[6],0.6)*(v[5])*np.power(v[5]*Molwtclay/(Molwtclay*v[3]+v[7]*Molwtsilica),0.1)-4*kalf2_hcl*multkfelspar*v[8]*np.power(v[0],1.7)*(v[4])*np.power(v[4]*Molwtkfelspar/(Molwtkfelspar*v[4]+v[7]*Molwtsilica),0.1)+4*(-kdolomite_hf*multdolomite*v[1]-kdolomite_hcl*multdolomite*np.power(v[0],1))*Densitydolomite/Molwtdolomite*multcalcite-2*ksiderite_hcl*multisiderite*np.power(v[0],1)*Densitysiderite/Molwtsiderite*multcalcite+4*(-kankerite_hf*multankerite*v[1]-kankerite_hcl*multankerite*np.power(v[0],1))*Densityankerite/Molwtankerite*multcalcite-nafelspar_hf*multnafelspar*(1+np.power(v[0],0.4))*np.power(v[1],1.2)*3-18*nafelsparh2sif6*multnafelspar*np.power(v[6],0.6)*(v[4])*np.power(v[4]*Molwtnafelspar/(Molwtnafelspar*v[4]+v[7]*Molwtsilica),0.1)-4*kalf2_hcl*multnafelspar*v[8]*np.power(v[0],1.7)*(v[4])*np.power(v[4]*Molwtnafelspar/(Molwtnafelspar*v[4]+v[7]*Molwtsilica),0.1)
          rhf = -kcalcite_hf*multcalcite*v[1]*2-kquartz_hf*multquartz*v[1]*6-kfelspar_hf*multkfelspar*(1+np.power(v[0],0.4))*np.power(v[1],1.2)*19-kclay_hf*v[1]*24-ksilica_hf*v[1]*6*Densitysilica/Molwtsilica-nafelspar_hf*multnafelspar*(1+np.power(v[0],0.4))*np.power(v[1],1.2)*19+(-32*kkaolinite_hf*multkaolinite*v[1]+6*kkaolinite_h2sif6*multkaolinite*np.power(v[6],1))*Densitykaolinite/Molwtkaolinite*multkaolinite+(-30*kchlorite_hf*multchlorite*v[1])*Densitychlorite/Molwtchlorite*multchlorite+(-36*kmuscovite_hf*multchlorite*v[1])*Densitymuscovite/Molwtmuscovite*multmuscovite+(-36*kmontmorillonite_hf*multmontmorillonite*v[1])*Densitymontmorillonite/Molwtmontmorillonite*multmontmorillonite
          rcalcite = -kcalcite_hf*multcalcite*v[1]-kcalcite_hcl*multcalcite*np.power(v[0],0.63)
          rquartz = -kquartz_hf*multquartz*v[1]
          rkfelspar = -kfelspar_hf*multkfelspar*(1+np.power(v[0],0.4))*np.power(v[1],1.2)-6*kfelsparh2sif6*multkfelspar*np.power(v[6],0.6)*(v[4])*np.power(v[4]*Molwtkfelspar/(Molwtkfelspar*v[4]+v[7]*Molwtsilica),0.1)-kalf2_hcl*multkfelspar*v[8]*np.power(v[0],1.7)*(v[4])*np.power(v[4]*Molwtkfelspar/(Molwtkfelspar*v[4]+v[7]*Molwtsilica),0.1)
          rclay = -kclay_hf*v[1]+kclayh2sif6*multclay*np.power(v[6],0.6)*(v[5])*np.power(v[5]*Molwtclay/(Molwtclay*v[3]+v[7]*Molwtsilica),0.1)
          rh2sif6felspar = -6*kfelsparh2sif6*multkfelspar*np.power(v[6],0.6)*(v[4])*np.power(v[4]*Molwtkfelspar/(Molwtkfelspar*v[4]+v[7]*Molwtsilica),0.1)
          rsilica = kfelsparh2sif6*multkfelspar*np.power(v[6],0.6)*(v[4])*np.power(v[4]*Molwtkfelspar/(Molwtkfelspar*v[4]+v[7]*Molwtsilica),0.1)+13*kclayh2sif6*multclay*np.power(v[6],0.6)*(v[5])*np.power(v[5]*Molwtclay/(Molwtclay*v[5]+v[7]*Molwtsilica),0.1)-ksilica_hf*v[1]*Densitysilica/Molwtsilica+2*kzeolite_h2sif6*multzeolite*np.power(v[6],1)*Densityzeolite/Molwtzeolite*multzeolite+nafelsparh2sif6*multnafelspar*np.power(v[6],0.6)*(v[4])*np.power(v[4]*Molwtnafelspar/(Molwtnafelspar*v[4]+v[7]*Molwtsilica),0.1)+6*kkaolinite_h2sif6*multkaolinite*np.power(v[6],1)*Densitykaolinite/Molwtkaolinite*multkaolinite+5*kchlorite_h2sif6*multchlorite*np.power(v[6],1)*Densitychlorite/Molwtchlorite*multchlorite+(6*kmuscovite_h2sif6*multchlorite*v[6])*Densitymuscovite/Molwtmuscovite*multmuscovite+6*kmontmorillonite_h2sif6*multmontmorillonite*np.power(v[6],1)*Densitymontmorillonite/Molwtmontmorillonite*multmontmorillonite
          rh2sif6clay = -kclayh2sif6*multclay*np.power(v[6],0.6)*(v[5])*np.power(v[5]*Molwtclay/(Molwtclay*v[5]+v[7]*Molwtsilica),0.1)
          
          rnafelspar = -nafelspar_hf*multnafelspar*(1+np.power(v[0],0.4))*np.power(v[1],1.2)-6*nafelsparh2sif6*multnafelspar*np.power(v[6],0.6)*(v[4])*np.power(v[4]*Molwtnafelspar/(Molwtnafelspar*v[4]+v[7]*Molwtsilica),0.1)-kalf2_hcl*multnafelspar*v[8]*np.power(v[0],1.7)*(v[4])*np.power(v[4]*Molwtnafelspar/(Molwtnafelspar*v[4]+v[7]*Molwtsilica),0.1)
          rh2sif6nafelspar = -6*nafelsparh2sif6*multnafelspar*np.power(v[6],0.6)*(v[4])*np.power(v[4]*Molwtnafelspar/(Molwtnafelspar*v[4]+v[7]*Molwtsilica),0.1)
          
          rh2sif6 = -rquartz-3*rclay-3*rkfelspar+rh2sif6felspar/6+12*rh2sif6clay+ksilica_hf*v[1]*Densitysilica/Molwtsilica-3*rnafelspar+rh2sif6nafelspar/6+kchlorite_h2sif6*multchlorite*np.power(v[6],1)*Densitychlorite/Molwtchlorite*multchlorite+2*kmontmorillonite_h2sif6*multmontmorillonite*np.power(v[6],1)*Densitymontmorillonite/Molwtmontmorillonite*multmontmorillonite
          ralf2 = 6*kfelsparh2sif6*multkfelspar*np.power(v[6],0.6)*(v[4])*np.power(v[4]*Molwtkfelspar/(Molwtkfelspar*v[4]+v[7]*Molwtsilica),0.1)+kfelspar_hf*multkfelspar*(1+np.power(v[0],0.4))*np.power(v[1],1.2)-kalf2_hcl*multkfelspar*v[8]*np.power(v[0],1.7)*(v[4])*np.power(v[4]*Molwtkfelspar/(Molwtkfelspar*v[4]+v[7]*Molwtsilica),0.1)+6*nafelsparh2sif6*multnafelspar*np.power(v[6],0.6)*(v[4])*np.power(v[4]*Molwtnafelspar/(Molwtnafelspar*v[4]+v[7]*Molwtsilica),0.1)+nafelspar_hf*multnafelspar*(1+np.power(v[0],0.4))*np.power(v[1],1.2)-kalf2_hcl*multnafelspar*v[8]*np.power(v[0],1.7)*(v[4])*np.power(v[4]*Molwtnafelspar/(Molwtnafelspar*v[4]+v[7]*Molwtsilica),0.1)+(4*kkaolinite_hf*multkaolinite*v[1]+3*kkaolinite_h2sif6*multkaolinite*np.power(v[6],1))*Densitykaolinite/Molwtkaolinite*multkaolinite
          rsilicagel = 3*kalf2_hcl*multkfelspar*v[8]*np.power(v[0],1.7)*(v[4])*np.power(v[4]*Molwtkfelspar/(Molwtkfelspar*v[4]+v[7]*Molwtsilica),0.1)+3*kalf2_hcl*multnafelspar*v[8]*np.power(v[0],1.7)*(v[4])*np.power(v[4]*Molwtnafelspar/(Molwtnafelspar*v[4]+v[7]*Molwtsilica),0.1)
          rcaf2 = kcalcite_hf*multcalcite*v[1]
         
          rdolomite = (-kdolomite_hf*multdolomite*v[1]-kdolomite_hcl*multdolomite*np.power(v[0],1))*Densitydolomite/Molwtdolomite*multcalcite
          rsiderite = -ksiderite_hcl*multisiderite*np.power(v[0],1)*Densitysiderite/Molwtsiderite*multcalcite
          rankerite = (-kankerite_hf*multankerite*v[1]-kankerite_hcl*multankerite*np.power(v[0],1))*Densityankerite/Molwtankerite*multcalcite
          rzeolite = (-kzeolite_hf*multzeolite*v[1]-kzeolite_h2sif6*multzeolite*np.power(v[6],1))*Densityzeolite/Molwtzeolite*multzeolite
          rkaolinite = (-kkaolinite_hf*multkaolinite*v[1]-kkaolinite_h2sif6*multkaolinite*np.power(v[6],1))*Densitykaolinite/Molwtkaolinite*multkaolinite
          rchlorite = (-kchlorite_hf*multchlorite*v[1]-kchlorite_h2sif6*multchlorite*np.power(v[6],1))*Densitychlorite/Molwtchlorite*multchlorite
          rmuscovite = (-kmuscovite_hf*multchlorite*v[1]-kmuscovite_h2sif6*multmuscovite*np.power(v[6],1))*Densitymuscovite/Molwtmuscovite*multmuscovite
          rmontmorillonite = (-kmontmorillonite_hf*multmontmorillonite*v[1]-kmontmorillonite_h2sif6*multmontmorillonite*np.power(v[6],1))*Densitymontmorillonite/Molwtmontmorillonite*multmontmorillonite
          
          rNa = nafelspar_hf*multnafelspar*(1+np.power(v[0],0.4))*np.power(v[1],1.2)+6*nafelsparh2sif6*multnafelspar*np.power(v[6],0.6)*(v[4])*np.power(v[4]*Molwtnafelspar/(Molwtnafelspar*v[4]+v[7]*Molwtsilica),0.1)+kalf2_hcl*multnafelspar*v[8]*np.power(v[0],1.7)*(v[4])*np.power(v[4]*Molwtnafelspar/(Molwtnafelspar*v[4]+v[7]*Molwtsilica),0.1)+(kzeolite_hf*multzeolite*v[1]+kzeolite_h2sif6*multzeolite*np.power(v[6],1))*Densityzeolite/Molwtzeolite*multzeolite
          Q1 = ((v[20] + v[6])+ np.sqrt((v[20] + v[6])*(v[20] + v[6])-4*(v[20]*v[6]-KspNa)))/2.
          Q2 = ((v[20] + v[6])- np.sqrt((v[20] + v[6])*(v[20] + v[6])-4*(v[20]*v[6]-KspNa)))/2. 
          Q = Q1
          if Q2<Q1:
              Q = Q2
          if v[21]<=Q:
              rNasif6 = 3864.556*np.exp(-2393.34/R*T)*np.power(Q,2)*(1+v[21])
          else:
              rNasif6 = 0.
              
          rK = kfelspar_hf*multkfelspar*(1+np.power(v[0],0.4))*np.power(v[1],1.2)+6*kfelsparh2sif6*multkfelspar*np.power(v[6],0.6)*(v[4])*np.power(v[4]*Molwtkfelspar/(Molwtkfelspar*v[4]+v[7]*Molwtsilica),0.1)+kalf2_hcl*multkfelspar*v[8]*np.power(v[0],1.7)*(v[4])*np.power(v[4]*Molwtkfelspar/(Molwtkfelspar*v[4]+v[7]*Molwtsilica),0.1)+(kmuscovite_hf*multchlorite*v[1]+kmuscovite_h2sif6*multmuscovite*np.power(v[6],1))*Densitymuscovite/Molwtmuscovite*multmuscovite
         
          Q1 = ((v[22] + v[6])+ np.sqrt((v[22] + v[6])*(v[22] + v[6])-4*(v[22]*v[6]-KspK)))/2.
          Q2 = ((v[22] + v[6])- np.sqrt((v[22] + v[6])*(v[22] + v[6])-4*(v[22]*v[6]-KspK)))/2. 
          Q = Q1
          if Q2<Q1:
              Q = Q2
          if v[23]<=Q:
              rKsif6 = 3864.556*np.exp(-2393.34/R*T)*np.power(Q,2)*(1+v[23])
          else:
              rKsif6 = 0.
          
          return [rhcl,rhf,rcalcite,rquartz,rkfelspar,rclay,rh2sif6,rsilica,ralf2,rsilicagel,rcaf2,rdolomite,rsiderite,rankerite,rzeolite,rnafelspar,rkaolinite,rchlorite,rmuscovite,rmontmorillonite,rNa,rNasif6,rK,rKsif6]
          #return [rcalcite,rquartz,rkfelspar,rclay,rh2sif6felspar,ralf2,rh2sif6clay,rcalcitehf]
      
    res = solve_ivp(rhs, (0, time), [conchcl, conchf, RockVolume*CalciteRockVolume/MolwtCalcite*DensityCalcite, RockVolume*QuartzRockVolume/MolwtQuartz*DensityQuartz,RockVolume*KfelsparRockVolume*Densitykfelspar/Molwtkfelspar,RockVolume*ClayRockVolume/Molwtclay*Densityclay,0.,0.,0.,0.,0., RockVolume*DolomiteRockVolume/Molwtdolomite*Densitydolomite,RockVolume*SideriteRockVolume/Molwtsiderite*Densitysiderite,RockVolume*AnkeriteRockVolume/Molwtankerite*Densityankerite,RockVolume*ZeoliteRockVolume/Molwtzeolite*Densityzeolite,RockVolume*NafelsparRockVolume*Densitynafelspar/Molwtnafelspar,RockVolume*KaoliniteRockVolume/Molwtkaolinite*Densitykaolinite,RockVolume*ChloriteRockVolume/Molwtchlorite*Densitychlorite,RockVolume*MuscoviteRockVolume/Molwtmuscovite*Densitymuscovite,RockVolume*MontmorilloniteRockVolume/Molwtmontmorillonite*Densitymontmorillonite,0.,0.,0.,0.],t_eval=np.linspace(0, time,10))
   
    #res = solve_ivp(rhs, (0, time), [RockVolume*CalciteRockVolume, RockVolume*QuartzRockVolume,RockVolume*KfelsparRockVolume,RockVolume*ClayRockVolume,0.,0.,0.,0.],t_eval=np.linspace(0, time,10))
    print("res0=",np.power(res.y[4],1.))
    print("res0=",np.power(res.y[0],0.63))
    print(list(np.linspace(0, time,10)))
    csilica = 0.
    csilicagel = 0.
    # ch2sif6 = res.y[4]*Molwth2sif6/Densityh2sif6*DensityQuartz/MolwtQuartz+res.y[4]*Molwth2sif6/Densityh2sif6*4*Densityclay/Molwtclay+res.y[4]*Molwth2sif6/Densityh2sif6*19*Densitykfelspar/Molwtkfelspar
    # csilica = -ch2sif6*Molwthsilica/Densitysilica*Densityh2sif6/Molwth2sif6*5
    # print("ch2sif6=",csilica)
    # csilicagel = -res.y[5]*Densityalf2/Molwthalf2*3*Molwtsilicagel/Densitysilicagel
    # ccaf2 = res.y[7]*Molwtcaf2/Densitycaf2*DensityCalcite/MolwtCalcite
    # print("ccaf2=",ccaf2)
    # print(res.y[0])
    df = pd.DataFrame({'Calcite RockVolume': res.y[2]*MolwtCalcite/DensityCalcite, 'time (seconds)': list(np.linspace(0, time,10))},columns=["Calcite RockVolume","time (seconds)"])
    st.line_chart(df,x="time (seconds)",y="Calcite RockVolume")
    df = pd.DataFrame({'Quartz RockVolume': res.y[3]*MolwtQuartz/DensityQuartz, 'time (seconds)': list(np.linspace(0, time,10))},columns=["Quartz RockVolume","time (seconds)"])
    st.line_chart(df,x="time (seconds)",y="Quartz RockVolume")

    df = pd.DataFrame({'K-Felspar RockVolume': res.y[4]/Densitykfelspar*Molwtkfelspar, 'time (seconds)': list(np.linspace(0, time,10))},columns=["K-Felspar RockVolume","time (seconds)"])
    st.line_chart(df,x="time (seconds)",y="K-Felspar RockVolume")
    df = pd.DataFrame({'Clay RockVolume': res.y[5]*Molwtclay/Densityclay, 'time (seconds)': list(np.linspace(0, time,10))},columns=["Clay RockVolume","time (seconds)"])
    st.line_chart(df,x="time (seconds)",y="Clay RockVolume")
    df = pd.DataFrame({'Concentration Silica': res.y[7]*Molwtsilica/Densitysilica, 'time (seconds)': list(np.linspace(0, time,10))},columns=["Concentration Silica","time (seconds)"])
    st.line_chart(df,x="time (seconds)",y="Concentration Silica")
    df = pd.DataFrame({'Silicagel': res.y[9]*Molwtsilicagel/Densitysilicagel, 'time (seconds)': list(np.linspace(0, time,10))},columns=["Silicagel","time (seconds)"])
    st.line_chart(df,x="time (seconds)",y="Silicagel")
    df = pd.DataFrame({'Calcium Fluoride': res.y[10]*Molwtcaf2/Densitycaf2, 'time (seconds)': list(np.linspace(0, time,10))},columns=["Calcium Fluoride","time (seconds)"])
    st.line_chart(df,x="time (seconds)",y="Calcium Fluoride")
    df = pd.DataFrame({'Porosity': (por+(-res.y[19]*Molwtmontmorillonite/Densitymontmorillonite+RockVolume*MontmorilloniteRockVolume-res.y[18]*Molwtmuscovite/Densitymuscovite+RockVolume*MuscoviteRockVolume-res.y[17]*Molwtchlorite/Densitychlorite+RockVolume*ChloriteRockVolume-res.y[16]*Molwtkaolinite/Densitykaolinite+RockVolume*KaoliniteRockVolume+RockVolume*NafelsparRockVolume-res.y[15]/Densitynafelspar*Molwtnafelspar-res.y[14]*Molwtzeolite/Densityzeolite+RockVolume*ZeoliteRockVolume-res.y[13]*Molwtankerite/Densityankerite+RockVolume*AnkeriteRockVolume-res.y[12]*Molwtsiderite/Densitysiderite+RockVolume*SideriteRockVolume-res.y[11]*Molwtdolomite/Densitydolomite+RockVolume*DolomiteRockVolume-res.y[2]*MolwtCalcite/DensityCalcite+RockVolume*CalciteRockVolume-res.y[3]*MolwtQuartz/DensityQuartz+RockVolume*KfelsparRockVolume+RockVolume*QuartzRockVolume-res.y[3]*MolwtQuartz/DensityQuartz+RockVolume*ClayRockVolume-res.y[4]/Densitykfelspar*Molwtkfelspar-res.y[7]*Molwtsilica/Densitysilica-res.y[9]*Molwtsilicagel/Densitysilicagel-res.y[10]*Molwtcaf2/Densitycaf2-res.y[21]*Molwtnasif6/Densitynasif6-res.y[23]*Molwtksif6/Densityksif6)/RockVolume), 'time (seconds)': list(np.linspace(0, time,10))},columns=["Porosity","time (seconds)"])
# # #print(np.vstack(((por-res.y[2]*aspmolwt/aspdens), np.linspace(0, time,10))))
    st.line_chart(df,x="time (seconds)",y="Porosity")
    df = pd.DataFrame({'Permeability': (kperm*np.power((por+(-res.y[19]*Molwtmontmorillonite/Densitymontmorillonite+RockVolume*MontmorilloniteRockVolume-res.y[18]*Molwtmuscovite/Densitymuscovite+RockVolume*MuscoviteRockVolume-res.y[17]*Molwtchlorite/Densitychlorite+RockVolume*ChloriteRockVolume-res.y[16]*Molwtkaolinite/Densitykaolinite+RockVolume*KaoliniteRockVolume+RockVolume*NafelsparRockVolume-res.y[15]/Densitynafelspar*Molwtnafelspar-res.y[14]*Molwtzeolite/Densityzeolite+RockVolume*ZeoliteRockVolume-res.y[13]*Molwtankerite/Densityankerite+RockVolume*AnkeriteRockVolume-res.y[12]*Molwtsiderite/Densitysiderite+RockVolume*SideriteRockVolume-res.y[11]*Molwtdolomite/Densitydolomite+RockVolume*DolomiteRockVolume-res.y[2]*MolwtCalcite/DensityCalcite+RockVolume*CalciteRockVolume-res.y[3]*MolwtQuartz/DensityQuartz+RockVolume*KfelsparRockVolume+RockVolume*QuartzRockVolume-res.y[3]*MolwtQuartz/DensityQuartz+RockVolume*ClayRockVolume-res.y[4]/Densitykfelspar*Molwtkfelspar-res.y[7]*Molwtsilica/Densitysilica-res.y[9]*Molwtsilicagel/Densitysilicagel-res.y[10]*Molwtcaf2/Densitycaf2-res.y[21]*Molwtnasif6/Densitynasif6-res.y[23]*Molwtksif6/Densityksif6)/RockVolume)/(por),g)), 'time (seconds)': list(np.linspace(0, time,10))},columns=["Permeability","time (seconds)"])
    st.line_chart(df,x="time (seconds)",y="Permeability")
    # st.line_chart(res.y[0])
    # st.line_chart(res.y[1])
    # st.line_chart(res.y[2])
    # print(res.y[2])
    # print(por-res.y[2]*aspmolwt/aspdens)
    # print(kperm*np.power((por-res.y[2]*aspmolwt/aspdens)/(por-asppor),g))
    # print(res.t)
    
    #plt.plot(res.t, res.y.T)