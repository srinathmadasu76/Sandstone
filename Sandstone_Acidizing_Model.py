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

aspdens = 75 #lb/ft3
aspmolwt = 1.54 #lb/mol
kcf = 12   #1/day
kinff = 20 #1/day
kinfd = 20 #1/day
ksf = 0.05 #lb-mol/ft3
ksd = 0.05 #lb-mol/ft3
alpha = 8#1/day
por = 0.1
concd = 2.5

K = 0.1/0.000694444#1/day
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

conchf = st.number_input('Enter HF concentration mol/Kgh20: ',min_value=0.,max_value=10000.,step=1e-6,format="%.5f", value=3.)

conchcl = st.number_input('Enter HCL concentration mol/Kgh20: ',min_value=0.,max_value=10000.,step=1e-6,format="%.5f", value=10.)

time = st.number_input('Enter time day: ', value=1)
Temperature = st.number_input('Enter Temperature F: ', value=100)

g = st.number_input('Enter permeability model exponent: ', value=3)
RockVolume = st.number_input('Enter rock volume cm3: ', value=1.)

CalciteRockVolume = st.number_input('Enter calcite rock fraction volume : ', value=0.1)
QuartzRockVolume = st.number_input('Enter quatrz rock fraction volume : ', value=0.2)
KfelsparRockVolume = st.number_input('Enter k-felspar rock fraction volume : ', value=0.4)
ClayRockVolume = st.number_input('Enter clay rock fraction volume : ', value=0.3)

R = 8.3144 #J/molK
T = 273+(Temperature-32)*5/9 #K

MolwtCalcite = 100.1 #g/mol
DensityCalcite = 2.71 #g/cm3
SurfaceareaCalcite = 0.38*10000 #(cm2/g)
kcoCalcite = 7.314*np.power(10,7)/10000 #mol/cm2s
EaCalcite = 62.8 #KJ/mol


MolwtQuartz =  60.1 #g/mol
DensityQuartz = 2.65 #g/cm3
SurfaceareaQuartz = 0.05*10000 #(cm2/g)
kcoQuartz = 2.32*np.float_power(10,-8)/10000 #mol/m2s
print("kco=",kcoQuartz)
#kcoQuartz = 1*np.float_power(10,-4)/10000 #mol/m2s
EaQuartz = 9.56 #KJ/mol
#EaQuartz = 1.8932 #KJ/mol

Molwtkfelspar =  278.4 #g/mol
Densitykfelspar = 2.5 #g/cm3
Surfaceareakfelspar = 1.79*10000 #(cm2/g)
kcokfelspar = 0.127 #mol/m2s
Eakfelspar = 38.9 #KJ/mol

Molwtclay = 398 #g/mol
Densityclay = 2.75 #g/cm3
Surfaceareaclay = 37*10000 #(cm2/g)
kcoclay = 0.0275 #mol/m2s
Eaclay = 54.4 #KJ/mol

Molwthcl = 36.5 #g/mol
Molwthf = 20 #g/mol
Densityhcl = 1.07 #g/cm3
Densityhf = 1.07 #g/cm3


Molwth2sif6 = 144.09 #g/mol
Densityh2sif6 = 1.22 #g/cm3

Molwthsilica = 70 #g/mol
Densitysilica = 2.1 #g/cm3

Molwtsilicagel = 192 #g/mol
Densitysilicagel = 2.2 #g/cm3


Molwthalf2 = 64.97 #g/mol
Densityalf2 = 3.1 #g/cm3

Molwtcaf2= 78
Densitycaf2= 3.18
rhcl = conchcl 
mult = 1
rquartz = SurfaceareaQuartz*kcoQuartz*np.exp(-EaQuartz*1000/(R*T))
hftest = (conchf-RockVolume*QuartzRockVolume*Molwthf/Densityhf*6*DensityQuartz/MolwtQuartz-RockVolume*QuartzRockVolume*Molwthf/Densityhf*36*Densityclay/Molwtclay-RockVolume*QuartzRockVolume*Molwthf/Densityhf*19*Densitykfelspar/Molwtkfelspar)*MolwtQuartz/DensityQuartz
#print(rquartz)
#print("hf=",hftest)
if st.button('Calculate'):
    asppor= conchf*aspmolwt/aspdens
    
    def rhs(s, v): 
         
          if(s==0):
              arg1=0
          else:
              arg1=1
          print("mult",arg1)
          #rcalcite = -SurfaceareaCalcite*kcoCalcite*np.exp(-EaCalcite*1000/(R*T))*np.power(v[2],0.63)*MolwtCalcite/DensityCalcite
          rhcl = conchcl -v[0]*Molwthcl/Densityhcl*2*DensityCalcite/MolwtCalcite*arg1
          rhf = conchf -v[1]*Molwthf/Densityhf*6*DensityQuartz/MolwtQuartz-v[1]*Molwthf/Densityhf*36*Densityclay/Molwtclay-v[1]*Molwthf/Densityhf*19*Densitykfelspar/Molwtkfelspar*arg1
          rcalcitehf = 0.001*100*np.exp(-40*1000/(8.314*T))*rhf
          rcalcite = -rcalcitehf-SurfaceareaCalcite*kcoCalcite*np.exp(-EaCalcite*1000/(R*T))*np.power(conchcl-v[0]*Molwthcl/Densityhcl*2*DensityCalcite/MolwtCalcite,0.63)*MolwtCalcite/DensityCalcite
         
         
          rquartz = -SurfaceareaQuartz*kcoQuartz*np.exp(-EaQuartz*1000/(R*T))*(conchf -v[1]*Molwthf/Densityhf*6*DensityQuartz/MolwtQuartz-v[1]*Molwthf/Densityhf*36*Densityclay/Molwtclay-v[1]*Molwthf/Densityhf*19*Densitykfelspar/Molwtkfelspar)*MolwtQuartz/DensityQuartz
          #print(rquartz)
          
         
          rkfelspar = -Surfaceareakfelspar*kcokfelspar*np.exp(-Eakfelspar*1000/(R*T))*(1+ 0.0566*np.exp(956/T))*(1+np.power(rhcl,0.4))*np.power(rhf,1.2)*Molwtkfelspar/Densitykfelspar
          rclay = -Surfaceareaclay*kcoclay*np.exp(-Eaclay*1000/(R*T))*rhf*Molwtclay/Densityclay
          ch2sif6 = v[1]*Molwth2sif6/Densityh2sif6*DensityQuartz/MolwtQuartz+v[1]*Molwth2sif6/Densityh2sif6*4*Densityclay/Molwtclay+v[1]*Molwth2sif6/Densityh2sif6*19*Densitykfelspar/Molwtkfelspar
          csilica = ch2sif6*Molwthsilica/Densitysilica*Densityh2sif6/Molwth2sif6*5
          print(csilica)
          csilicagel = v[5]*Densityalf2/Molwthalf2*3*Molwtsilicagel/Densitysilicagel
          calf2 = v[2]*Molwthalf2/Densityalf2*Densitykfelspar/Molwtkfelspar+v[3]*Molwthalf2/Densityalf2*3*Densityclay/Molwtclay
          
          rh2sif6clay = 7.82*np.power(10,22)*np.exp(-37.8*4184/(8.314*T))*np.power(ch2sif6,0.6)*(v[3]-ch2sif6*Densityh2sif6/Molwth2sif6*1*Molwtclay/Densityclay)*np.power(v[3]*Densityclay/(Densityclay*v[3]+csilica*Densitysilica),0.1)
          rh2sif6felspar = 122775*np.exp(-14*4184/(8.314*T))*np.power(ch2sif6,0.6)*(v[2]-ch2sif6*Densityh2sif6/Molwth2sif6*6*Molwtkfelspar/Densitykfelspar)*np.power(v[2]*Densitykfelspar/(Densitykfelspar*v[2]+csilica*Densitysilica),0.1)
          ralf2 = 5.525*np.exp(-8.12*4184/(8.314*T))*calf2*np.power(rhcl,1.7)*(v[3]-ch2sif6*Densityh2sif6/Molwth2sif6*6*Molwtkfelspar/Densitykfelspar-v[5]*Densityalf2/Molwthalf2*1*Molwtkfelspar/Densitykfelspar)*np.power(v[2]*Densitykfelspar/(Densitykfelspar*v[2]+csilica*Densitysilica),0.1)
          
          
          #return [rcalcite, rquartz]
          return [rcalcite,rquartz,rkfelspar,rclay,rh2sif6felspar,ralf2,rh2sif6clay,rcalcitehf]
      
    #res = solve_ivp(rhs, (0, time), [RockVolume*CalciteRockVolume, RockVolume*QuartzRockVolume,RockVolume*KfelsparRockVolume,RockVolume*ClayRockVolume,0,0],t_eval=np.linspace(0, time,10))
   
    res = solve_ivp(rhs, (0, time), [RockVolume*CalciteRockVolume, RockVolume*QuartzRockVolume,RockVolume*KfelsparRockVolume,RockVolume*ClayRockVolume,0.,0.,0.,0.],t_eval=np.linspace(0, time,10))
    #print("res=",res.y[0])
    csilica = 0.
    csilicagel = 0.
    ch2sif6 = res.y[4]*Molwth2sif6/Densityh2sif6*DensityQuartz/MolwtQuartz+res.y[4]*Molwth2sif6/Densityh2sif6*4*Densityclay/Molwtclay+res.y[4]*Molwth2sif6/Densityh2sif6*19*Densitykfelspar/Molwtkfelspar
    csilica = -ch2sif6*Molwthsilica/Densitysilica*Densityh2sif6/Molwth2sif6*5
    print("ch2sif6=",csilica)
    csilicagel = -res.y[5]*Densityalf2/Molwthalf2*3*Molwtsilicagel/Densitysilicagel
    ccaf2 = res.y[7]*Molwtcaf2/Densitycaf2*DensityCalcite/MolwtCalcite
    print("ccaf2=",ccaf2)
    print(res.y[0])
    df = pd.DataFrame({'Calcite RockVolume': res.y[0], 'time (hrs)': list(np.linspace(0, time*24,10))},columns=["Calcite RockVolume","time (hrs)"])
    st.line_chart(df,x="time (hrs)",y="Calcite RockVolume")
    df = pd.DataFrame({'Quartz RockVolume': res.y[1], 'time (hrs)': list(np.linspace(0, time*24,10))},columns=["Quartz RockVolume","time (hrs)"])
    st.line_chart(df,x="time (hrs)",y="Quartz RockVolume")
    df = pd.DataFrame({'Concentration Silica': csilica, 'time (hrs)': list(np.linspace(0, time*24,10))},columns=["Concentration Silica","time (hrs)"])
    st.line_chart(df,x="time (hrs)",y="Concentration Silica")
    df = pd.DataFrame({'Silicagel': csilicagel, 'time (hrs)': list(np.linspace(0, time*24,10))},columns=["Silicagel","time (hrs)"])
    st.line_chart(df,x="time (hrs)",y="Silicagel")
    df = pd.DataFrame({'K-Felspar RockVolume': res.y[2], 'time (hrs)': list(np.linspace(0, time*24,10))},columns=["K-Felspar RockVolume","time (hrs)"])
    st.line_chart(df,x="time (hrs)",y="K-Felspar RockVolume")
    df = pd.DataFrame({'Clay RockVolume': res.y[3], 'time (hrs)': list(np.linspace(0, time*24,10))},columns=["Clay RockVolume","time (hrs)"])
    st.line_chart(df,x="time (hrs)",y="Clay RockVolume")
    df = pd.DataFrame({'Porosity': (por+(-res.y[0]-res.y[1]+RockVolume*CalciteRockVolume-res.y[2]+RockVolume*KfelsparRockVolume+RockVolume*QuartzRockVolume-res.y[3]+RockVolume*ClayRockVolume-csilica-csilicagel-ccaf2)/RockVolume), 'time (hrs)': list(np.linspace(0, time*24,10))},columns=["Porosity","time (hrs)"])
# #print(np.vstack(((por-res.y[2]*aspmolwt/aspdens), np.linspace(0, time,10))))
    st.line_chart(df,x="time (hrs)",y="Porosity")
    df = pd.DataFrame({'Permeability': (kperm*np.power((por+(-res.y[0]-res.y[1]-res.y[2]-res.y[3]-csilica-csilicagel-ccaf2+RockVolume*CalciteRockVolume+RockVolume*KfelsparRockVolume+RockVolume*QuartzRockVolume+RockVolume*ClayRockVolume))/(por),g)), 'time (hrs)': list(np.linspace(0, time*24,10))},columns=["Permeability","time (hrs)"])
    st.line_chart(df,x="time (hrs)",y="Permeability")
    # st.line_chart(res.y[0])
    # st.line_chart(res.y[1])
    # st.line_chart(res.y[2])
    # print(res.y[2])
    # print(por-res.y[2]*aspmolwt/aspdens)
    # print(kperm*np.power((por-res.y[2]*aspmolwt/aspdens)/(por-asppor),g))
    # print(res.t)
    
    #plt.plot(res.t, res.y.T)