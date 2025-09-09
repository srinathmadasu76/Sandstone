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

CalciteRockVolume = st.number_input('Enter calcite rock fraction volume : ', value=0.1)
QuartzRockVolume = st.number_input('Enter quatrz rock fraction volume : ', value=0.2)
KfelsparRockVolume = st.number_input('Enter k-felspar rock fraction volume : ', value=0.4)
ClayRockVolume = st.number_input('Enter clay rock fraction volume : ', value=0.3)

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
rhcl = conchcl 
mult = 1
rquartz = SurfaceareaQuartz*kcoQuartz*np.exp(-EaQuartz*1000/(R*T))
hftest = (conchf-RockVolume*QuartzRockVolume*Molwthf/Densityhf*6*DensityQuartz/MolwtQuartz-RockVolume*QuartzRockVolume*Molwthf/Densityhf*36*Densityclay/Molwtclay-RockVolume*QuartzRockVolume*Molwthf/Densityhf*19*Densitykfelspar/Molwtkfelspar)*MolwtQuartz/DensityQuartz
#print(rquartz)
#print("hf=",hftest)
rsilica = 0.
if st.button('Calculate'):
    
    def rhs(s, v): 
              
          #rcalcite = -SurfaceareaCalcite*kcoCalcite*np.exp(-EaCalcite*1000/(R*T))*np.power(v[2],0.63)*MolwtCalcite/DensityCalcite
          rhcl = -SurfaceareaCalcite*kcoCalcite*np.exp(-EaCalcite*1000/(R*T))*np.power(v[0],0.63)*2-Surfaceareakfelspar*kcokfelspar*np.exp(-Eakfelspar*1000/(R*T))*(1+ 0.0566*np.exp(956/T))*(1+np.power(v[0],0.4))*np.power(v[1],1.2)*3-18*122775*np.exp(-14*4184/(8.314*T))*np.power(v[6],0.6)*(v[4])*np.power(v[4]*Molwtkfelspar/(Molwtkfelspar*v[4]+v[7]*Molwtsilica),0.1)-6*7.82*np.power(10,22)*np.exp(-37.8*4184/(8.314*T))*np.power(v[6],0.6)*(v[5])*np.power(v[5]*Molwtclay/(Molwtclay*v[3]+v[7]*Molwtsilica),0.1)-4*5.525*np.exp(-8.12*4184/(8.314*T))*v[8]*np.power(v[0],1.7)*(v[4])*np.power(v[4]*Molwtkfelspar/(Molwtkfelspar*v[4]+v[7]*Molwtsilica),0.1)
          rhf = -SurfaceareaCalcite*0.001*np.exp(-40*1000/(8.314*T))*v[1]*2-SurfaceareaQuartz*kcoQuartz*np.exp(-EaQuartz*1000/(R*T))*v[1]*6 -Surfaceareakfelspar*kcokfelspar*np.exp(-Eakfelspar*1000/(R*T))*(1+ 0.0566*np.exp(956/T))*(1+np.power(v[0],0.4))*np.power(v[1],1.2)*19-Surfaceareaclay*kcoclay*np.exp(-Eaclay*1000/(R*T))*v[1]*36
          rcalcite = -SurfaceareaCalcite*0.001*np.exp(-40*1000/(8.314*T))*v[1]-SurfaceareaCalcite*kcoCalcite*np.exp(-EaCalcite*1000/(R*T))*np.power(v[0],0.63)
          rquartz = -SurfaceareaQuartz*kcoQuartz*np.exp(-EaQuartz*1000/(R*T))*v[1]
          rkfelspar = -Surfaceareakfelspar*kcokfelspar*np.exp(-Eakfelspar*1000/(R*T))*(1+ 0.0566*np.exp(956/T))*(1+np.power(v[0],0.4))*np.power(v[1],1.2)-6*122775*np.exp(-14*4184/(8.314*T))*np.power(v[6],0.6)*(v[4])*np.power(v[4]*Molwtkfelspar/(Molwtkfelspar*v[4]+v[7]*Molwtsilica),0.1)-5.525*np.exp(-8.12*4184/(8.314*T))*v[8]*np.power(v[0],1.7)*(v[4])*np.power(v[4]*Molwtkfelspar/(Molwtkfelspar*v[4]+v[7]*Molwtsilica),0.1)
          rclay = -Surfaceareaclay*kcoclay*np.exp(-Eaclay*1000/(R*T))*v[1]+7.82*np.power(10,22)*np.exp(-37.8*4184/(8.314*T))*np.power(v[6],0.6)*(v[5])*np.power(v[5]*Molwtclay/(Molwtclay*v[3]+v[7]*Molwtsilica),0.1)
          rh2sif6felspar = -6*122775*np.exp(-14*4184/(8.314*T))*np.power(v[6],0.6)*(v[4])*np.power(v[4]*Molwtkfelspar/(Molwtkfelspar*v[4]+v[7]*Molwtsilica),0.1)
          rsilica = 122775*np.exp(-14*4184/(8.314*T))*np.power(v[6],0.6)*(v[4])*np.power(v[4]*Molwtkfelspar/(Molwtkfelspar*v[4]+v[7]*Molwtsilica),0.1)+4*7.82*np.power(10,22)*np.exp(-37.8*4184/(8.314*T))*np.power(v[6],0.6)*(v[5])*np.power(v[5]*Molwtclay/(Molwtclay*v[5]+v[7]*Molwtsilica),0.1)
          rh2sif6clay = -7.82*np.power(10,22)*np.exp(-37.8*4184/(8.314*T))*np.power(v[6],0.6)*(v[5])*np.power(v[5]*Molwtclay/(Molwtclay*v[5]+v[7]*Molwtsilica),0.1)
          rh2sif6 = -rquartz-4*rclay-3*rkfelspar+rh2sif6felspar/6+rh2sif6clay
          ralf2 = 6*122775*np.exp(-14*4184/(8.314*T))*np.power(v[6],0.6)*(v[4])*np.power(v[4]*Molwtkfelspar/(Molwtkfelspar*v[4]+v[7]*Molwtsilica),0.1)+Surfaceareakfelspar*kcokfelspar*np.exp(-Eakfelspar*1000/(R*T))*(1+ 0.0566*np.exp(956/T))*(1+np.power(v[0],0.4))*np.power(v[1],1.2)-5.525*np.exp(-8.12*4184/(8.314*T))*v[8]*np.power(v[0],1.7)*(v[4])*np.power(v[4]*Molwtkfelspar/(Molwtkfelspar*v[4]+v[7]*Molwtsilica),0.1)
          rsilicagel = 3*5.525*np.exp(-8.12*4184/(8.314*T))*v[8]*np.power(v[0],1.7)*(v[4])*np.power(v[4]*Molwtkfelspar/(Molwtkfelspar*v[4]+v[7]*Molwtsilica),0.1)
          rcaf2 = SurfaceareaCalcite*0.001*np.exp(-40*1000/(8.314*T))*v[1]
          # rquartz = -SurfaceareaQuartz*kcoQuartz*np.exp(-EaQuartz*1000/(R*T))*(conchf -v[1]*Molwthf/Densityhf*6*DensityQuartz/MolwtQuartz-v[1]*Molwthf/Densityhf*36*Densityclay/Molwtclay-v[1]*Molwthf/Densityhf*19*Densitykfelspar/Molwtkfelspar)*MolwtQuartz/DensityQuartz
          # #print(rquartz)
          
         
          # rkfelspar = -Surfaceareakfelspar*kcokfelspar*np.exp(-Eakfelspar*1000/(R*T))*(1+ 0.0566*np.exp(956/T))*(1+np.power(rhcl,0.4))*np.power(rhf,1.2)*Molwtkfelspar/Densitykfelspar
          # rclay = -Surfaceareaclay*kcoclay*np.exp(-Eaclay*1000/(R*T))*rhf*Molwtclay/Densityclay
          # ch2sif6 = v[1]*Molwth2sif6/Densityh2sif6*DensityQuartz/MolwtQuartz+v[1]*Molwth2sif6/Densityh2sif6*4*Densityclay/Molwtclay+v[1]*Molwth2sif6/Densityh2sif6*19*Densitykfelspar/Molwtkfelspar
          # csilica = ch2sif6*Molwthsilica/Densitysilica*Densityh2sif6/Molwth2sif6*5
          # print(csilica)
          # csilicagel = v[5]*Densityalf2/Molwthalf2*3*Molwtsilicagel/Densitysilicagel
          # calf2 = v[2]*Molwthalf2/Densityalf2*Densitykfelspar/Molwtkfelspar+v[3]*Molwthalf2/Densityalf2*3*Densityclay/Molwtclay
          
          # rh2sif6clay = 7.82*np.power(10,22)*np.exp(-37.8*4184/(8.314*T))*np.power(ch2sif6,0.6)*(v[3]-ch2sif6*Densityh2sif6/Molwth2sif6*1*Molwtclay/Densityclay)*np.power(v[3]*Densityclay/(Densityclay*v[3]+csilica*Densitysilica),0.1)
          # rh2sif6felspar = 122775*np.exp(-14*4184/(8.314*T))*np.power(ch2sif6,0.6)*(v[2]-ch2sif6*Densityh2sif6/Molwth2sif6*6*Molwtkfelspar/Densitykfelspar)*np.power(v[2]*Densitykfelspar/(Densitykfelspar*v[2]+csilica*Densitysilica),0.1)
          # ralf2 = 5.525*np.exp(-8.12*4184/(8.314*T))*calf2*np.power(rhcl,1.7)*(v[3]-ch2sif6*Densityh2sif6/Molwth2sif6*6*Molwtkfelspar/Densitykfelspar-v[5]*Densityalf2/Molwthalf2*1*Molwtkfelspar/Densitykfelspar)*np.power(v[2]*Densitykfelspar/(Densitykfelspar*v[2]+csilica*Densitysilica),0.1)
          
          
          return [rhcl,rhf,rcalcite,rquartz,rkfelspar,rclay,rh2sif6,rsilica,ralf2,rsilicagel,rcaf2]
          #return [rcalcite,rquartz,rkfelspar,rclay,rh2sif6felspar,ralf2,rh2sif6clay,rcalcitehf]
      
    res = solve_ivp(rhs, (0, time), [conchcl, conchf, RockVolume*CalciteRockVolume/MolwtCalcite*DensityCalcite, RockVolume*QuartzRockVolume/MolwtQuartz*DensityQuartz,RockVolume*KfelsparRockVolume*Densitykfelspar/Molwtkfelspar,RockVolume*ClayRockVolume/Molwtclay*Densityclay,0.,0.,0.,0.,0.],t_eval=np.linspace(0, time,10))
   
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
    df = pd.DataFrame({'Porosity': (por+(-res.y[2]*MolwtCalcite/DensityCalcite+RockVolume*CalciteRockVolume-res.y[3]*MolwtQuartz/DensityQuartz+RockVolume*KfelsparRockVolume+RockVolume*QuartzRockVolume-res.y[3]*MolwtQuartz/DensityQuartz+RockVolume*ClayRockVolume-res.y[4]/Densitykfelspar*Molwtkfelspar-res.y[7]*Molwtsilica/Densitysilica-res.y[9]*Molwtsilicagel/Densitysilicagel-res.y[10]*Molwtcaf2/Densitycaf2)/RockVolume), 'time (seconds)': list(np.linspace(0, time,10))},columns=["Porosity","time (seconds)"])
# # #print(np.vstack(((por-res.y[2]*aspmolwt/aspdens), np.linspace(0, time,10))))
    st.line_chart(df,x="time (seconds)",y="Porosity")
    df = pd.DataFrame({'Permeability': (kperm*np.power((por+(-res.y[2]*MolwtCalcite/DensityCalcite+RockVolume*CalciteRockVolume-res.y[3]*MolwtQuartz/DensityQuartz+RockVolume*KfelsparRockVolume+RockVolume*QuartzRockVolume-res.y[3]*MolwtQuartz/DensityQuartz+RockVolume*ClayRockVolume-res.y[4]/Densitykfelspar*Molwtkfelspar-res.y[7]*Molwtsilica/Densitysilica-res.y[9]*Molwtsilicagel/Densitysilicagel-res.y[10]*Molwtcaf2/Densitycaf2)/RockVolume)/(por),g)), 'time (seconds)': list(np.linspace(0, time,10))},columns=["Permeability","time (seconds)"])
    st.line_chart(df,x="time (seconds)",y="Permeability")
    # st.line_chart(res.y[0])
    # st.line_chart(res.y[1])
    # st.line_chart(res.y[2])
    # print(res.y[2])
    # print(por-res.y[2]*aspmolwt/aspdens)
    # print(kperm*np.power((por-res.y[2]*aspmolwt/aspdens)/(por-asppor),g))
    # print(res.t)
    
    #plt.plot(res.t, res.y.T)