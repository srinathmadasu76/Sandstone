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



st.header('Sandstone Acidizing Simulation Expert System- Linqx')   
      
por = st.number_input('Enter initial porosity: ', value=0.2)
kperm = st.number_input('Enter inital permeability mD: ', value=100)

Temperature = st.number_input('Enter Temperature F: ', value=100)


RockVolume = st.number_input('Enter rock volume cm3: ', value=1000.)

CalciteRockVolume = st.number_input('Enter calcite rock fraction volume : ', value=0.3)
QuartzRockVolume = st.number_input('Enter quatrz rock fraction volume : ', value=0.0002)
KfelsparRockVolume = st.number_input('Enter k-felspar rock fraction volume : ', value=0.0002)
ClayRockVolume = st.number_input('Enter clay rock fraction volume : ', value=0.1)

DolomiteRockVolume = st.number_input('Enter dolomite rock fraction volume : ', value=0.01)
SideriteRockVolume = st.number_input('Enter siderite rock fraction volume : ', value=0.01)
AnkeriteRockVolume = st.number_input('Enter ankerite rock fraction volume : ', value=0.01)
ZeoliteRockVolume = st.number_input('Enter zeolite rock fraction volume : ', value=0.01)
NafelsparRockVolume = st.number_input('Enter Na-felspar rock fraction volume : ', value=0.01)
KaoliniteRockVolume = st.number_input('Enter kaolinite rock fraction volume : ', value=0.01)
ChloriteRockVolume = st.number_input('Enter chlorite rock fraction volume : ', value=0.005)
MuscoviteRockVolume = st.number_input('Enter muscovite rock fraction volume : ', value=0.01)
MontmorilloniteRockVolume = st.number_input('Enter Montmorillonite rock fraction volume : ', value=0.005)

Nonreactiverockvolume = RockVolume*(1-(CalciteRockVolume+QuartzRockVolume+KfelsparRockVolume+ClayRockVolume+DolomiteRockVolume+SideriteRockVolume+AnkeriteRockVolume+ZeoliteRockVolume+NafelsparRockVolume+KaoliniteRockVolume+ChloriteRockVolume+MuscoviteRockVolume+MontmorilloniteRockVolume))
RockVolumereactive = RockVolume-Nonreactiverockvolume
ClayExchangeMineralVolume = ClayRockVolume+SideriteRockVolume+AnkeriteRockVolume+ZeoliteRockVolume+NafelsparRockVolume+KaoliniteRockVolume+ChloriteRockVolume+MuscoviteRockVolume+MontmorilloniteRockVolume

T = 273+(Temperature-32)*5/9 #K
Iron_Presence = st.checkbox('Iron-bearing formations / pyrite presence')
Organic_Presence = st.checkbox('Organic blockage (oil, asphaltenes, polymers)')
Mechanical_Presence = st.checkbox('Mechanical fines / compacted solids')



if st.button('Advice Expert System'):
    if Iron_Presence:
        st.write('Preflush Stge:HCl/HF with strong iron chelants (EDTA, citric)')
    if Organic_Presence:
        st.write('Preflush Stge:Solvent or mutual solvent preflush')
    if Mechanical_Presence:
        st.write('Preflush Stge:Preflush with chelants and solvents; acid to dissolve cements')
    if(ClayRockVolume>0. or ZeoliteRockVolume>0.):
        st.text_area("Preflush Stge:", "Condition with CLAY-SAFE 5")
        if(ClayExchangeMineralVolume>0.15):
            if(MontmorilloniteRockVolume>0.15):
                st.text_area("Preflush Stge:", "Condition with CLAYSAFE 5 ahead of CLAY-SAFE 5")
            if (CalciteRockVolume > 0.05):
                     st.text_area("Preflush Stge:", "Conditioning Acid Volume = HF Treatment Volume")
            else:
                    st.text_area("Preflush Stge:", "Conditioning Acid Volume = 75% HF Treatment Volume")     
    else:
        st.text_area("Preflush Stge:", "Condition with HCL")
        if(CalciteRockVolume > 0.05):
                 st.text_area("Preflush Stge:", "Conditioning Acid Volume = HF Treatment Volume")
        else:
                 st.text_area("Preflush Stge:", "Conditioning Acid Volume = 75% HF Treatment Volume")
                 
    if(KfelsparRockVolume>0. or NafelsparRockVolume>0.):
        if(KfelsparRockVolume>0.2):
            if(T>175):
                st.text_area("Main Stge:", "Main Stage: Sandstone Completion Acid, Fines Control Acid")
            else:
                 st.text_area("Main Stge:", "Main Stage: K-spar, Fines Control Acid")
        elif(NafelsparRockVolume>0.2):
            if(T>250):
                st.text_area("Main Stge:", "Main Stage: Sandstone Completion Acid, Fines Control Acid")
            else:
                 st.text_area("Main Stge:", "Main Stage: K-spar, Fines Control Acid")
    
        if(ChloriteRockVolume>0. or ZeoliteRockVolume>0.):
            st.text_area("Main Stge:", "Treat with Volcanic Acid, Choose most compatable or weakest HF blend from all applicable categories")
        elif(ClayRockVolume>0.25):
                 if(T>200):
                     st.text_area("Main Stge:", "Main Stage: Sandstone Completion, Fines Control Acid")
                 else:
                    st.text_area("Main Stge:", "Main Stage: K-spar, Fines Control Acid")
            
        elif(ZeoliteRockVolume>0.2):
            st.text_area("Main Stge:", "Treat with Volcanic Acid, Choose most compatable or weakest HF blend from all applicable categories")
        else:
            st.text_area("Main Stge:", "Treat with Sandstone Completion Acid, Choose most compatable or weakest HF blend from all applicable categories")
    
    if(MontmorilloniteRockVolume>0.2):
                st.text_area("Main Stge:", "Treat with HF+Acetic Acid/Formic Acid")
    st.text_area("Overflush Stge:", "Treat with KCL brine or NH4CL brine")
    st.text_area("Injection rate: Below fracture pressure for matrix acidizing.")
    st.text_area("Contact time: 30–60 min typical for HF systems.")
    st.text_area("Temperature control: Monitor to avoid accelerated HF reactivity.")
    st.text_area("Volume guidelines:Preflush: 30–50 gal/ft,Main acid: 50–100 gal/ftOverflush: 20–40 gal/ft")
