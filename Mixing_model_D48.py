# -*- coding: utf-8 -*-
"""
Created on Mon Aug  2 13:31:30 2021

@author: Jacko
"""
#Import packages that are not Built-in functions

import pandas as pd
from math import e

#Create an object for the tabular dataset.
#In the case of an Excel spreadsheet file, we are interested in a particular
#sheet (D48) and we are not interested in the indexing on the left.

df= pd.read_excel(r"C:\Users\Jacko\Desktop\Modelling\Clumped Mixing Line feat D48.xlsx"
                  ,sheet_name="D48", index_col=False)

# convert d18O (VPDB) and d13C (VPDB) data to absolute ratios (Ri)

df_d18O = df.iloc[3, 0:3]
df_d13C = df.iloc[4, 0:3]

# Eq. 1
def R13(d13C_value):
    R13_value = ((d13C_value/1000)+1)*0.01118
    return(R13_value)

R13_list = []

for i in range(1, len(df_d13C)):
    R13_value = R13(df_d13C[i])
    R13_list.append(R13_value)

R13_tuple = tuple(R13_list)

# Eq. 2
def R18(d18O_value):
    R18_value = ((d18O_value/1000)+1)*0.0020052
    return(R18_value)

R18_list = []

for i in range(1, len(df_d18O)):
    R18O_value = R18(df_d18O[i])
    R18_list.append(R18O_value)

R18_tuple = tuple(R18_list)

# Eq. 3
def R17(R18_value):
    R17_value = ((R18_value/0.0020052)**0.5164)*0.0003799
    return(R17_value)

R17_list = []

for i in range(len(R18_list)):
    R17_value = R17(R18_list[i])
    R17_list.append(R17_value)

R17_tuple = tuple(R17_list)

# Eq. 4
def C12(R13_value):
    C12_value = 1/(1+R13_value)
    return C12_value

C12_list = []

for i in range(len(R13_list)):
    C12_value = C12(R13_list[i])
    C12_list.append(C12_value)

C12_tuple = tuple(C12_list)

# Eq. 5
def C13(R13_value):
    C13_value = C12(R13_value)*R13_value
    return C13_value

C13_list = []

for i in range(len(R13_list)):
    C13_value = C13(R13_list[i])
    C13_list.append(C13_value)
    
C13_tuple = tuple(C13_list)

# Eq. 6
def O16(R17_value, R18_value):
    O16_value = 1/(1+R17_value+R18_value)
    return O16_value

O16_list = []

for i in range(min([len(R17_tuple), len(R18_tuple)])):
    O16_value = O16(R17_tuple[i], R18_tuple[i])
    O16_list.append(O16_value)

O16_tuple = tuple(O16_list)

O17_list = []

for i in range(min([len(O16_tuple), len(R17_tuple)])):
    O17_value = O16_tuple[i]*R17_tuple[i]
    O17_list.append(O17_value)

O17_tuple = tuple(O17_list)

O18_list = []

for i in range(min([len(O16_tuple), len(R18_tuple)])):
    O18_value = O16_tuple[i]*R18_tuple[i]
    O18_list.append(O18_value)

O18_tuple = tuple(O18_list)

R45_sto_list = []

for i in range(min([len(C12_tuple), len(C13_tuple), len(O16_tuple), len(O17_tuple)])):
    R45_sto_value = (C13_tuple[i]*O16_tuple[i]**2+2*C12_tuple[i]*O16_tuple[i]*O17_tuple[i])/(C12_tuple[i]*O16_tuple[i]**2)
    R45_sto_list.append(R45_sto_value)

R45_sto_tuple = tuple(R45_sto_list)

R46_sto_list = []

for i in range(min([len(C12_tuple), len(C13_tuple), len(O16_tuple), len(O17_tuple), len(O18_tuple)])):
    R46_sto_value = (2*C12_tuple[i]*O16_tuple[i]*O18_tuple[i] + C12_tuple[i]*O17_tuple[i]**2 + 2*C13_tuple[i]*O16_tuple[i]*O17_tuple[i])/(C12_tuple[i]*O16_tuple[i]**2)
    R46_sto_list.append(R46_sto_value)

R46_sto_tuple = tuple(R46_sto_list)

R47_sto_list = []

for i in range(min([len(C12_tuple), len(C13_tuple), len(O16_tuple), len(O17_tuple), len(O18_tuple)])):
    R47_sto_value = (2*C13_tuple[i]*O16_tuple[i]*O18_tuple[i] + C13_tuple[i]*O17_tuple[i]**2 + 2*C12_tuple[i]*O17_tuple[i]*O18_tuple[i])/(C12_tuple[i]*O16_tuple[i]**2)
    R47_sto_list.append(R47_sto_value)

R47_sto_tuple = tuple(R47_sto_list)

# Eq. 12
#def d45 (R45*, R45*_WG):
#    y = ((R45*/R45*_WG)-1)*1000
#    return y

# Eq. 13
#def d46 (R46*, R46*_WG):
#    y = ((R46*/R46*_WG)-1)*1000
#    return y

#transfer = df.iloc[9:13, 0:3]
#transfer_function=transfer.rename(columns={"Component 1":"Transfer Function"                                         ,"Component 2":"Heated Gas Line"})

#def d18O_alpha(T=25):
#    y = float(1.00397+((5.25*10**2)/(273.15+T)**2))
#    print(y)

#def acid_frac_alpha(T=25):
#    y = float(e**((22.434/((273.15+T)**2))-0.0002524))
#    return(y)
    
# Eq. 14
#def D47_RF (D47EndMemberValue, acid_frac_alpha):
#    y = D47EndMemberValue - acid_frac_alpha
#    return y

# Eq. 15
#def D47_SGvsWG_0 (D47_RF, ETF_int, ETF_slope):
#    y = (D47_RF - ETF_int)/ETF_slope
#    return y

# Eq. 16
#def d47 (D47_SGvsWG_0, R47, R47*, R47*_WG, EGL_slope):
#    y = ((D47_SGvsWG_0() - 1000)*R47() - 1000*R47*_WG)/#(R47*_WG - EGL_slope*R47*)
#    return y

# Eq. 17
#def dmix (x, d):
#    y = sum(x*d)
#    return y

# Eq. 18
#def R (d, R_WG):
#    y = ((d/1000) + 1)*R_WG
#    return y

# Eq. 19
#def D47_SGvsWG_mix (R45, R45*, R46, R46*, R47, R47*):
#    y = ((R47/R47* - 1) - (R46/R46* - 1) - (R45/R45* - 1))*1000
#    return y

# Eq. 20
#def D47_SGvsWG_0_mix (D47_SGvsWG_mix, d47_mix, EGL_slope):
#    y = D47_SGvsWG_mix - d47_mix*EGL_slope
#    return y 

# Eq. 21
#def D47mix (D47_SGvsWG_0_mix, ETF_slope, ETF_int):
#    y = D47_SGvsWG_0_mix*ETF_slope + ETF_int
#    return y
