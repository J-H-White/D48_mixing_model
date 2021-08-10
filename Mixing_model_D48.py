# -*- coding: utf-8 -*-
"""
Created on Mon Aug  2 13:31:30 2021

@author: Jacko
"""
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

df= pd.read_excel(r"C:\Users\Jacko\Desktop\Modelling\Clumped Mixing Line feat D48.xlsx"\
                  ,sheet_name="D48", index_col=False)


    
components_data = df.iloc[0:5, 0:3]
transfer = df.iloc[9:13, 0:3]
transfer_function=transfer.rename(columns={"Component 1":"Transfer Function"\
                                           ,"Component 2":"Heated Gas Line"})


def d18O_alpha(T=25):
    y = float(1.00397+((5.25*10**2)/(273.15+T)**2))
    print(y)

def acid_frac_alpha(T=25):
    y = float(e**((22.434/((273.15+T)**2))-0.0002524))
    return(y)

# Eq. 1
def R_13(d13C):
    y = float(((d13C/1000)+1)*0.0112372)
    return(y)

# Eq. 2
def R18(d18O_VSMOW):
    y = float(((d18O_VSMOW/1000)+1)*0.0020052)
    return(y)

# Eq. 3
def R17(R18):
    y = float(((R18/0.0020052)**0.5164)*0.0003799)
    return(y)

# Eq. 4
def 12C(R13):
    y = float(1/(1+R13))
    return y

# Eq. 5
def 13C(R13):
    y = 12C(R13)*R13
    return y

# Eq. 6
def 16O(R17, R18):
    y = 1/(1+R17+R18)
    return y

# Eq. 7
def 17O(16O):
    y = 16O*R17()
    return y

# Eq. 8
def 18O(16O):
    y = 16O*R18()
    return y

# Eq. 9
def R45*(12C, 13C, 16O, 17O):
    y = (13C*16O**2+2*12C*16O*17O)/(12C*16O**2)
    return y

# Eq. 10
def R46*(12C, 13C, 16O, 17O, 18O):
    y = (2*12C*18O+12C*17O**2+2*13C*16O*17O)/(12C*16O**2)
    return y

# Eq. 11
def R47*(12C, 13C, 16O, 17O, 18O):
    y = (2*13C*18O+12C*17O**2+2*13C*17O*18O)/(12C*16O**2)
    return y

# Eq. 12
def d45 (R45*, R45*_WG):
    y = ((R45*/R45*_WG)-1)*1000
    return y

# Eq. 13
def d46 (R46*, R46*_WG):
    y = ((R46*/R46*_WG)-1)*1000
    return y

# Eq. 14
def D47_RF (D47EndMemberValue, acid_frac_alpha):
    y = D47EndMemberValue - acid_frac_alpha
    return y

# Eq. 15
def D47_SGvsWG_0 (D47_RF, ETF_int, ETF_slope):
    y = (D47_RF - ETF_int)/ETF_slope
    return y

# Eq. 16
def d47 (D47_SGvsWG_0, R47, R47*, R47*_WG, EGL_slope):
    y = ((D47_SGvsWG_0() - 1000)*R47() - 1000*R47*_WG)/(R47*_WG - EGL_slope*R47*)
    return y

# Eq. 17
def dmix (x, d):
    y = sum(x*d)
    return y

# Eq. 18
def R (d, R_WG):
    y = ((d/1000) + 1)*R_WG
    return y

# Eq. 19
def D47_SGvsWG_mix (R45, R45*, R46, R46*, R47, R47*):
    y = ((R47/R47* - 1) - (R46/R46* - 1) - (R45/R45* - 1))*1000
    return y

# Eq. 20
def D47_SGvsWG_0_mix (D47_SGvsWG_mix, d47_mix, EGL_slope):
    y = D47_SGvsWG_mix - d47_mix*EGL_slope
    return y 

# Eq. 21
def D47mix (D47_SGvsWG_0_mix, ETF_slope, ETF_int):
    y = D47_SGvsWG_0_mix*ETF_slope + ETF_int
    return y

