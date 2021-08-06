# -*- coding: utf-8 -*-
"""
Created on Mon Aug  2 13:31:30 2021

@author: Jacko
"""
import pandas as pd
from math import e
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
    y=float(e**((22.434/((273.15+T)**2))-0.0002524))
    print(y)

acid_frac_alpha()
