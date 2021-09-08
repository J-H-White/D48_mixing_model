 # -*- coding: utf-8 -*-
"""
Created on Mon Aug  2 13:31:30 2021

@author: Jacko
This model describes the non-linear mixing effects on mass-47 CO2 clumped isotope data. This code is written based on Defliese & Lohmann (2015)'s work. Here, we try to extend this model to the non-linear mixing effects on mass-48 CO2 clumped isotope data. To make reading easier, we recommend wrapping the code lines in Spyder by going to Preferences -> Editor -> select the Display tab -> tick wrap lines.

References

Defliese, W.F., Lohmann, K.C., 2015. Non-linear mixing effects on mass-47 CO2clumped isotope thermometry: Patterns and implications. Rapid Communications in Mass Spectrometry 29, 901–909.. doi:10.1002/rcm.7175
"""

#Import packages that are not Built-in functions

from math import e

""" In this section, the input the values for the isotope geochemistry data for the two components (component_1_D47, component_1_D48, component_1_d18O_VPDB, component_1_d13C_VPDB, component_1_contribution, component_2_D47, component_2_D48, component_2_d18O_VPDB, component_2_d13C_VPDB); transfer function slopes and intercepts (transfer_function_slope_D47, transfer_function_intercept_D47, transfer_function_slope_D48, transfer_function_intercept_D48); heated gas line slopes (heated_gas_line_slope_D47 and heated_gas_line_slope_D48); working gas delta values (working_gas_d13C_VPDB and working_gas_d18O_VSMOW); and the acid reaction temperature (T). """  

component_1_D47 = 0.214
component_1_D48 = 0.138
component_1_d18O_VPDB = -2.19
component_1_d13C_VPDB = 2.02
component_1_contribution = 0

component_2_D47 = 0.215
component_2_D48 = 0.138
component_2_d18O_VPDB = -18.69
component_2_d13C_VPDB = -10.17

# Note there is no component 2 contribution input, this is because this model is a mixture between two components. Therefore, the component 2 contribution is the compliment of the component 1 contribution (1 - component_1_contribution). Contributions have values between 0 and 1.

T = 25

transfer_function_slope_D47 = 1.0391
transfer_function_intercept_D47 = 0.9499
transfer_function_slope_D48 = 1.0005
transfer_function_intercept_D48 = 0.3367

heated_gas_line_slope_D47 = 0.0274
heated_gas_line_slope_D48 = 0

reference_gas_d13C_VPDB = -3.700
reference_gas_d18O_VSMOW = 34.990

""" d18O_alpha_value is defnied as the oxygen fractionation from the conversion from CaCO3 to CO2 at a temperature T. """
d18O_alpha = 1.00397+((5.25*10**2)/(273.15+T)**2)

""" acid_frac_alpha is fractionation due to acid digestion. """
acid_frac_alpha = e**((22.434/((273.15 + T)**2))-0.0002524)

acid_frac_offset = (1000/acid_frac_alpha)-1000

def mixture (C1, C2, contribution_1):
    """ mixture defines the weighted average given the inputs. """
    mixture_value = (C1*contribution_1) + (C2*(1-contribution_1))
    return mixture_value

# d13C mixture value from component 1 and component 2
mixture_d13C_VPDB = mixture(component_1_d13C_VPDB, component_2_d13C_VPDB, component_1_contribution)

# d18O mixture value from component 1 and component 2
mixture_d18O_VPDB = mixture(component_1_d18O_VPDB, component_2_d18O_VPDB, component_1_contribution)

""" What the hell is going on here Will? """

# R18 from VSMOW scale to VPDB calcite ((30.92/1000 + 1) * 0.0020052) and scaled to VPDB CO2 by multiplying by a fractionation factor of 1.01025
    
#VPDB_calcite_R18 = ((30.92/1000+1)*0.0020052)*1.01025

# R17 for VSMOW to VPDB calcite

#VPDB_calcite_R17 = ((VPDB_calcite_R18/(0.0020052*1.01025)**0.528)*0.00038475

# Convert d18O values from the VPDB scale to VSMOW scale and corrected for CaCO3 to CO2 fractionation (d18O_alpha_value).
def CO3VPDB_to_CO2VSMOW(d18O_VPDB_value):
    CO3_VSMOW_value = ((d18O_VPDB_value + 29.98)/0.97001)
    CO2_VSMOW_value = ((1000+CO3_VSMOW_value)*d18O_alpha)-1000
    return CO2_VSMOW_value

mixture_d18O_CO2VSMOW = CO3VPDB_to_CO2VSMOW(mixture_d18O_VPDB)

component_1_d18O_CO2VSMOW = CO3VPDB_to_CO2VSMOW(component_1_d18O_VPDB)

component_2_d18O_CO2VSMOW = CO3VPDB_to_CO2VSMOW(component_2_d18O_VPDB)

""" Lists have been created to perform calculations on the equations will preform calculations on. Python has different data structures of which lists have the property that they can be changed, which is what our intention is by calculating with these equations below. After these calculations are done, the 'keys' list will be merged with the two 'values' lists to create a dictionary. Dictionaries are data structures which have a key, which can be called, and results in a output. """

keys = ['Reference_gas', 'Composite_gas', 'Component_1', 'Component_2']

values_d13C = [reference_gas_d13C_VPDB, mixture_d13C_VPDB, component_1_d13C_VPDB, component_2_d13C_VPDB]

values_d18O = [reference_gas_d18O_VSMOW, mixture_d18O_CO2VSMOW, component_1_d18O_CO2VSMOW, component_2_d18O_CO2VSMOW]

""" Here we define functions which convert the d13C and d18O data into R13, R18 and R17 data. The equation numbers (Eq.) correspond to the same as in Defliese & Lohmann (2015). """

R13_list = []

for i in range(len(values_d13C)):
    # Eq. 1
    R13_values = ((values_d13C[i]/1000)+1)*0.01118
    R13_list.append(R13_values)

R18_list = []

for i in range(len(values_d18O)):
    # Eq. 2
    R18_values = ((values_d18O[i]/1000)+1)*0.0020052
    R18_list.append(R18_values)

R17_list = []

for i in range(len(R18_list)):
    # Eq. 3
    R17_values = ((R18_list[i]/0.0020052)**0.5164)*0.00038475
    R17_list.append(R17_values)

""" Based on the R13, R18 and R17 data, the concentrations for each isotopes can be estimated. That is 12C, 13C, 16O, 17O, and 18O. The values and letters had to be switched because 12C, for exmaple, can produce later syntax errors when calculations are being preformed. """

C12_list = []

for i in range(len(R13_list)):
    # Eq. 4
    C12_value = 1/(1+R13_list[i])
    C12_list.append(C12_value)

C13_list = []

for i in range(min([len(R13_list), len(C12_list)])):
    # Eq. 5
    C13_value = C12_list[i]*R13_list[i]
    C13_list.append(C13_value)
    
O16_list = []

for i in range(min([len(R17_list), len(R18_list)])):
    # Eq. 6
    O16_value = 1/(1+R17_list[i]+R18_list[i])
    O16_list.append(O16_value)

O17_list = []

for i in range(min([len(O16_list), len(R17_list)])):
    # Eq. 7
    O17_value = O16_list[i]*R17_list[i]
    O17_list.append(O17_value)

O18_list = []

for i in range(min([len(O16_list), len(R18_list)])):
    # Eq. 8
    O18_value = O16_list[i]*R18_list[i]
    O18_list.append(O18_value)

""" In this section the stochastic distributions are calculated for R45, R46 and R47, based on the concentration data from the previous section. """

R45_sto_list = []

for i in range(min([len(C12_list), len(C13_list), len(O16_list), len(O17_list)])):
    # Eq. 9
    R45_sto_value = (C13_list[i]*O16_list[i]**2+2*C12_list[i]*O16_list[i]*O17_list[i])/(C12_list[i]*O16_list[i]**2)
    R45_sto_list.append(R45_sto_value)

R46_sto_list = []

for i in range(min([len(C12_list), len(C13_list), len(O16_list), len(O17_list), len(O18_list)])):
    # Eq. 10
    R46_sto_value = (2*C12_list[i]*O16_list[i]*O18_list[i] + C12_list[i]*O17_list[i]**2 + 2*C13_list[i]*O16_list[i]*O17_list[i])/(C12_list[i]*O16_list[i]**2)
    R46_sto_list.append(R46_sto_value)

R47_sto_list = []

for i in range(min([len(C12_list), len(C13_list), len(O16_list), len(O17_list), len(O18_list)])):
    # Eq. 11
    R47_sto_value = (2*C13_list[i]*O16_list[i]*O18_list[i] + C13_list[i]*O17_list[i]**2 + 2*C12_list[i]*O17_list[i]*O18_list[i])/(C12_list[i]*O16_list[i]**2)
    R47_sto_list.append(R47_sto_value)

R48_sto_list = []

for i in range(min([len(C12_list), len(C13_list), len(O16_list), len(O17_list)])):
    # Eq. 9
    R48_sto_value = (C12_list[i]*O18_list[i]**2+2*C13_list[i]*O17_list[i]*O18_list[i])/(C12_list[i]*O16_list[i]**2)
    R48_sto_list.append(R48_sto_value)

R49_sto_list = []

for i in range(min([len(C12_list), len(C13_list), len(O16_list)])):
    # Eq. 9
    R49_sto_value = (C13_list[i]*O18_list[i]**2)/(C12_list[i]*O16_list[i]**2)
    R49_sto_list.append(R49_sto_value)

#Scrambled gas isotopologue abundances

# Eq. 12
d45_component_1 = ((R45_sto_list[2]/R45_sto_list[0])-1)*1000

d45_component_2 = ((R45_sto_list[3]/R45_sto_list[0])-1)*1000

# Eq. 13
d46_component_1 = ((R46_sto_list[2]/R46_sto_list[0])-1)*1000

d46_component_2 = ((R46_sto_list[3]/R46_sto_list[0])-1)*1000

# Eq. 14
#def D47_RF (D47EndMemberValue, acid_frac_alpha):
#    y = D47EndMemberValue - acid_frac_alpha
#    return y

# Eq. 15
#def D47_SGvsWG_0 (D47_RF, ETF_int, ETF_slope):
#    y = (D47_RF - ETF_int)/ETF_slope
#    return y

def remove_acidfrac_tf (componentD47_value):
    #Eq. 14
    removed_acidfrac_offset = componentD47_value - acid_frac_offset
    #Eq. 15
    removed_transferfunction = (removed_acidfrac_offset - transfer_function_intercept_D47)/transfer_function_slope_D47
    return removed_transferfunction

removed_acid_tf_c1 = remove_acidfrac_tf(component_1_D47)
removed_acid_tf_c2 = remove_acidfrac_tf(component_2_D47)

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
