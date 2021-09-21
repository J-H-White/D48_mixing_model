 # -*- coding: utf-8 -*-
"""
Created on Mon Aug  2 13:31:30 2021

@author: Jackson White
This model describes the non-linear mixing effects on mass-47 CO2 clumped isotope data. This code is written based on Defliese & Lohmann (2015)'s work. Here, we try to extend this model to the non-linear mixing effects on mass-48 CO2 clumped isotope data. To make reading easier, we recommend wrapping the code lines in Spyder by going to Preferences -> Editor -> select the Display tab -> tick wrap lines.

References

Defliese, W.F., Lohmann, K.C., 2015. Non-linear mixing effects on mass-47 CO2clumped isotope thermometry: Patterns and implications. Rapid Communications in Mass Spectrometry 29, 901–909.. doi:10.1002/rcm.7175
"""

#Import packages that are not Built-in functions

from math import e

""" In this section, the input the values for the isotope geochemistry data for the two components (component_1_D47, component_1_D48, component_1_d18O_VPDB, component_1_d13C_VPDB, component_1_contribution, component_2_D47, component_2_D48, component_2_d18O_VPDB, component_2_d13C_VPDB); transfer function slopes and intercepts (empirical_transfer_function_slope_D47, empirical_transfer_function_intercept_D47, empirical_transfer_function_slope_D48, empirical_transfer_function_intercept_D48); heated gas line slopes (heated_gas_line_slope_D47 and heated_gas_line_slope_D48); working gas delta values (working_gas_d13C_VPDB and working_gas_d18O_VSMOW); and the acid reaction temperature (T). """  

def mixing_model(component_1_D47, component_1_D48, component_1_d18O_VPDB, component_1_d13C_VPDB, component_1_contribution, component_2_D47, component_2_D48, component_2_d18O_VPDB, component_2_d13C_VPDB, T, empirical_transfer_function_slope_D47, empirical_transfer_function_intercept_D47, empirical_transfer_function_slope_D48, empirical_transfer_function_intercept_D48, heated_gas_line_slope_D47, heated_gas_line_slope_D48, reference_gas_d13C_VPDB, reference_gas_d18O_VSMOW):
    """ d18O_alpha_value is defnied as the oxygen fractionation from the conversion from CaCO3 to CO2 at a temperature T. """
    d18O_alpha = 1.00397+((5.25*10**2)/(273.15+T)**2)

    """ acid_frac_alpha is fractionation due to acid digestion. """
    acid_frac_alpha = e**((22.434/((273.15 + T)**2))-0.0002524)

    acid_frac_offset = (1000/acid_frac_alpha)-1000

    def mixture (C1, C2, contribution_1):
    
        """ mixture defines the weighted average given the inputs. """
    
        # Eq. 17
        mixture_value = (C1*contribution_1) + (C2*(1-contribution_1))
        return mixture_value

    # d13C mixture value from component 1 and component 2
    mixture_d13C_VPDB = mixture(component_1_d13C_VPDB, component_2_d13C_VPDB, component_1_contribution)

    # d18O mixture value from component 1 and component 2
    mixture_d18O_VPDB = mixture(component_1_d18O_VPDB, component_2_d18O_VPDB, component_1_contribution)

    # Convert d18O values from the VPDB scale to VSMOW scale and corrected for CaCO3 to CO2 fractionation (d18O_alpha_value).
    def CO3VPDB_to_CO2VSMOW(d18O_VPDB_value):
        CO3_VSMOW_value = ((d18O_VPDB_value + 29.98)/0.97001)
        CO2_VSMOW_value = ((1000+CO3_VSMOW_value)*d18O_alpha)-1000
        return CO2_VSMOW_value

    mixture_d18O_CO2VSMOW = CO3VPDB_to_CO2VSMOW(mixture_d18O_VPDB)

    component_1_d18O_CO2VSMOW = CO3VPDB_to_CO2VSMOW(component_1_d18O_VPDB)

    component_2_d18O_CO2VSMOW = CO3VPDB_to_CO2VSMOW(component_2_d18O_VPDB)

    """ Lists have been created to perform calculations on the equations will preform calculations on. Python has different data structures of which lists have the property that they can be changed, which is what our intention is by calculating with these equations below. After these calculations are done, the 'keys' list will be merged with the two 'values' lists to create a dictionary. Dictionaries are data structures which have a key, which can be called, and results in a output. """

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

    def remove_acidfrac_etf (componentD47_value, remove_acidfrac = True):
    
        """ This function will correct the input value for the acid fractionation offset and empirical transfer function. By default remove_acidfrac = True when there is a 18O/16O fractionation, hence when using the function for D47 calculations. However, remove_acidfrac = False, when calculating D48 or D49 since the values measured do not fractionate since they have to loose the 16O atom. """
    
        if remove_acidfrac == True:
            #Eq. 14
            removed_acidfrac_offset = componentD47_value - acid_frac_offset
            #Eq. 15
            removed_empiricaltransferfunction = (removed_acidfrac_offset - empirical_transfer_function_intercept_D47)/empirical_transfer_function_slope_D47
        
        elif remove_acidfrac == False:
            #Eq. 15
            removed_empiricaltransferfunction = (componentD47_value - empirical_transfer_function_intercept_D48)/empirical_transfer_function_slope_D48
        
        return removed_empiricaltransferfunction

    def d4i(c, c_R4i, ref_R4i, heated_gas_line_slope = 47):
    
        """ This function calculates d4i (d47 or d48) from the component c which has had an acid correction and transfer function correction performed on previously. It also builds on the R47 and R48 stochastic values for a given component, c_R4i (R47_sto_list[2] or R48_sto_list[2]; and R47_sto_list[3] or R48_sto_list[3]) with its respective reference gas value, ref_R4i (R47_sto_list[0] or R48_sto_list[0]). The heated_gas_line_slope_D47 is a default value defined by the user at the begining of the script. However, when the heated_gas_line_slope = 48, the heated_gas_line_slope_D48 value will be used instead, granted the the approapriate D48 inputs have been specified. """
    
        if heated_gas_line_slope == 47:
           # Eq. 16
           numerator = ((c+1000)*c_R4i/(1000*ref_R4i))-1
           denominator = (1/1000)-(heated_gas_line_slope_D47*c_R4i/(1000*ref_R4i))
    
        elif heated_gas_line_slope == 48:
           numerator = ((c+1000)*c_R4i/(1000*ref_R4i))-1
           denominator = (1/1000)-(heated_gas_line_slope_D48*c_R4i/(1000*ref_R4i))
       
        return numerator/denominator 

    d47_component_1 = d4i(remove_acidfrac_etf(component_1_D47), R47_sto_list[2], R47_sto_list[0])
    d47_component_2 = d4i(remove_acidfrac_etf(component_2_D47), R47_sto_list[3], R47_sto_list[0])

    d48_component_1 = d4i(remove_acidfrac_etf(component_1_D48, remove_acidfrac = False), R48_sto_list[2], R48_sto_list[0], heated_gas_line_slope=48) 
    d48_component_2 = d4i(remove_acidfrac_etf(component_2_D48, remove_acidfrac = False), R48_sto_list[3], R48_sto_list[0], heated_gas_line_slope=48) 

    d45_mixture = mixture(d45_component_1, d45_component_2, component_1_contribution)
    d46_mixture = mixture(d46_component_1, d46_component_2, component_1_contribution)
    d47_mixture = mixture(d47_component_1, d47_component_2, component_1_contribution)
    d48_mixture = mixture(d48_component_1, d48_component_2, component_1_contribution)
    
    def R4i(d4i_mixture, R4i_sto_list):
        # Eq. 18
        R4i_mixture = ((d4i_mixture/1000)+1)*R4i_sto_list[0]
        return R4i_mixture
        
    R45_mixture = R4i(d45_mixture, R45_sto_list)
    R46_mixture = R4i(d46_mixture, R46_sto_list)
    R47_mixture = R4i(d47_mixture, R47_sto_list)
    R48_mixture = R4i(d48_mixture, R48_sto_list)
    
    def D4i(R4i_mixture, R4i_sto_list):
        D4i_mixture = ((R4i_mixture/R4i_sto_list[1])-1)*1000
        return D4i_mixture
    
    D45_mixing = D4i(R45_mixture, R45_sto_list)
    D46_mixing = D4i(R46_mixture, R46_sto_list)
    D47_mixing = D4i(R47_mixture, R47_sto_list)
    D48_mixing = D4i(R48_mixture, R48_sto_list)
    
    # Eq. 19
    D47_full = D47_mixing - D46_mixing - D45_mixing

    # Eq. 20
    D48_full = D48_mixing - 2*D46_mixing

    # Eq. 21
    def add_etf_acidfrac(D47 = True):
    
        """If D47 = True, then the according D47 parameters are passed. If D47 = False, then the according D48 parameters are passed. Note no acid_frac_offset parameter is required as there is no acid fractionation for mass 48 isotopologues, since mass 48 CO2 can only be measured by removing the 16O atom from the CO3 group."""
    
        if D47 == True:
            ETF_step_1 = D47_full - heated_gas_line_slope_D47 * d47_mixture
            ETF_step_2 = (ETF_step_1)*empirical_transfer_function_slope_D47 + empirical_transfer_function_intercept_D47
            acid_frac = ETF_step_2 + acid_frac_offset
            return acid_frac
    
        elif D47 == False:
            ETF_step_1 = D48_full - heated_gas_line_slope_D48 * d48_mixture
            ETF_step_2 = (ETF_step_1)*empirical_transfer_function_slope_D48 + empirical_transfer_function_intercept_D48
            return ETF_step_2
    
    export_line = {"component 1 contribution": component_1_contribution, "component 2 contribution": 1-component_1_contribution, "d18O mix": mixture_d18O_VPDB, "d13C mix": mixture_d13C_VPDB, "D47 model": add_etf_acidfrac(D47 = True), "D47 linear": mixture(component_1_D47, component_2_D47, component_1_contribution), "D48 model": add_etf_acidfrac(D47 = False), "D48 linear": mixture(component_1_D48, component_2_D48, component_1_contribution), "d45 mix": d45_mixture, "d46 mix": d46_mixture, "d47 mix": d47_mixture, "d48 mix": d48_mixture, "G47 mix": add_etf_acidfrac(D47 = True) - mixture(component_1_D47, component_2_D47, component_1_contribution), "G48 mix": add_etf_acidfrac(D47 = False) - mixture(component_1_D48, component_2_D48, component_1_contribution)}
    
    return print(export_line)
