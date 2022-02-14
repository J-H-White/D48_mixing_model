"""
Created on Mon Aug  2 13:31:30 2021

@author: Jackson White
This model describes the non-linear mixing effects on mass-48 CO2 clumped isotope data. This code is written based on Defliese & Lohmann (2015)'s work. Here, we try to extend this model to the non-linear mixing effects on mass-48 CO2 clumped isotope data. To make reading easier, we recommend wrapping the code lines in Spyder by going to Preferences -> Editor -> select the Display tab -> tick wrap lines.

Assumptions

This mixing model does not include a Δ48 acid fractionation factor, which is assumed to be negligible.

References

Defliese, W.F., Lohmann, K.C., 2015. Non-linear mixing effects on mass-47 CO2clumped isotope thermometry: Patterns and implications. Rapid Communications in Mass Spectrometry 29, 901–909.. doi:10.1002/rcm.7175
""" 

def mixing_model(component_1_Δ47, component_1_Δ48, component_1_δ18O_VPDB, component_1_δ13C_VPDB, component_1_contribution, component_2_Δ47, component_2_Δ48, component_2_δ18O_VPDB, component_2_δ13C_VPDB, temperature, empirical_transfer_function_slope_Δ47, empirical_transfer_function_intercept_Δ47, empirical_transfer_function_slope_Δ48, empirical_transfer_function_intercept_Δ48, standard_gas_line_slope_Δ47, standard_gas_line_slope_Δ48, reference_gas_δ13C_VPDB, reference_gas_δ18O_VSMOW):
    
    """
    This model describes the non-linear mixing effects on mass-47 CO2 clumped isotope data. This code was written based on Defliese & Lohmann (2015)'s work. Here, we try to extend this model to the non-linear  mixing effects on mass-48 CO2 clumped isotope data. To make the script reading easier, we recommend wrapping the code lines in Spyder by going to Preferences -> Editor -> select the Display tab -> tick wrap lines.
    
    For this model to run, it will require the following data:
        
        1. D47, for component_1_D47 and component_2_D47.
        2. D48, for component_1_D48 and component_2_D48.
        3. δ13C, for component_1_δ13C_VPDB, component_2_δ13C_VPDB.
        4. δ18O, for component_1_δ18O_VPDB, component_2_δ18O_VPDB.
        5. the contribution for component 1, (component_1_contribution) which is a value between 0 and 1, and contribution of component 2 is computed based on the difference between 1 and component_1_contribution.
        6. temperature (in degrees Celsius).
        7. the slope and interecept from the empirical transfer function for D47 (empirical_transfer_function_slope_D47, empirical_transfer_function_intercept_D47 respectively). 
        8. the slope and interecept from the empirical transfer function for D48 (empirical_transfer_function_slope_D48, empirical_transfer_function_intercept_D48 respectively). 
        9. the standard gas line slope for D47 and D48 (heated_gas_line_slope_D47, heated_gas_line_slope_D48).
        10. the reference gas δ13C and δ18O composition (reference_gas_δ13C_VPDB, reference_gas_δ18O_VSMOW). 

    All mentioned references are found at the end of this script.
    """
    
    def mixing_model_calculation(component_1_Δ47, component_1_Δ48, component_1_δ18O_VPDB, component_1_δ13C_VPDB, component_1_contribution, component_2_Δ47, component_2_Δ48, component_2_δ18O_VPDB, component_2_δ13C_VPDB, temperature, empirical_transfer_function_slope_Δ47, empirical_transfer_function_intercept_Δ47, empirical_transfer_function_slope_Δ48, empirical_transfer_function_intercept_Δ48, standard_gas_line_slope_Δ47, standard_gas_line_slope_Δ48, reference_gas_δ13C_VPDB, reference_gas_δ18O_VSMOW):

        """ This function defines all the calculations necessary compute Δ48 and Δ47 based on sample heterogeneities between two samples. """        

        # The variable "δ18O_acid_value" is defined as the oxygen fractionation coeffcient between converting CaCO3 to CO2 at a given temperature. This equation is defined in equation (6) of Swart et al., (1991).
        
        δ18O_acid_α = 1.00397 + ((5.25 * 10**2) / (273.15 + temperature) ** 2)
        
        # The variable "Δ47_acid_α" is the alpha fractionation coefficient (α) for Δ47 during acid digestion which is temperature dependant (Petersen et al., 2019). Δ47_acid_α is defined by the difference between a given carbonate sample not loosing any of its 13C–18O bonds assuming a carbonate group of interest is the isotopologue 13C18O16O2.   
        
        Δ47_acid_α = (0.0383 * 10**-6) * (10**6 / ((273.15 + temperature) ** 2)) + (0.258 * 10**-5)

        # The variable "acid_frac_offset" is the difference between the sample's true Δ47 value and that measured after the carbonate has undergone acid digestion. acid_frac_offset is calculated from the ratio version.
        
        acid_frac_offset = (1 / Δ47_acid_α - 1) * 1000

        def mixture (C1, C2, contribution_1):
    
            """ The "mixture" function defines a weighted average based on the fraction of a corresponding endmember, that is a value between 0 and 1. For example, component 1 (C1) may encompase 0.6 of the sample, therefore the contribution of C1 to the overall sample is contribution_1 = 0.6. The contribution of the other endmember, component 2 (C2) can be calculated as 1 - contribution_1, or 1 - 0.6 = 0.4. Hence, for this function only necessary variables are C1, C2 and contribution_1 are necesarry as inputs. """
            
            mixture_value = (C1 * contribution_1) + (C2 * (1 - contribution_1))
            return mixture_value

        # The variable "mixture_δ13C_VPDB" is calculated from the "mixture" function defined previously by using the d13C inputs for component 1 and 2 (i.e. component_1_δ13C_VPDB and component_2_δ13C_VPDB, respectively) and the contribution to the overall sample from component 1 (component_1_contribution).
        
        mixture_δ13C_VPDB = mixture(component_1_δ13C_VPDB, component_2_δ13C_VPDB, component_1_contribution)

        # The variable "mixture_δ18O_VPDB" is used in the same way as the previous comment, but takes the corresponding δ18O values instead.
        
        mixture_δ18O_VPDB = mixture(component_1_δ18O_VPDB, component_2_δ18O_VPDB, component_1_contribution)

        # "VPDB_to_acidVSMOW" converts δ18O values from the VPDB scale to VSMOW scale and corrects the output for the acid fractionation associated at the given temperature (δ18O_acid_α).
        
        def VPDB_to_acidVSMOW(δ18O_VPDB_value):
            
            # Scale conversion from VPDB to VSMOW (Kim et al. (2015)).
            
            VSMOW_value = 30.92 + 1.03092 * δ18O_VPDB_value
            
            # Acid fractionation correction determined by the specified temperature, and redefined in terms of the δ definition (Hoefs, 2018).
            
            acidVSMOW_value = ((1000 + VSMOW_value) * δ18O_acid_α) - 1000
            return acidVSMOW_value

        mixture_δ18O_CO2VSMOW = VPDB_to_acidVSMOW(mixture_δ18O_VPDB)

        component_1_δ18O_CO2VSMOW = VPDB_to_acidVSMOW(component_1_δ18O_VPDB)

        component_2_δ18O_CO2VSMOW = VPDB_to_acidVSMOW(component_2_δ18O_VPDB)

        # The lists created below —δ13C_values and δ18O_values— save the output of the equations defined above. Python has different data structures of which lists have the property that they modified.

        δ13C_values = [reference_gas_δ13C_VPDB, mixture_δ13C_VPDB, component_1_δ13C_VPDB, component_2_δ13C_VPDB]

        δ18O_values = [reference_gas_δ18O_VSMOW, mixture_δ18O_CO2VSMOW, component_1_δ18O_CO2VSMOW, component_2_δ18O_CO2VSMOW]

        # In the following section the δ13C and δ18O values are converted into R13, R18 and R17 values. To perform the conversion, the absolute abundance of heavy isotopes for the 13C standard, Vienna Pee Dee Belemnite (VPDB), and 18O standard, Vienna Standard Mean Ocean Water (VSMOW), must be defined. The absolute abundances are expressed in ratios, R13(VPDB), for 13C/12C, R17(VSMOW) and R18(VSMOW), for 18O/16O; in addition to the triple oxygen line, λ, to convert R18(VSMOW) to R17(VSMOW) (Petersen et al., (2019)). The absolute values are the following:
        
        # 1. R13(VPDB) = 0.011180
        # 2. R17(VSMOW) = 0.038475
        # 3. λ = 0.528
        # 4. R18(VSMOW) = 0.0020052
        # (Brand et al.(2010))

        R13_VPDB = 0.011180
        R17_VSMOW = 0.038475
        λ = 0.528
        R18_VSMOW = 0.0020052 

        R13_list = []

        for i in range(len(δ13C_values)):
            R13_values = ((δ13C_values[i] / 1000) + 1) * R13_VPDB
            R13_list.append(R13_values)

        R17_list = []

        for i in range(len(R17_list)):
            R17_values = (((δ18O_values[i] / 1000) + 1) ** λ) * R17_VSMOW
            R17_list.append(R17_values)

        R18_list = []

        for i in range(len(δ18O_values)):
            R18_values = ((δ18O_values[i] / 1000) + 1) * R18_VSMOW
            R18_list.append(R18_values)

        # Based on the R13, R18 and R17 data, the concentrations for each isotopes can be estimated. That is 12C, 13C, 16O, 17O, and 18O. The values and letters had to be switched because 12C, for exmaple, can produce syntax errors when calculations are being preformed.

        C12_list = []

        for i in range(len(R13_list)):
            C12_value = 1 / (1 + R13_list[i])
            C12_list.append(C12_value)

        C13_list = []

        for i in range(min([len(R13_list), len(C12_list)])):
            C13_value = C12_list[i] * R13_list[i]
            C13_list.append(C13_value)
    
        O16_list = []

        for i in range(min([len(R17_list), len(R18_list)])):
            O16_value = 1 / (1 + R17_list[i] + R18_list[i])
            O16_list.append(O16_value)

        O17_list = []

        for i in range(min([len(O16_list), len(R17_list)])):
            O17_value = O16_list[i] * R17_list[i]
            O17_list.append(O17_value)

        O18_list = []

        for i in range(min([len(O16_list), len(R18_list)])):
            O18_value = O16_list[i] * R18_list[i]
            O18_list.append(O18_value)

        # In this section, the expected concentrations for each isotopologue are estimated probabilistically. These are calculated for R45, R46 and R47, based on the relative concentration data from the previous section.

        R45_sto_list = []

        for i in range(min([len(C12_list), len(C13_list), len(O16_list), len(O17_list)])):
            
            # A mass 45 CO2 isotopologue can result from the isotope combination 16O-13C-16O or 17O-12C-16O. These isotopologues can be determined by multiplying the corresponding values for each isotope from the previous code. The multiplication can be performed since all isotope values are normalized. Hence, the results from the multiplications of these isotope concentrations are somewhere in the range between zero and one, and are stochastic or probabilistic.
            
            # The '2' multiplying 17O-12C-16O arises from the fact that the 17O in a carbonate group can be in one of two possible positions when CO2 is formed from acid digestion. By 'positions' it is implied that a 16O was released during acid digestion, leaving two posible potential sites in the carbonate group. A factor '2' multiplication is found multiplying 16O-12C-18O and 17O-13C-16O in R46, 16O-13C-18O and 17O-12C-18O in R47, and 17O-13C-18O in R48. See appendix for further details.   
            
            R45_sto_value = (C13_list[i] * O16_list[i] ** 2 + 2 * C12_list[i] * O16_list[i] * O17_list[i]) /(C12_list[i] * O16_list[i] ** 2)
            R45_sto_list.append(R45_sto_value)

        R46_sto_list = []

        for i in range(min([len(C12_list), len(C13_list), len(O16_list), len(O17_list), len(O18_list)])):
            
            #  A mass 46 CO2 isotopologue can result from the combinations: 16O-12C-18O, 17O-12C-17O or 17O-13C-16O.
            
            R46_sto_value = (2 * C12_list[i] * O16_list[i] * O18_list[i] + C12_list[i] * O17_list[i] ** 2 + 2 * C13_list[i] * O16_list[i] * O17_list[i]) / (C12_list[i] * O16_list[i] ** 2)
            R46_sto_list.append(R46_sto_value)

        R47_sto_list = []

        for i in range(min([len(C12_list), len(C13_list), len(O16_list), len(O17_list), len(O18_list)])):
            
            # A mass 47 CO2 isotopologue can result from the combinations: 16O-13C-18O, 17O-13C-17O or 17O-12C-18O.
            
            R47_sto_value = (2 * C13_list[i] * O16_list[i] * O18_list[i] + C13_list[i] * O17_list[i] ** 2 + 2 * C12_list[i] * O17_list[i] * O18_list[i]) / (C12_list[i] * O16_list[i] ** 2)
            R47_sto_list.append(R47_sto_value)

        R48_sto_list = []

        for i in range(min([len(C12_list), len(C13_list), len(O16_list), len(O17_list)])):
            
            # A mass 48 CO2 isotopologue can result from the combinations: 18O-12C-18O or 17O-13C-18O.
            
            R48_sto_value = (C12_list[i] * O18_list[i] ** 2 + 2 * C13_list[i] * O17_list[i] * O18_list[i])/(C12_list[i] * O16_list[i] ** 2)
            R48_sto_list.append(R48_sto_value)

        R49_sto_list = []

        for i in range(min([len(C12_list), len(C13_list), len(O16_list)])):
            
            # A mass 48 CO2 isotopologue can result from the combination 18O-13C-18O.
            
            R49_sto_value = (C13_list[i] * O18_list[i] ** 2) / (C12_list[i] * O16_list[i] ** 2)
            R49_sto_list.append(R49_sto_value)

        # [0], [2] and [3] are the reference gas, component 1 and component 2 values from the according R45 and R46 stochastic values.
        
        δ45_component_1 = ((R45_sto_list[2] / R45_sto_list[0]) - 1) * 1000

        δ45_component_2 = ((R45_sto_list[3] / R45_sto_list[0]) - 1) * 1000

        δ46_component_1 = ((R46_sto_list[2] / R46_sto_list[0]) - 1) * 1000

        δ46_component_2 = ((R46_sto_list[3] / R46_sto_list[0]) - 1) * 1000

        def remove_acidfrac_etf (componentΔ4i_value, remove_acidfrac = 47):
    
            """ This function corrects the input value for the offset caused by the acid fractionation during carbonate digestion. In addition, the will account for mass spectrometric artifacts by applying the empirical transfer function. By default, the acid fractionation correction, remove_acidfrac = 47, is required when Δ47 calculations are performed. Although fractionation does occur when calculating Δ48 or Δ49, remove_acidfrac = 48, we assume it is neglegible; however we do not know this for sure. """
    
            if remove_acidfrac == 47:
                
                # A raw measurement for a component is corrected by the fractionation occured during acid digestion. 
                
                removed_acidfrac_offset = componentΔ4i_value - acid_frac_offset
                
                # The corrected value, that is remove_acidfrac_offset, is standardized to the absolute reference frame (Dennis et al., 2011).
                
                removed_etf = (removed_acidfrac_offset - empirical_transfer_function_intercept_Δ47)/empirical_transfer_function_slope_Δ47
        
            elif remove_acidfrac == 48:
                
                # This conditional assumes that no acid fractionation is needed here. Therefore, we procede to standardize into to the absolute reference frame (Dennis et al., 2011).
                
                removed_etf = (componentΔ4i_value - empirical_transfer_function_intercept_Δ48)/empirical_transfer_function_slope_Δ48
        
            return removed_etf

        def δ4i(c, c_R4i, ref_R4i, standard_gas_line_slope = 47):
    
            """ The δ4i function which correlates the Δ47 or Δ48 sample measurements, to a temperature-Δ47 or Δ48 calibration based on standards. These standards are heated CO2 gas. This step allow to determine the final clumped measurement. The function calculates δ4i (δ47 or δ48) from the component c which has had an acid correction and an empirical transfer function correction performed on previously. The function builds on the R47 and R48 stochastic values for a given component, c_R4i. c_R4i corresponds to component 1 (R47_sto_list[2] or R48_sto_list[2]) and component 2 (R47_sto_list[3] or R48_sto_list[3]) values with the respective reference gas value, ref_R4i (R47_sto_list[0] or R48_sto_list[0]). The standard_gas_line_slope_D47 is a default value defined by the user at the begining of the script. However, when the standard_gas_line_slope = 48, the heated_gas_line_slope_D48 value will be used instead, granted the the approapriate D48 inputs have been specified. """
    
            if standard_gas_line_slope == 47:
               
               numerator = ((c + 1000) * c_R4i / (1000 * ref_R4i)) - 1
               denominator = (1 / 1000) - (standard_gas_line_slope_Δ47 * c_R4i / (1000 * ref_R4i))
    
            elif standard_gas_line_slope == 48:
               numerator = ((c + 1000) * c_R4i / (1000 * ref_R4i)) - 1
               denominator = (1 / 1000) - (standard_gas_line_slope_Δ48 * c_R4i / (1000 * ref_R4i))
       
            return numerator/denominator 

        # δ47 and δ48 are the component values once the acid fractionation has been applied and then ploted on the absolute reference frame. This allows the user to mix sample values which have been performed on different mass spectrometers, or different calibrations of the empirical transfer function.

        δ47_component_1 = δ4i(remove_acidfrac_etf(component_1_Δ47), R47_sto_list[2], R47_sto_list[0])
        δ47_component_2 = δ4i(remove_acidfrac_etf(component_2_Δ47), R47_sto_list[3], R47_sto_list[0])
        δ48_component_1 = δ4i(remove_acidfrac_etf(component_1_Δ48, remove_acidfrac = 48), R48_sto_list[2], R48_sto_list[0], heated_gas_line_slope = 48) 
        δ48_component_2 = δ4i(remove_acidfrac_etf(component_2_Δ48, remove_acidfrac = 48), R48_sto_list[3], R48_sto_list[0], heated_gas_line_slope = 48) 

        # These mixture values are the computations by a weighted average.

        δ45_mixture = mixture(δ45_component_1, δ45_component_2, component_1_contribution)
        δ46_mixture = mixture(δ46_component_1, δ46_component_2, component_1_contribution)
        δ47_mixture = mixture(δ47_component_1, δ47_component_2, component_1_contribution)
        δ48_mixture = mixture(δ48_component_1, δ48_component_2, component_1_contribution)
    
        def R4i(δ4i_mixture, R4i_sto_list):
            
            # This function coverts δ45, δ46, δ47, and δ48  mixture values to the equivalent ratios. Those are the values R45_mixture, R46_mixture, R47_mixture, and R48_mixture
            
            R4i_mixture = ((δ4i_mixture / 1000) + 1) * R4i_sto_list[0]
            return R4i_mixture
        
        R45_mixture = R4i(δ45_mixture, R45_sto_list)
        R46_mixture = R4i(δ46_mixture, R46_sto_list)
        R47_mixture = R4i(δ47_mixture, R47_sto_list)
        R48_mixture = R4i(δ48_mixture, R48_sto_list)
    
        def Δ4i(R4i_mixture, R4i_sto_list):
            
            # This function takes the corresponding ratio of a mixture (R45_mixture, R46_mixture, R47_mixture or R48_mixture), substracts the expected stochastic value, normalizes the difference and then expresses the final result in per mille.
            
            Δ4i_mixture = ((R4i_mixture / R4i_sto_list[1]) - 1) * 1000
            return Δ4i_mixture
    
        Δ45_mixing = Δ4i(R45_mixture, R45_sto_list)
        Δ46_mixing = Δ4i(R46_mixture, R46_sto_list)
        Δ47_mixing = Δ4i(R47_mixture, R47_sto_list)
        Δ48_mixing = Δ4i(R48_mixture, R48_sto_list)
    
        # The Δ47_mixing, is substracted by Δ46_mixing and Δ45_mixing to determine the Δ47 value.
        
        Δ47_no_reference_frame = Δ47_mixing - Δ46_mixing - Δ45_mixing

        # The Δ48_mixing, is substracted by two times Δ46_mixing.
        
        Δ48_no_reference_frame = Δ48_mixing - 2*Δ46_mixing

        def add_etf_acidfrac(Δ47 = 47):
    
            """If D47 = True, then the according D47 parameters are passed. If D47 = False, then the according D48 parameters are passed. Note no acid_frac_offset parameter as it is assumed that no acid fractionation for mass 48 isotopologues, since mass 48 CO2 can only be measured by removing the 16O atom from the CO3 group."""
    
            if Δ47 == 47: #This is the set default value
                ETF_step_1 = Δ47_no_reference_frame - standard_gas_line_slope_Δ47 * δ47_mixture
                ETF_step_2 = (ETF_step_1) * empirical_transfer_function_slope_Δ47 + empirical_transfer_function_intercept_Δ47
                acid_frac = ETF_step_2 + acid_frac_offset
                return acid_frac
    
            elif Δ47 == 48:
                ETF_step_1 = Δ48_no_reference_frame - standard_gas_line_slope_Δ48 * δ48_mixture
                ETF_step_2 = (ETF_step_1) * empirical_transfer_function_slope_Δ48 + empirical_transfer_function_intercept_Δ48
                return ETF_step_2
        
        # The export line is a dictionary, where the "keys" or names in "" are what the number represents at after the colon. The export line will be taking into a DataFrame, which will then be converted to an Excel file. 
        
        export_line = {"component 1 δ13C (VPDB)": component_1_δ13C_VPDB, "component 1 δ18O (VPDB)":  component_1_δ18O_VPDB, "component 1 δ18O (VSMOW)": component_1_δ18O_CO2VSMOW, "component 1 contribution": component_1_contribution, "component 2 δ13C (VPDB)": component_2_δ13C_VPDB, "component 2 δ18O (VPDB)": component_2_δ18O_VPDB, "component 2 δ18O (VSMOW)": component_2_δ18O_CO2VSMOW, "component 2 contribution": 1 - component_1_contribution, "δ13C mix": mixture_δ13C_VPDB, "δ18O mix": mixture_δ18O_VPDB, "Δ47 model": add_etf_acidfrac(Δ47 = 47), "Δ47 linear": mixture(component_1_Δ47, component_2_Δ47, component_1_contribution),   "Δ48 model": add_etf_acidfrac(Δ47 = 48), "Δ48 linear": mixture(component_1_Δ48, component_2_Δ48, component_1_contribution), "δ45 mix": δ45_mixture, "δ46 mix": δ46_mixture, "δ47 mix": δ47_mixture, "δ48 mix": δ48_mixture, "G47 mix": add_etf_acidfrac(Δ47 = 47) - mixture(component_1_Δ47, component_2_Δ47, component_1_contribution), "G48 mix": add_etf_acidfrac(Δ47 = 48) - mixture(component_1_Δ48, component_2_Δ48, component_1_contribution)}
    
        return export_line

    def mixing_model_data_simulation():
        
        df_list = []
        
        for i in range(len(component_1_contribution)):
            loop = mixing_model_calculation(component_1_Δ47, component_1_Δ48, component_1_δ18O_VPDB, component_1_δ13C_VPDB, component_1_contribution[i], component_2_Δ47, component_2_Δ48, component_2_δ18O_VPDB, component_2_δ13C_VPDB, temperature, empirical_transfer_function_slope_Δ47, empirical_transfer_function_intercept_Δ47, empirical_transfer_function_slope_Δ48, empirical_transfer_function_intercept_Δ48, standard_gas_line_slope_Δ47, standard_gas_line_slope_Δ48, reference_gas_δ13C_VPDB, reference_gas_δ18O_VSMOW)
            df_list.append(loop)
        
        # Import pandas
        import pandas as pd
        
        # Convert list of dictionaries to a pandas dataframe
        df = pd.DataFrame(df_list)
    
        # Export DataFrame as an excel spreadsheet
        return df.to_excel("output.xlsx")
    return mixing_model_data_simulation()

#References

# Brand, W. A., Assonov, S. S., & Coplen, T. B. (2010). Correction for the 17O interference in δ(13C) measurements when analyzing CO2 with stable isotope mass spectrometry (IUPAC Technical Report). Pure and Applied Chemistry, 82(8), 1719–1733. https://doi.org/10.1351/PAC-REP-09-01-05

# Defliese, W. F., Hren, M. T., & Lohmann, K. C. (2015). Compositional and temperature effects of phosphoric acid fractionation on Δ47 analysis and implications for discrepant calibrations. Chemical Geology, 396, 51–60. https://doi.org/10.1016/j.chemgeo.2014.12.018

# Defliese, W. F., Lohmann, K. C., (2015). Non-linear mixing effects on mass-47 CO2clumped isotope thermometry: Patterns and implications. Rapid Communications in Mass Spectrometry 29, 901–909.. doi:10.1002/rcm.7175

# Dennis, K. J., Affek, H. P., Passey, B. H., Schrag, D. P., & Eiler, J. M. (2011). Defining an absolute reference frame for 'clumped' isotope studies of CO2. Geochimica et Cosmochimica Acta, 75(22), 7117–7131. https://doi.org/10.1016/j.gca.2011.09.025

# Hoefs, J. (2018). Stable isotope geochemistry (Eighth edition ed.). Springer.

# Kim, S.-T., Coplen, T. B., & Horita, J. (2015). Normalization of stable isotope data for carbonate minerals; implementation of IUPAC guidelines. Geochimica et Cosmochimica Acta, 158, 276–289. https://doi.org/10.1016/j.gca.2015.02.011

# Petersen, S. V., Defliese, W. F., Saenger, C., Daëron, M., Huntington, K. W., John, C. M., et al (2019). Effects of Improved 17O Correction on Interlaboratory Agreement in Clumped Isotope Calibrations, Estimates of Mineral‐Specific Offsets, and Temperature Dependence of Acid Digestion Fractionation. Geochemistry, geophysics, geosystems : G3, 20(7), 3495-3519. https://doi.org/10.1029/2018GC008127

# Swart, P. K., Burns, S., & Leder, J. (1991). Fractionation of the stable isotops of oxygen and carbon in carbon dioxide during the reaction of calcite with phosphoric acid as a function of temperature and technique. Chemical Geology, 86(2), 89–96. https://doi.org/10.1016/0168-9622(91)90055-2
