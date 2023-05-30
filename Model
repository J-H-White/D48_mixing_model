 # -*- coding: utf-8 -*-
"""
Created on Mon Aug  2 13:31:30 2021

@author: Jackson White

This model describes the non-linear mixing effects on mass-48 CO2 clumped isotope data and the difference between the non-linear model values and  those from a linear, or weighted average, approach. This code is written based on Defliese & Lohmann (2015)'s work. Here, we try to extend this model to the non-linear mixing effects on mass-48 CO2 clumped isotope data. To make reading easier, we recommend wrapping the code lines in Spyder by going to Preferences -> Editor -> select the Display tab -> tick wrap lines.

Assumptions:

1. Δ47 and Δ48 inputs are in the absolute reference frame (Dennis et al., 2011).
2. The stochastic ratio values for the samples and working gas are assumed to be the same to those ratios directly measured by mass spectrometry for δ45 and δ46.

Limitations:

1. Due to a lack in studies, this model does not include a Δ48 acid fractionation coefficient, therefore the Δ48 acid fractionation coefficient is assumed to be negligible.
2. This model is limited to two end-members.
 
References can be found at the end of this script.

""" 

def mixing_model(endmember_1_Δ47, endmember_1_Δ48, endmember_1_δ18O_VPDB, endmember_1_δ13C_VPDB, endmember_1_contribution, endmember_2_Δ47, endmember_2_Δ48, endmember_2_δ18O_VPDB, endmember_2_δ13C_VPDB, temperature, empirical_transfer_function_slope_Δ47, empirical_transfer_function_intercept_Δ47, empirical_transfer_function_slope_Δ48, empirical_transfer_function_intercept_Δ48, heated_gas_line_slope_Δ47, heated_gas_line_slope_Δ48, working_gas_δ13C_VPDB, working_gas_δ18O_VSMOW):
    
    """
    This model describes the non-linear mixing effects on mass 47 and mass 48 CO2 clumped isotope data. This code was written based on Defliese & Lohmann (2015)'s work. Here, we try to extend their mass 47 mixing model to the non-linear mixing effects on mass-48 CO2 clumped isotope data. To make the script reading easier, we recommend wrapping the code lines in Spyder by going to Preferences -> Editor -> select the Display tab -> tick wrap lines.
    
    For this model to run, it will require the following data for both end-members or endmembers:
        
        1. Δ47 (endmember_1_Δ47 and endmember_2_Δ47).
        2. Δ48 (endmember_1_Δ48 and endmember_2_Δ48).
        3. δ13C (VPDB) (endmember_1_δ13C_VPDB, endmember_2_δ13C_VPDB).
        4. δ18O (VPDB) (endmember_1_δ18O_VPDB, endmember_2_δ18O_VPDB).
        5. The fraction of the sample size for endmember 1 (endmember_1_contribution)—expressed as a decimal—and the fractional contribution of endmember 2 is calculated by substracting endmember_1_contribution from 1.
        6. Temperature of acid digestion in degrees Celsius.
        7. The slopes and interecepts from the empirical transfer functions for Δ47 (empirical_transfer_function_slope_Δ47 and empirical_transfer_function_intercept_Δ47) and Δ48 (empirical_transfer_function_slope_Δ48 and empirical_transfer_function_intercept_Δ48). 
        8. The heated gas line slope for Δ47 and Δ48 (heated_gas_line_slope_Δ47, heated_gas_line_slope_Δ48).
        9. The working gas δ13C and δ18O composition (working_gas_δ13C_VPDB, working_gas_δ18O_VSMOW). 

    All mentioned references are found at the end of this script.
    """
    
    def mixing_model_calculation(endmember_1_Δ47, endmember_1_Δ48, endmember_1_δ18O_VPDB, endmember_1_δ13C_VPDB, endmember_1_contribution, endmember_2_Δ47, endmember_2_Δ48, endmember_2_δ18O_VPDB, endmember_2_δ13C_VPDB, temperature, empirical_transfer_function_slope_Δ47, empirical_transfer_function_intercept_Δ47, empirical_transfer_function_slope_Δ48, empirical_transfer_function_intercept_Δ48, heated_gas_line_slope_Δ47, heated_gas_line_slope_Δ48, working_gas_δ13C_VPDB, working_gas_δ18O_VSMOW):

        """ This function defines all the calculations necessary compute Δ48 and Δ47 based on sample heterogeneities between two samples. """        

        def mixture (C1, C2, contribution_1):
    
            """ The "mixture" function defines a weighted average based on the fraction of a corresponding endmember, that is a value between 0 and 1. For example, endmember 1 (C1) may encompase 0.6 of the sample, therefore the contribution of C1 to the overall sample is contribution_1 = 0.6. The contribution of the other endmember, endmember 2 (C2) can be calculated as 1 - contribution_1, or 1 - 0.6 = 0.4. Hence, for this function only necessary variables are C1, C2 and contribution_1 are necesarry as inputs. """
            
            mixture_value = (C1 * contribution_1) + (C2 * (1 - contribution_1)) # Equation (1)
            
            return mixture_value

        # The variable "mixture_δ13C_VPDB" is calculated from the "mixture" function defined previously by using the δ13C inputs for endmember 1 and 2 (i.e. endmember_1_δ13C_VPDB and endmember_2_δ13C_VPDB, respectively) and the contribution to the overall sample from endmember 1 (endmember_1_contribution).
        
        mixture_δ13C_VPDB = mixture(endmember_1_δ13C_VPDB, endmember_2_δ13C_VPDB, endmember_1_contribution)

        # The variable "mixture_δ18O_VPDB" is used in the same way as the previous comment, but takes the corresponding δ18O values instead.
        
        mixture_δ18O_VPDB = mixture(endmember_1_δ18O_VPDB, endmember_2_δ18O_VPDB, endmember_1_contribution)

        # "VPDB_to_acidVSMOW" converts δ18O values from the VPDB scale to VSMOW scale and corrects the output for the acid fractionation associated at the given temperature (δ18O_acid_α).
        
        def VPDB_to_acidVSMOW(δ18O_VPDB_value):
            
            # Scale conversion from VPDB to VSMOW (equation (3), Kim et al. (2015)).
            
            VSMOW_value = 1.03092 * δ18O_VPDB_value + 30.92 # Equation (2)
            
            # The variable "δ18O_acid_value" is defined as the oxygen fractionation coeffcient between converting CaCO3 to CO2 at a given temperature. This equation is defined in equation (6) of Swart et al., (1991).
            
            δ18O_acid_α = ((5.25 * 10**2) / (273.15 + temperature) ** 2) + 1.00397 # Equation (3) 
            
            # VSMOW_value is corrected by the acid fractionation coefficient, δ18O_acid_α, and the sample ratio is redefined in terms of the δ definition (Hoefs, 2018).
            
            acidVSMOW_value = (((VSMOW_value + 1000) * δ18O_acid_α) - 1000) # Equation (4)
            
            return acidVSMOW_value

        mixture_δ18O_CO2VSMOW = VPDB_to_acidVSMOW(mixture_δ18O_VPDB)

        endmember_1_δ18O_CO2VSMOW = VPDB_to_acidVSMOW(endmember_1_δ18O_VPDB)

        endmember_2_δ18O_CO2VSMOW = VPDB_to_acidVSMOW(endmember_2_δ18O_VPDB)

        # The lists created below—δ13C_values and δ18O_values—save the output of the equations defined above. Python has different data structures of which lists have the property that they modified.

        δ13C_values = [working_gas_δ13C_VPDB, mixture_δ13C_VPDB, endmember_1_δ13C_VPDB, endmember_2_δ13C_VPDB]

        δ18O_values = [working_gas_δ18O_VSMOW, mixture_δ18O_CO2VSMOW, endmember_1_δ18O_CO2VSMOW, endmember_2_δ18O_CO2VSMOW]

        # In the following section the δ13C and δ18O values are converted into R13, R18 and R17 values, to later estimate the stochastic concentration of CO2 isotopologues. To perform the conversion, the absolute abundance of heavy isotopes for the 13C standard, Vienna Pee Dee Belemnite (VPDB), and 18O standard, Vienna Standard Mean Ocean Water (VSMOW), must be defined. The absolute abundances are expressed in ratios, R13(VPDB), for 13C/12C, R17(VSMOW) and R18(VSMOW), for 18O/16O; in addition to the triple oxygen line, λ, to convert R18(VSMOW) to R17(VSMOW) (Petersen et al., (2019)). The absolute values are taken from Table 2 in Brand et al. (2010) for the 17O correction:

        R13_VPDB = 0.011180 # Ratio value of the NBS19 standard, evolved from 25⁰C digestion with 100% H3PO4 (Zhang et al., 1990). 
        R17_VSMOW = 0.0003931 # Is the 17R value proposed by 17O-correction algorithms and graphite combustion experiments with greater than 99.99% 16O-16O, O2 gas  (Assonov & Brenninkmeijer, 2003a)
        λ = 0.528 # Best estimate for the general use in 17O correction of CO2 originating from the relevant, most abundand, natural oxygen pools, where CO2 rapidly exchanges with liquid water (Assonov & Brenninkmeijer, 2003b). 
        R18_VSMOW = 0.00208835 # Ratio of 18O-to-16O of the VSMOW standard, but converted to the VPDB scale, taking into account the conversion of calcite to CO2 at 25⁰C (As referenced in Brand et al., 2010)

        R13_list = []

        for i in range(len(δ13C_values)):
            R13_values = ((δ13C_values[i] / 1000) + 1) * R13_VPDB
            R13_list.append(R13_values)

        R18_list = []

        for i in range(len(δ18O_values)):
            R18_values = ((δ18O_values[i] / 1000) + 1) * R18_VSMOW
            R18_list.append(R18_values)
            
        R17_list = []

        for i in range(len(R18_list)):
            R17_values = ((R18_list[i] / R18_VSMOW) ** λ) * R17_VSMOW
            R17_list.append(R17_values)

        # Based on the R13, R18 and R17 data, the normalized concentrations for each isotopes to the respective element concentration can be estimated. That is 12C, 13C, 16O, 17O, and 18O. The values and letters had to be switched because 12C, for example, can produce syntax errors when calculations are being performed.

        C12_list = []

        for i in range(len(R13_list)):
            C12_value = 1 / (1 + R13_list[i])
            C12_list.append(C12_value)

        C13_list = []

        for i in range(min([len(R13_list), len(C12_list)])):
            C13_value = R13_list[i] / (1 + R13_list[i])
            C13_list.append(C13_value)
    
        O16_list = []

        for i in range(min([len(R17_list), len(R18_list)])):
            O16_value = 1 / (1 + R17_list[i] + R18_list[i])
            O16_list.append(O16_value)

        O17_list = []

        for i in range(len(R17_list)):
            O17_value = R17_list[i] / (1 + R17_list[i] + R18_list[i])
            O17_list.append(O17_value)

        O18_list = []

        for i in range(len(R18_list)):
            O18_value = R18_list[i] / (1 + R17_list[i] + R18_list[i])
            O18_list.append(O18_value)

        # In this section, the expected concentrations for each isotopologue are estimated probabilistically. These are calculated for R45, R46 and R47, based on the relative concentration data from the previous section.

        R45_sto_list = []

        for i in range(min([len(C12_list), len(C13_list), len(O16_list), len(O17_list)])):
            
            # A mass 45 CO2 isotopologue can result from the isotope combination 16O-13C-16O or 17O-12C-16O. These isotopologues can be determined by multiplying the corresponding values for each isotope from the previous code. The multiplication can be performed since all isotope values are normalized. Hence, the results from the multiplications of these isotope concentrations are somewhere in the range between zero and one, and are stochastic or probabilistic.
            
            # The '2' multiplying 17O-12C-16O arises from the fact that the 17O in a carbonate group can be in one of two possible positions when CO2 is formed from acid digestion. By 'positions' it is implied that a 16O was released during acid digestion, leaving two posible potential sites in the carbonate group of a crystal lattice. A factor '2' multiplication is found multiplying 16O-12C-18O and 17O-13C-16O in R46, 16O-13C-18O and 17O-12C-18O in R47, and 17O-13C-18O in R48. See appendix for further details.   
            
            R45_sto_value = ((C13_list[i] * O16_list[i] ** 2) + (2 * C12_list[i] * O16_list[i] * O17_list[i])) /(C12_list[i] * O16_list[i] ** 2)
            R45_sto_list.append(R45_sto_value)

        R46_sto_list = []

        for i in range(min([len(C12_list), len(C13_list), len(O16_list), len(O17_list), len(O18_list)])):
            
            #  A mass 46 CO2 isotopologue can result from the combinations: 16O-12C-18O, 17O-12C-17O or 17O-13C-16O.
            
            R46_sto_value = ((2 * C12_list[i] * O16_list[i] * O18_list[i]) + (C12_list[i] * O17_list[i] ** 2) + (2 * C13_list[i] * O16_list[i] * O17_list[i])) / (C12_list[i] * O16_list[i] ** 2)
            R46_sto_list.append(R46_sto_value)

        R47_sto_list = []

        for i in range(min([len(C12_list), len(C13_list), len(O16_list), len(O17_list), len(O18_list)])):
            
            # A mass 47 CO2 isotopologue can result from the combinations: 16O-13C-18O, 17O-13C-17O or 17O-12C-18O.
            
            R47_sto_value = ((2 * C13_list[i] * O16_list[i] * O18_list[i]) + (C13_list[i] * O17_list[i] ** 2) + (2 * C12_list[i] * O17_list[i] * O18_list[i])) / (C12_list[i] * O16_list[i] ** 2)
            R47_sto_list.append(R47_sto_value)

        R48_sto_list = []

        for i in range(min([len(C12_list), len(C13_list), len(O16_list), len(O17_list)])):
            
            # A mass 48 CO2 isotopologue can result from the combinations: 18O-12C-18O or 17O-13C-18O.
            
            R48_sto_value = ((C12_list[i] * O18_list[i] ** 2) + (2 * C13_list[i] * O17_list[i] * O18_list[i]))/(C12_list[i] * O16_list[i] ** 2)
            R48_sto_list.append(R48_sto_value)

        R49_sto_list = []

        for i in range(min([len(C12_list), len(C13_list), len(O16_list)])):
            
            # A mass 49 CO2 isotopologue can result from the combination 18O-13C-18O.
            
            R49_sto_value = (C13_list[i] * O18_list[i] ** 2) / (C12_list[i] * O16_list[i] ** 2)
            R49_sto_list.append(R49_sto_value)

        # The stochastic R45 and R46 values are redefined in terms of the δ definition (Hoefs, 2018). These stochastic values are the same to those measured directly by mass spectrometry. [0], [2] and [3] are the working gas, endmember 1 and endmember 2 values from the according R45 and R46 stochastic value lists.
        
        δ45_endmember_1 = ((R45_sto_list[2] / R45_sto_list[0]) - 1) * 1000

        δ45_endmember_2 = ((R45_sto_list[3] / R45_sto_list[0]) - 1) * 1000

        δ46_endmember_1 = ((R46_sto_list[2] / R46_sto_list[0]) - 1) * 1000

        δ46_endmember_2 = ((R46_sto_list[3] / R46_sto_list[0]) - 1) * 1000
        
        # The variable "Δ47_acid_α" is the acid correction for Δ47 during carbonate's digestion which is temperature dependant (Petersen et al., 2019). Δ47_acid_α is defined by the difference between a given carbonate sample not loosing any of its 13C–18O bonds assuming a carbonate group of interest is the isotopologue 13C-18O-16O2.   
        
        Δ47_acid_α = (0.0383 * (10**6 / ((273.15 + temperature) ** 2))) + 0.258
        
        # As the "Δ47_acid_α" variable, "Δ48_acid_α" is the acid correction for mass 48 isotopologues. Δ48_acid_α is assumed to be 0.138, as determined at 90 °C by Bajnai et al. (2020). However, once a relationship is established, the constant can be replaced by an equation.  
        
        Δ48_acid_α = 0.138
        
        def remove_arf (endmemberΔi_value, Δ = 47):
    
            """ This function undoes the value transformations to the absolute reference frame (remove_arf) as defined by Dennis et al., (2011). The afr function corrects the value for the acid fractionation and accounts for mass spectrometric artifacts by applying the empirical transfer function (Dennis et al. (2011)). By default, the acid fractionation correction, Δ = 47, is required when Δ47 calculations are performed. Although fractionation does occur when calculating Δ48 or Δ49, Δ = 48, we assume it is neglegible. """
    
            if Δ == 47:
                
                # This step removes the acid fractionation correction, Δ47_acid_α. 
                
                ΔiRF = endmemberΔi_value - Δ47_acid_α
                
                # The corrected value, that is ΔiRF, is standardized to the absolute reference frame (Dennis et al., 2011), which is defined by the empirical transfer function of the instrument used. The empirical transfer function is defined by measuring  This stndardization has to be undone
                
                ΔiSGvsWG0 = (ΔiRF - empirical_transfer_function_intercept_Δ47) / empirical_transfer_function_slope_Δ47

            elif Δ == 48:
                
                # These steps are the absolute reference frame (Dennis et al., 2011) equivalent steps for Δ48 measurements.
                
                ΔiRF = endmemberΔi_value - Δ48_acid_α
                
                ΔiSGvsWG0 = (ΔiRF - empirical_transfer_function_intercept_Δ48) / empirical_transfer_function_slope_Δ48
                
            return ΔiSGvsWG0

        def δi(ΔiSGvsWG0, sto_Ri, WG_Ri, standard_gas_line_slope = 47):
    
            """ The δi function which correlates the Δ47 or Δ48 sample measurements, to a temperature-Δ47 or Δ48 calibration based on standards. These standards are heated CO2, or CO2 gas from carbonate standards. This step allow to determine the final clumped measurement. The function calculates δ47 or δ48 from the endmember "raw_Ri" which has had an acid correction and an empirical transfer function correction performed by the previous step "remove_arf". The function builds on the R47 and R48 stochastic values for a given endmember, sto_Ri. sto_Ri corresponds to endmember 1 (R47_sto_list[2] or R48_sto_list[2]) and endmember 2 (R47_sto_list[3] or R48_sto_list[3]) values with the respective working gas value, WG_Ri (R47_sto_list[0] or R48_sto_list[0]). The standard_gas_line_slope_D47 is a default value defined by the user at the begining of the script. However, when the standard_gas_line_slope = 48, the standard_gas_line_slope_D48 value will be used instead. """
    
            if standard_gas_line_slope == 47:
               
               m =  ((sto_Ri) / (WG_Ri - heated_gas_line_slope_Δ47 * sto_Ri))
               
               n = (1000 * (sto_Ri - WG_Ri) / (WG_Ri - heated_gas_line_slope_Δ47 * sto_Ri))
               
               δ47 = m * ΔiSGvsWG0 + n
               
               return δ47
    
            elif standard_gas_line_slope == 48:
               
               m =  ((sto_Ri) / (WG_Ri - heated_gas_line_slope_Δ48 * sto_Ri))
               
               n = (1000 * (sto_Ri - WG_Ri) / (WG_Ri - heated_gas_line_slope_Δ48 * sto_Ri))
               
               δ48 = m * ΔiSGvsWG0 + n
               
               return δ48 

        # δ47 and δ48 are the endmember values once the acid fractionation has been applied and then ploted on the absolute reference frame. This allows the user to mix sample values which have been performed on different mass spectrometers, or different calibrations of the empirical transfer function.

        δ47_endmember_1 = δi(remove_arf(endmember_1_Δ47), R47_sto_list[2], R47_sto_list[0])
        
        δ47_endmember_2 = δi(remove_arf(endmember_2_Δ47), R47_sto_list[3], R47_sto_list[0])
        
        δ48_endmember_1 = δi(remove_arf(endmember_1_Δ48, Δ = 48), R48_sto_list[2], R48_sto_list[0], standard_gas_line_slope = 48) 
        
        δ48_endmember_2 = δi(remove_arf(endmember_2_Δ48, Δ = 48), R48_sto_list[3], R48_sto_list[0], standard_gas_line_slope = 48) 

        # These mixture values corresponding to the non-linear mixing model and are computated by a weighted average.

        δ45_mixture = mixture(δ45_endmember_1, δ45_endmember_2, endmember_1_contribution)
        
        δ46_mixture = mixture(δ46_endmember_1, δ46_endmember_2, endmember_1_contribution)
        
        δ47_mixture = mixture(δ47_endmember_1, δ47_endmember_2, endmember_1_contribution)
        
        δ48_mixture = mixture(δ48_endmember_1, δ48_endmember_2, endmember_1_contribution)
    
        def Ri(δi_mixture, Ri_sto_list):
            
            # This function is a rearragement of the δ definition (Hoefs, 2018), which solves for the sample ratio (R). Therefore, the Ri function takes the δ45, δ46, δ47, and δ48 mixture values and converts them to their respective R45, R46, R47 and R48 mixture values; as suggested by Affek & Eiler (2006), Appendix B. This step achives the cancelation of the working gas.
            
            Ri_mixture = ((δi_mixture / 1000) + 1) * Ri_sto_list[0]
            return Ri_mixture
        
        R45_mixture = Ri(δ45_mixture, R45_sto_list)
        
        R46_mixture = Ri(δ46_mixture, R46_sto_list)
        
        R47_mixture = Ri(δ47_mixture, R47_sto_list)
        
        R48_mixture = Ri(δ48_mixture, R48_sto_list)
    
        def Δi(Ri_mixture, Ri_sto_list):
            
            # This function takes the corresponding ratio of a mixture (R45_mixture, R46_mixture, R47_mixture or R48_mixture) and defines it according to the δ definition (Hoefs, 2018), but the sample is normalized to the respective mixed value, Ri_sto_list[1].
            
            Δi_mixture = ((Ri_mixture / Ri_sto_list[1]) - 1)
            return Δi_mixture
    
        Δ45_mixing = Δi(R45_mixture, R45_sto_list)
        
        Δ46_mixing = Δi(R46_mixture, R46_sto_list)
        
        Δ47_mixing = Δi(R47_mixture, R47_sto_list)
        
        Δ48_mixing = Δi(R48_mixture, R48_sto_list)
    
        # The Δ47_mixing value, as defined in equation (3) of Huntington et al. (2009).
        
        Δ47_no_arf = (Δ47_mixing - Δ46_mixing - Δ45_mixing) * 1000

        # The Δ48_mixing value, as defined in equation (4) of Huntington et al. (2009).
        
        Δ48_no_arf = (Δ48_mixing - 2 * Δ46_mixing) * 1000

        def add_arf(Δ47 = 47):
    
            """If D47 = True, then the according D47 parameters are passed. If D47 = False, then the according D48 parameters are passed. The acid fractionation associated with Δ48 is assumed to be negligible."""
    
            if Δ47 == 47: # This is the set default value
            
                Δ47SGvsWG0 = Δ47_no_arf - heated_gas_line_slope_Δ47 * δ47_mixture
                
                Δ47RF =  empirical_transfer_function_slope_Δ47 * (Δ47SGvsWG0) + empirical_transfer_function_intercept_Δ47
                
                Δ47AC = Δ47RF + Δ47_acid_α
                
                return Δ47AC
    
            elif Δ47 == 48:
                
               Δ48SGvsWG0 = Δ48_no_arf - heated_gas_line_slope_Δ48 * δ48_mixture
               
               Δ48RF = (Δ48SGvsWG0) * empirical_transfer_function_slope_Δ48 + empirical_transfer_function_intercept_Δ48
               
               Δ48AC = Δ48RF + Δ48_acid_α
               
               return Δ48AC
        
        # The export line is a dictionary, where the "keys" or names in quatation marks are what the number represents at after the colon. The export line will be taking into a DataFrame, which will then be converted to an Excel file. 
        
        export_line = {'endmember 1 Δ47 (CDES)': endmember_1_Δ47, 'endmember 2 Δ47 (CDES)': endmember_2_Δ47,
                       'endmember 1 Δ48 (CDES)': endmember_1_Δ48, 'endmember 2 Δ48 (CDES)': endmember_2_Δ48,
                       
                       'endmember 1 contribution': endmember_1_contribution,'endmember 2 contribution': 1 - endmember_1_contribution,
                       
                       'endmember 1 δ13C (VPDB)': round(endmember_1_δ13C_VPDB, 2), 'endmember 2 δ13C (VPDB)': round(endmember_2_δ13C_VPDB, 2), 
                       
                       'endmember 1 δ18O (VPDB)':  round(endmember_1_δ18O_VPDB, 2), 'endmember 2 δ18O (VPDB)': round(endmember_2_δ18O_VPDB, 2),
                       'endmember 1 δ18O (VSMOW)': round(endmember_1_δ18O_CO2VSMOW, 2), 'endmember 2 δ18O (VSMOW)': round(endmember_2_δ18O_CO2VSMOW, 2), 
                       
                       'δ13C (VPDB) mix': round(mixture_δ13C_VPDB, 2), 'δ18O (VPDB) mix': round(mixture_δ18O_VPDB, 2),
                       'Δ47 model (CDES)': add_arf(Δ47 = 47), 'Δ48 model (CDES)': add_arf(Δ47 = 48),
                       
                       'Γ47': add_arf(Δ47 = 47) - mixture(endmember_1_Δ47, endmember_2_Δ47, endmember_1_contribution),
                       'Γ48': add_arf(Δ47 = 48) - mixture(endmember_1_Δ48, endmember_2_Δ48, endmember_1_contribution)
                       }
    
        return export_line
    
    return mixing_model_calculation(endmember_1_Δ47, endmember_1_Δ48, endmember_1_δ18O_VPDB, endmember_1_δ13C_VPDB, endmember_1_contribution, endmember_2_Δ47, endmember_2_Δ48, endmember_2_δ18O_VPDB, endmember_2_δ13C_VPDB, temperature, empirical_transfer_function_slope_Δ47, empirical_transfer_function_intercept_Δ47, empirical_transfer_function_slope_Δ48, empirical_transfer_function_intercept_Δ48, heated_gas_line_slope_Δ47, heated_gas_line_slope_Δ48, working_gas_δ13C_VPDB, working_gas_δ18O_VSMOW) # Delete once paper analysis is ready

# Remove all "#" in the lines with arrows (<-) below if you rather work on an excel sheet.

    #def mixing_model_data_simulation(): <-
        
        # Import pandas module 
        
        #import pandas as pd <-
        
        #df_list = [] <-
        
        #for i in range(len(endmember_1_contribution)): <-
                 
            #loop = mixing_model_calculation(endmember_1_Δ47, endmember_1_Δ48, endmember_1_δ18O_VPDB, endmember_1_δ13C_VPDB, endmember_1_contribution[i], endmember_2_Δ47, endmember_2_Δ48, endmember_2_δ18O_VPDB, endmember_2_δ13C_VPDB, temperature, empirical_transfer_function_slope_Δ47, empirical_transfer_function_intercept_Δ47, empirical_transfer_function_slope_Δ48, empirical_transfer_function_intercept_Δ48, heated_gas_line_slope_Δ47, heated_gas_line_slope_Δ48, working_gas_δ13C_VPDB, working_gas_δ18O_VSMOW) <-
                 
            #df_list.append(loop) <-
        
         # Convert list of dictionaries to a pandas dataframe
        #df = pd.DataFrame(df_list) <-
    
        # Export DataFrame as an excel spreadsheet
        #return df.to_excel("output.xlsx") <-
    #return mixing_model_data_simulation() <-

#References

# Affek, & Eiler, J. M. (2006). Abundance of mass 47 CO2 in urban air, car exhaust, and human breath. Geochimica et Cosmochimica Acta, 70(1), 1–12. https://doi.org/10.1016/j.gca.2005.08.021

# Assonov, S. S., & Brenninkmeijer, C. A. M. (2003a). A redetermination of absolute values for 17RVPDB-CO2 and 17RVSMOW. Rapid Communications in Mass Spectrometry, 17(10), 1017–1029. https://doi.org/10.1002/rcm.1011

# Assonov, S. S., & Brenninkmeijer, C. A. M. (2003b). On the 17O correction for CO2 mass spectrometric isotopic analysis. Rapid Communications in Mass Spectrometry, 17(10), 1007–1016. https://doi.org/10.1002/rcm.1012

# Bajnai, D., Guo, W., Spötl, C., Coplen, T. B., Methner, K., Löffler, N., Krsnik, E., Gischler, E., Hansen, M., Henkel, D., Price, G. D., Raddatz, J., Scholz, D., & Fiebig, J. (2020). Dual clumped isotope thermometry resolves kinetic biases in carbonate formation temperatures. Nature Communications, 11(1), 4005–4005. https://doi.org/10.1038/s41467-020-17501-0

# Brand, W. A., Assonov, S. S., & Coplen, T. B. (2010). Correction for the 17O interference in δ(13C) measurements when analyzing CO2 with stable isotope mass spectrometry (IUPAC Technical Report). Pure and Applied Chemistry, 82(8), 1719–1733. https://doi.org/10.1351/PAC-REP-09-01-05

# Defliese, W. F., Hren, M. T., & Lohmann, K. C. (2015). Compositional and temperature effects of phosphoric acid fractionation on Δ47 analysis and implications for discrepant calibrations. Chemical Geology, 396, 51–60. https://doi.org/10.1016/j.chemgeo.2014.12.018

# Defliese, W. F., Lohmann, K. C., (2015). Non-linear mixing effects on mass-47 CO2clumped isotope thermometry: Patterns and implications. Rapid Communications in Mass Spectrometry 29, 901–909.. doi:10.1002/rcm.7175

# Dennis, K. J., Affek, H. P., Passey, B. H., Schrag, D. P., & Eiler, J. M. (2011). Defining an absolute reference frame for 'clumped' isotope studies of CO2. Geochimica et Cosmochimica Acta, 75(22), 7117–7131. https://doi.org/10.1016/j.gca.2011.09.025

# Hoefs, J. (2018). Stable isotope geochemistry (Eighth edition ed.). Springer.

# Huntington, K. W., Eiler, J. M., Affek, H. P., Guo, W., Bonifacie, M., Yeung, L. Y., Thiagarajan, N., Passey, B., Tripati, A., Daëron, M., & Came, R. (2009). Methods and limitations of ‘clumped’ CO2 isotope (Δ47) analysis by gas-source isotope ratio mass spectrometry. Journal of Mass Spectrometry, 44(9), 1318-1329. 

# Kim, S.-T., Coplen, T. B., & Horita, J. (2015). Normalization of stable isotope data for carbonate minerals; implementation of IUPAC guidelines. Geochimica et Cosmochimica Acta, 158, 276–289. https://doi.org/10.1016/j.gca.2015.02.011

# Petersen, S. V., Defliese, W. F., Saenger, C., Daëron, M., Huntington, K. W., John, C. M., et al (2019). Effects of Improved 17O Correction on Interlaboratory Agreement in Clumped Isotope Calibrations, Estimates of Mineral‐Specific Offsets, and Temperature Dependence of Acid Digestion Fractionation. Geochemistry, geophysics, geosystems : G3, 20(7), 3495-3519. https://doi.org/10.1029/2018GC008127

# Swart, P. K., Burns, S., & Leder, J. (1991). Fractionation of the stable isotops of oxygen and carbon in carbon dioxide during the reaction of calcite with phosphoric acid as a function of temperature and technique. Chemical Geology, 86(2), 89–96. https://doi.org/10.1016/0168-9622(91)90055-2

# Zhang, Ql, Chang, Tl, & Li, Wj. (1990). A CALIBRATED MEASUREMENT OF THE ATOMIC-WEIGHT OF CARBON. Science Bulletin., 35(4), 290–296.
