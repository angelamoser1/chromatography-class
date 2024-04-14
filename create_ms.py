# -*- coding: utf-8 -*-
"""
This script contains the function that pulls in user entered inputs and sets
other model configuration parameters that are needed for the CADET simulation.

The number of components is fixed at three proteins and the residence time is 
fixed at a single value for all steps of the process.
"""
import math


def create_ms(run_conditions):    
    # update path to CADET bin folder
    cadet_bin_path = run_conditions['cadet_bin_path']
           
    load_salt = run_conditions['load salt [mM]']
    wash_salt = run_conditions['wash salt [mM]']
    step1_salt = run_conditions['step 1 salt [mM]']
    step2_salt = run_conditions['step 2 salt [mM]']
    step3_salt = run_conditions['step 3 salt [mM]']
    
    load_CV = run_conditions['load CVs']
    wash_CV = run_conditions['wash CVs']
    step1_CV = run_conditions['step 1 CVs']
    step2_CV = run_conditions['step 2 CVs']
    step3_CV = run_conditions['step 3 CVs']
    
    residence_time = run_conditions['residence time [min]']
            
    comp1_conc = run_conditions['impurity A concentration [mg/mL]']
    comp2_conc = run_conditions['product concentration [mg/mL]']
    comp3_conc = run_conditions['impurity B concentration [mg/mL]']
       
    feed_concs = [comp1_conc, comp2_conc, comp3_conc]
    
    MW_1 = run_conditions['impurity A MW']
    MW_2 = run_conditions['product MW']
    MW_3 = run_conditions['impurity B MW']
    
    protein_MWs = [MW_1, MW_2, MW_3] # kDa
     
    # column information
    col_id = run_conditions['column ID [cm]']  # [cm], Column diameter
    col_length = run_conditions['column length [cm]']  # [cm], Column length   
    Ee = run_conditions['Ee']
    Ep = run_conditions['Ep']
    particle_diameter = run_conditions['particle diameter [um]']
    pore_radius = run_conditions['pore radius [nm]']
    ionic_cap_mM = run_conditions['ionic capacity [mol / L column]'] # [mol / L column]
    
    Keq_1 = run_conditions['Keq 1']
    steric_factor_1 = run_conditions['steric factor 1']
    characteristic_charge_1 = run_conditions['characteristic charge 1']
    
    Keq_2 = run_conditions['Keq 2']
    steric_factor_2 = run_conditions['steric factor 2']
    characteristic_charge_2 = run_conditions['characteristic charge 2']
    
    Keq_3 = run_conditions['Keq 3']
    steric_factor_3 = run_conditions['steric factor 3']
    characteristic_charge_3 = run_conditions['characteristic charge 3']
    
    Keqs = [Keq_1, Keq_2, Keq_3]
    steric_factors = [steric_factor_1, steric_factor_2, steric_factor_3]
    characteristic_charges = [characteristic_charge_1, characteristic_charge_2, characteristic_charge_3]
    
    # Elution pooling criteria
    cutoff = run_conditions['collection cutoff [mg/mL]'] # [mg/ml], Start pool concentration
    
    # choose the x axis for output plots
    # xaxis_plot_set = run_conditions['x axis']
    xaxis_plot_set = 'CV'
                             
    # [1/s], Adsorption rate constant
    k_kin = [1e8] * (len(protein_MWs)) 

###############################################################################    
###############################################################################
    
    # Constants for diffusivity calculation   
    Kb              = 1.38064e-23                       # [m^2*kg/(s^2*K)], Boltzmann constant
    T               = 297                               # [K], Temperature
    visc            = 0.001                             # [kg/(m*s)], Viscosity                 
    rho             = 1000                              # [kg/m^3], Density

    # UNIT CONVERSIONS (all to SI - m, mol, s) - VERY IMPORTANT FOR CADET!!
    col_id          = col_id/100                        # [m]
    col_length      = col_length/100                    # [m]
    particle_diameter = particle_diameter*1e-6          # [m]
    r_pore = pore_radius*1e-9                           # [m]

    # More calculations    
    Ac_col          = 0.25*math.pi*col_id**2
    col_lin_vel = col_length/(residence_time*60) # [m/s]
    int_vel = col_lin_vel/Ee # [m/s]            
    flow_rate = Ac_col*col_lin_vel # [m^3/s]                                    
    col_vol = Ac_col*col_length # [m^3]          
    Et = Ee + Ep*(1-Ee) # [-], Total porosity
    ionic_capacity = ionic_cap_mM / (1 - Et)  # [mol / m^3 solid phase]
    HUV = col_vol*Et # [m^3]
    
###############################################################################
                    
    # Combine all molecular weights in list for further calculations
    MW = protein_MWs
    feed_conc_mM = [(feed_concs[i] / MW[i]) for i in range(len(feed_concs))]

    R_Na = 1.94e-10 # m
    # Calculate hydrodynamic radius from globular ptotein correlation [nm] 
    R_hyd = [R_Na] + [1e-9*0.7429*i**0.3599 for i in MW]

    # Calculate molecular diffusivities [m^2/s], 
    Do = [Kb*T/(6*math.pi*visc*i) for i in R_hyd]   

    #### Transport parameter calculations #####################################       
    # Calculate one Reynolds number for each section
    Re = rho*col_lin_vel*particle_diameter/visc
    # Calculate one Scanlon number for each component
    Sc = [visc/(rho*i) for i in Do]

    # Calculate one Sherwood number for each component
    # (same for each section because flow rates are the same)
    # Sherwood number, mass transfer - Wilson & Geankopolis 1966
    Sh = [(1.09*Re**0.33*i**0.33)/Ee for i in Sc]
    
    # Calculate one film mass transfer coefficient for each component [m/s]
    # (same for each section because flow rates are the same)
    k_film = [(Sh[i]*Do[i]/particle_diameter) for i in range(len(Do))]
           
    # Column peclet number, dispersion - Rastegar & Gu 2017
    # (same for each section because flow rates are the same)
    Pe = (0.7*min(Do)/(col_length*col_lin_vel) +
         (particle_diameter/col_length)*(Ee/(0.18 + 0.008*Re**0.59)))**-1
    
    # Column axial dispersion [m^2/s]
    # (same for each section because flow rates are the same)
    Dax_col = col_lin_vel*col_length/Pe
    
    # Effective pore diffusivity (for salt), Mackie-Meares correlation
    salt_pore_diff = (Ep/(2-Ep))**2*Do[0]
    
    # Pore diffusivity for protein based on tortuosity estimate
    tortuosity = 4 # literally a guess but what can you do ¯\_(ツ)_/¯
    Ep_mAb = [0.75] * len(MW) # estimate based on Ep from mAb3 with Eshmuno A in Pabst paper
    lambda_m = [Rh / r_pore for Rh in R_hyd[1:]]

    # diffusional hinderance coefficient from Jungbauer and Carta
    # psi_p = [1 + 1.125 * l * math.log(l) - 1.539 * l for l in lambda_m]
    
    # diffusional hinderance coefficient correlation from Dechadilok and Deen 2006
    diff_hind = [1 + 1.125 * l * math.log(l)
                 - 1.5604 * l 
                 + 0.528155 * l**2
                 + 1.91521 * l**3
                 - 2.81903 * l**4
                 + 0.270788 * l**5
                 + 1.10115 * l**6
                 - 0.435933 * l**7
                 for l in lambda_m]
    
    pore_diff = [Ep_mAb[i] * diff_hind[i] * Do[i+1] / tortuosity for i in range(len(diff_hind))]

###############################################################################
    # create the model structure (dictionary) to contain data for passing into functions
    ms = {}
    
    # Assign all to model structure (ms)
    ms['cadet_bin_path'] = cadet_bin_path   
    ms['xaxis_plot_set'] = xaxis_plot_set 
    ms['cutoff'] = cutoff
    ms['col_id'] = col_id 
    ms['col_length'] = col_length
    ms['Ac_col'] = Ac_col
    ms['particle_diameter'] = particle_diameter 
    ms['ionic_capacity'] = ionic_capacity
    ms['Dax_col'] = Dax_col
    
    ms['Ee'] = Ee 
    ms['Ep'] = Ep 
    ms['Et'] = Et
    ms['flow_rate'] = flow_rate 
    ms['int_vel'] = int_vel
    ms['k_film'] = k_film 
    ms['col_lin_vel'] = col_lin_vel
    ms['MW'] = MW 
    ms['feed_conc'] = feed_concs
    ms['feed_conc_mM'] = feed_conc_mM
    ms['col_vol'] = col_vol 
    ms['r_hyd'] = R_hyd 
    ms['k_kin'] = k_kin 

    ms['HUV'] = HUV
    
    ms['load_CV'] = load_CV
    ms['wash_CV'] = wash_CV
    ms['step1_CV'] = step1_CV
    ms['step2_CV'] = step2_CV
    ms['step3_CV'] = step3_CV
    
    ms['load_salt'] = load_salt
    ms['wash_salt'] = wash_salt
    ms['step1_salt'] = step1_salt
    ms['step2_salt'] = step2_salt
    ms['step3_salt'] = step3_salt
    
    ms['residence_time'] = residence_time # [min]
    ms['flow_rate'] = flow_rate # [m^3/sec]
    
    # Calculate the end time for each section in the simulation [s]   
    ms['load_time'] = ms['residence_time']*60*ms['load_CV']
    ms['wash_time'] = ms['residence_time']*60*ms['wash_CV'] + ms['load_time']
    ms['step1_time'] = ms['residence_time']*60*ms['step1_CV'] + ms['wash_time']
    ms['step2_time'] = ms['residence_time']*60*ms['step2_CV'] + ms['step1_time']
    ms['step3_time'] = ms['residence_time']*60*ms['step3_CV'] + ms['step2_time']
    
    # Calculate the end volume for each section in the simulation [m^3]
    ms['load_vol'] = ms['col_vol']*ms['load_CV']
    ms['wash_vol'] = ms['col_vol']*ms['wash_CV'] + ms['load_vol']
    ms['step1_vol'] = ms['col_vol']*ms['step1_CV'] + ms['wash_vol']
    ms['step2_vol'] = ms['col_vol']*ms['step2_CV'] + ms['step1_vol']
    ms['step3_vol'] = ms['col_vol']*ms['step3_CV'] + ms['step2_vol']
    
    ms['salt_pore_diff'] = salt_pore_diff
    ms['pore_diff'] = pore_diff # [m^2/s] pore diffusivity
    
    ms['Keq'] = Keqs
    ms['steric_factor'] = steric_factors
    ms['characteristic_charge'] = characteristic_charges
       
    return ms
