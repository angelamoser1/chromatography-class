# -*- coding: utf-8 -*-
"""
This script contains the function that pulls in user entered inputs and sets
other model configuration parameters that are needed for the CADET simulation.

The number of components is fixed at three proteins and the residence time is 
fixed at a single value for all steps of the process.
"""
from pathlib import Path
import math
import pandas as pd


def create_ms():    
    # update path to CADET bin folder
    cadet_bin_path = Path() / r'/Users/angelamoser/opt/anaconda3/envs/CADET-env/bin'
    
    # pull in experiment specific information
    # read inputs from spreadsheet ###########################################
    run_conditions_file = Path() / r'/Users/angelamoser/Documents/05_Academic/RPI/07_Chromatography_Class/modeling_assignment_code/run_conditions.xlsx'
    df = pd.read_excel(run_conditions_file, header=None, usecols=[0, 1])
    run_conditions = df.dropna()
    run_conditions = run_conditions.T # transpose
    run_conditions.columns = run_conditions.iloc[0]  # Set the first row as the header
    
    # Step 4: Remove the first row from the DataFrame
    run_conditions = run_conditions.drop(run_conditions.index[0])
    
    column = 1
       
    load_salt = run_conditions['load salt [mM]'].at[column]
    wash_salt = run_conditions['wash salt [mM]'].at[column]
    step1_salt = run_conditions['step 1 salt [mM]'].at[column]
    step2_salt = run_conditions['step 2 salt [mM]'].at[column]
    step3_salt = run_conditions['step 3 salt [mM]'].at[column]
    
    load_CV = run_conditions['load CVs'].at[column]
    wash_CV = run_conditions['wash CVs'].at[column]
    step1_CV = run_conditions['step 1 CVs'].at[column]
    step2_CV = run_conditions['step 2 CVs'].at[column]
    step3_CV = run_conditions['step 3 CVs'].at[column]
    
    residence_time = run_conditions['residence time [min]'].at[column]
            
    comp1_conc = run_conditions['impurity A concentration [mg/mL]'].at[column]
    comp2_conc = run_conditions['product concentration [mg/mL]'].at[column]
    comp3_conc = run_conditions['impurity B concentration [mg/mL]'].at[column]
       
    feed_concs = [comp1_conc, comp2_conc, comp3_conc]
    
    MW_1 = run_conditions['impurity A MW'].at[column]
    MW_2 = run_conditions['product MW'].at[column]
    MW_3 = run_conditions['impurity B MW'].at[column]
    
    protein_MWs = [MW_1, MW_2, MW_3] # kDa
     
    # column information
    col_id = run_conditions['column ID [cm]'].at[column]  # [cm], Column diameter
    col_length = run_conditions['column length [cm]'].at[column]  # [cm], Column length   
    Ee = run_conditions['Ee'].at[column]
    Ep = run_conditions['Ep'].at[column]
    particle_diameter = run_conditions['particle diameter [um]'].at[column]
    pore_radius = run_conditions['pore radius [nm]'].at[column]
    ionic_cap_mM = run_conditions['ionic capacity [mol / L column]'].at[column] # [mol / L column]
    
    Keq_1 = run_conditions['Keq 1'].at[column]
    steric_factor_1 = run_conditions['steric factor 1'].at[column]
    characteristic_charge_1 = run_conditions['characteristic charge 1'].at[column]
    
    Keq_2 = run_conditions['Keq 2'].at[column]
    steric_factor_2 = run_conditions['steric factor 2'].at[column]
    characteristic_charge_2 = run_conditions['characteristic charge 2'].at[column]
    
    Keq_3 = run_conditions['Keq 3'].at[column]
    steric_factor_3 = run_conditions['steric factor 3'].at[column]
    characteristic_charge_3 = run_conditions['characteristic charge 3'].at[column]
    
    Keqs = [Keq_1, Keq_2, Keq_3]
    steric_factors = [steric_factor_1, steric_factor_2, steric_factor_3]
    characteristic_charges = [characteristic_charge_1, characteristic_charge_2, characteristic_charge_3]
    
    # Elution pooling criteria
    cutoff = run_conditions['collection cutoff [mg/mL]'].at[column] # [mg/ml], Start pool concentration
    
    # choose the x axis for output plots
    # xaxis_plot_set = run_conditions['x axis'].at[column]
    xaxis_plot_set = 'CV'

    # Computational inputs
    n_threads = 4   # [-], Number of CPU threads (max is 2*CPU cores), increases simulation speed
    abstol = 1e-6   # usually 1e-6. absolute tolerance for error during time integration
    reltol = 1e-6   # usually 1e-6. relative tolerance for error during time integration
    axial_nodes = 50    # [1/cm], Number of axial discretization nodes - # PER 2.5 cm OF COLUMN!!
    particle_nodes = 25    # [-], Number of particle discretization nodes
                             
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

    ms['n_threads'] = n_threads
    ms['abstol'] = abstol
    ms['reltol'] = reltol
    ms['axial_nodes'] = axial_nodes 
    ms['particle_nodes'] = particle_nodes
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