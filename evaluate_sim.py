# -*- coding: utf-8 -*-
"""
This script contains the functions used to process the simulation results
"""
import numpy as np


def process_simulation_data(sim_solution, ms):    
    time = sim_solution.solution_times
    solution_array = sim_solution.unit_003.solution_outlet
    inlet_array = sim_solution.unit_000.solution_inlet

    inlet = {}
    s = {}
    s['time_sec'] = time # [s]
        
    inlet['prot_c_mM'] = []
    s['prot_c_mM'] = []
    inlet['prot_c_mg'] = []
    s['prot_c_mg'] = []
    
    n_comp = len(solution_array[1,:])-1
    s['salt'] = solution_array[:,0] # [mM] 
               
    for c in range(n_comp):
        idx = c+1

        inlet_component_c = [inlet_array[i][idx] for i in range(len(inlet_array))]    # [mM]           
        component_c = [solution_array[i][idx] for i in range(len(solution_array))]    # [mM]

        # Convert concentration from mM to mg/ml
        inlet_component_c_mg = [mM*ms['MW'][c] for mM in inlet_component_c]
        component_c_mg = [mM*ms['MW'][c] for mM in component_c]
        
        inlet['prot_c_mM'].append(inlet_component_c)
        s['prot_c_mM'].append(component_c)
        inlet['prot_c_mg'].append(inlet_component_c_mg)
        s['prot_c_mg'].append(component_c_mg)
    
    # Convert from time (s) to time (min)
    s['time_min'] = s['time_sec']/60
    
    # Create a volume column [mL] from time column [min] using flow rate [m^3/s]
    s['volume'] = s['time_sec'] * ms['flow_rate']*1e6
    
    # Convert from volume [mL] to CV
    s['CV'] = s['volume']/(ms['col_vol']*1e6)
        
    # calculate cumulative mass for all components
    cum_masses = []
    for comp in s['prot_c_mg']:
        cum_mass = [sum(comp[:i]) for i in range(len(comp))]
        cum_masses.append(cum_mass)

    # calculate cumulative mass loaded (from inlet profile) for each component
    in_masses = []
    for comp in inlet['prot_c_mg']:
        in_mass = [sum(comp[:i]) for i in range(len(comp))]
        in_masses.append(in_mass)
        
    total_mass = [sum(j) for j in zip(*cum_masses)]

    purity = []
    cum_yield = []
    for k in range(len(cum_masses[0])):
        if total_mass[k] > 0 and cum_masses[0][k] > 0:
            p = cum_masses[0][k] / total_mass[k]
            y = cum_masses[0][k] / in_masses[0][k]
        else:
            p = 0
            y = 0
        purity.append(p)
        cum_yield.append(y)
             
    s['purity'] = purity
    s['yield'] = cum_yield
    
    return s
 

def same_dimensions(x_data, y_data):
    x_dim = len(np.shape(x_data))
    y_dim = len(np.shape(y_data))
    if x_dim == y_dim:
        return x_data, y_data
    new_x_data = []
    new_y_data = []
    if x_dim == 1 and y_dim ==2:
        for i in y_data:
            new_x_data.append(x_data)
            new_y_data.append(np.array(i))
        np.array(new_x_data)
        np.array(new_y_data)
        return new_x_data, new_y_data
    

def select_section(start_cutoff, end_cutoff, x_data, y_data):
    # select data from a section for a pair of x and y data
    start_idx = (np.abs(x_data - start_cutoff)).argmin()
    end_idx = (np.abs(x_data - end_cutoff)).argmin()
    x_section = x_data[start_idx:end_idx+1]
    y_section = y_data[start_idx:end_idx+1] 
               
    return x_section, y_section


def select_section_multi(start_cutoff, end_cutoff, x_data, y_data):
    # select data from a section for multiple x and y data pairs
    x_section = []
    y_section = []
    for i in range(len(x_data)):
        x_sec, y_sec = select_section(start_cutoff, end_cutoff, x_data[i], y_data[i])
        x_section.append(x_sec)
        y_section.append(y_sec)
    
    return x_section, y_section
          
# s is the processed simulation solution as a dictionary
def evaluate_simulation(s, ms):
    vol_sim, conc_sim = same_dimensions(s['volume'], s['prot_c_mg'])
        
    # Calculate recovery
    load_masses = []
    mass_bal = []
    for idx, comp_feed_c in enumerate(ms['feed_conc']):
        load_mass = comp_feed_c*ms['load_vol']*1e6
        load_masses.append(load_mass)
        total_mass = np.trapz(s['prot_c_mg'][idx], s['volume'])    
        mass_bal_comp = 100*total_mass/load_mass
        mass_bal.append(mass_bal_comp)

    # Find start and ends of each step
    start_vols = [ms['wash_vol']*1e6, ms['step1_vol']*1e6, ms['step2_vol']*1e6]
    end_vols = [ms['step1_vol']*1e6, ms['step2_vol']*1e6, ms['step3_vol']*1e6]
    
    pool_vols = []
    pool_concs = []
    
    # calculate pool volume and concentration for each component for each step
    for i in range(len(start_vols)):
        start_vol = start_vols[i]
        end_vol = end_vols[i]
        cutoff = ms['cutoff']
        # Select data from the desired section
        elution_vol, elution_conc = select_section_multi(start_vol, end_vol, vol_sim, conc_sim)
        # check if the cutoff was exceeded in the elution section
        main_comp_c = elution_conc[1]
        main_vol = elution_vol[1]

        if max(main_comp_c) >= cutoff:        
            # Find pool indices
            index_max = np.argmax(main_comp_c)
            pool_start_idx = (np.abs(np.array(main_comp_c[:index_max]) - cutoff)).argmin()
            pool_end_idx = (np.abs(np.array(main_comp_c[index_max:]) - cutoff)).argmin() + index_max
            # Calculate pool volume [mL]
            pool_volume = main_vol[pool_end_idx] - main_vol[pool_start_idx]
            pool_conc = []
            for idx in range(len(ms['feed_conc'])):
                pool_mass = np.trapz(elution_conc[idx][pool_start_idx:pool_end_idx], elution_vol[idx][pool_start_idx:pool_end_idx])
                # Calculate elution pool concentration
                pool_conc_comp = pool_mass/pool_volume
                pool_conc.append(pool_conc_comp)
                
            # Convert pool volume to CV
            pool_volume = [pool_volume/(ms['col_vol']*1e6)] * len(ms['feed_conc'])
        else:
            pool_volume = [0] * len(ms['feed_conc'])
            pool_conc = [0] * len(ms['feed_conc'])
            
        # add result for each step in the method to the list
        pool_vols.append(pool_volume)
        pool_concs.append(pool_conc)
            
    # calculate yields and purities
    yields = []
    purities = []
    for idx in range(len(pool_vols)):
        vol = [v*ms['col_vol']*1e6 for v in pool_vols[idx]] # convert from CV to mL
        conc = pool_concs[idx]
        step_masses = [vol[j]*conc[j] for j in range(len(vol))]
        step_yields = [step_masses[k] * 100 / load_masses[k] for k in range(len(step_masses))]
        yields.append(step_yields)
        try:
            step_purities = [comp_mass * 100 / sum(step_masses) for comp_mass in step_masses]
            purities.append(step_purities)
        except ZeroDivisionError:
            purities.append(['N/A', 'N/A', 'N/A'])
            
               
    # save results in dictionary   
    results = {
               'pool vol 1':pool_vols[0],
               'pool vol 2':pool_vols[1],
               'pool vol 3':pool_vols[2],
               'pool C 1':pool_concs[0], 
               'pool C 2':pool_concs[1],
               'pool C 3':pool_concs[2],
               'yield 1':yields[0],
               'yield 2':yields[1],
               'yield 3':yields[2],
               'purity 1':purities[0],
               'purity 2':purities[1],
               'purity 3':purities[2],
               'Keq':ms['Keq'],
               'sigma':ms['steric_factor'],
               'nu':ms['characteristic_charge'], 
               'Dp':ms['pore_diff'],
               'mass_bal':mass_bal 
                }

    return results