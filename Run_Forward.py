# -*- coding: utf-8 -*-
'''
This script is for running forward simulations based on parameters 
specified by the user. 
CADET License: https://cadet.github.io/v4.3.0/license.html
'''

import time
import create_ms
import CADET_sim
import evaluate_sim
import plotting
import pandas as pd


def main():
    # start timer
    tic = time.perf_counter()
        
    # pull in inputs
    ms = create_ms.create_ms()
            
    # run simulation with parameters entered by user in data_entry_forward
    # Model structure is sent to the simulation function  
    sim_solution = CADET_sim.run_cadet_sim(ms)

    # convert time from [sec] to [min], create a [CV] column, and convert [H+] to pH
    sim_solution = evaluate_sim.process_simulation_data(sim_solution, ms) 
                  
    # calculates DBC, EBC, pool volume, recovery, etc. and saves results as a dictionary
    results = evaluate_sim.evaluate_simulation(sim_solution, ms)
    
    print('Forward Simulation Results:')    
    # Convert from list of dicts to DataFrame
    results = pd.DataFrame(results)
    # Print results table to console
    with pd.option_context('display.max_rows', 5,
                           'display.max_columns', None,
                           'display.precision', 3,
                           ):
        print(results)
                 
    # Generate plots
    plot = plotting.generate_plot(sim_solution, ms)  
                   
    # end timer and print run time
    toc = time.perf_counter()
    print(f'\n Completed in {toc - tic:.4f} seconds')
    

###############################################################################       
if __name__ == '__main__':
    main()