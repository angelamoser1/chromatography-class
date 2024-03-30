# -*- coding: utf-8 -*-
"""
This script contains the CADET simulation code for
the 3 component step elution CADET model.

CADET License: https://cadet.github.io/v4.3.0/license.html

This script is the front-end for the CADET simulation software
      
Copyright (C) 2008-2022: The CADET Authors
           Please see the AUTHORS and CONTRIBUTORS file.
  
All rights reserved. This program and the accompanying materials
are made available under the terms of the GNU Public License v3.0 (or, at
your option, any later version) which accompanies this distribution, and
is available at http://www.gnu.org/licenses/gpl.html
"""
import platform
from cadet import Cadet
import numpy as np


def run_cadet_sim(ms):
    
   cadet_bin_path = ms['cadet_bin_path']
   
   if platform.system() == 'Windows':
       cadet_path = cadet_bin_path / 'cadet-cli.exe'
       lwe_path = cadet_bin_path / 'createLWE.exe'
       
   else:
        cadet_path = cadet_bin_path / 'cadet-cli'
        lwe_path = cadet_bin_path / 'createLWE'
     
   if cadet_path.exists() and lwe_path.exists():
       Cadet.cadet_path = cadet_path.as_posix()
   elif cadet_path.exists() and not lwe_path.exists():
       print('CADET was found but createLWE.exe was not found. Please make sure that none of the files have been moved.')
   else:
       print('CADET could not be found. Please check the bin path')
       
   print('\nSimulation Running!\n')

   # Create model object
   model = Cadet()

   # Number of unit operations
   model.root.input.model.nunits = 4

   # Inlet for SMA ############################################################
   model.root.input.model.unit_000.unit_type = 'INLET'
   model.root.input.model.unit_000.ncomp = len(ms['feed_conc_mM']) + 1
   model.root.input.model.unit_000.inlet_type = 'PIECEWISE_CUBIC_POLY'
   
   # Inlet sections
   # Section 1: Loading phase
   model.root.input.model.unit_000.sec_000.const_coeff = [ms['load_salt']] + ms['feed_conc_mM'] # [mol / m^3] = [mM]
   model.root.input.model.unit_000.sec_000.lin_coeff = [0] + [0 for c in ms['feed_conc']]
   model.root.input.model.unit_000.sec_000.quad_coeff = [0] + [0 for c in ms['feed_conc']]
   model.root.input.model.unit_000.sec_000.cube_coeff = [0] + [0 for c in ms['feed_conc']]
   
   # Section 2: Washing phase [no protein feed]
   model.root.input.model.unit_000.sec_001.const_coeff = [ms['wash_salt']] + [0 for c in ms['feed_conc']] # [mol / m^3] component 1, component 2

   # Section 3: step 1
   model.root.input.model.unit_000.sec_002.const_coeff = [ms['step1_salt']] + [0 for c in ms['feed_conc']] # [mol / m^3] component 1, component 2
   
   # Section 4: step 2
   model.root.input.model.unit_000.sec_003.const_coeff = [ms['step2_salt']] + [0 for c in ms['feed_conc']] # [mol / m^3] component 1, component 2

   # Section 5: step 3
   model.root.input.model.unit_000.sec_004.const_coeff = [ms['step3_salt']] + [0 for c in ms['feed_conc']] # [mol / m^3] component 1, component 2
   
   ## Construct CSTR (mixer) unit operation (tank model) ######################
   model.root.input.model.unit_001.unit_type = 'CSTR'
   model.root.input.model.unit_001.ncomp = len(ms['feed_conc_mM']) + 1

   model.root.input.model.unit_001.init_c = [ms['load_salt']] + [0 for c in ms['feed_conc']] # [mol / m^3], Initial concentrations for each component in liquid phase
   model.root.input.model.unit_001.init_q = [] # [mol / m^3], Initial concentrations for each component in solid phase
                      
   model.root.input.model.unit_001.init_volume = 1e-6 # CSTR volume in [m^3] 
   
   # General Rate Model #######################################################
   model.root.input.model.unit_002.unit_type = 'GENERAL_RATE_MODEL'
   
   model.root.input.model.unit_002.ncomp = len(ms['feed_conc_mM']) + 1

   ## Geometry
   model.root.input.model.unit_002.col_length = ms['col_length']                       # m
   model.root.input.model.unit_002.cross_section_area = ms['Ac_col']                   # m
   model.root.input.model.unit_002.col_porosity = ms['Ee']                             # -
   model.root.input.model.unit_002.par_porosity = ms['Ep']                             # -
   model.root.input.model.unit_002.par_radius = 0.5*ms['particle_diameter']            # m
                                                                   
   ## Transport
   model.root.input.model.unit_002.col_dispersion = ms['Dax_col'] # m^2 / s (interstitial volume)
   model.root.input.model.unit_002.film_diffusion = ms['k_film'] # m / s
   model.root.input.model.unit_002.par_diffusion = [ms['salt_pore_diff']] + list(ms['pore_diff']) # m^2 / s (mobile phase)  
       
   model.root.input.model.unit_002.adsorption.is_kinetic = True    # Kinetic binding

   ## Isotherm ################################################################
   model.root.input.model.unit_002.adsorption_model = 'STERIC_MASS_ACTION'
   model.root.input.model.unit_002.adsorption.SMA_KKIN = [0,] + ms['k_kin']
   model.root.input.model.unit_002.adsorption.SMA_KA = [0,] + list(ms['Keq'])
   model.root.input.model.unit_002.adsorption.SMA_KD = [0,] + [1] * len(ms['Keq'])
   model.root.input.model.unit_002.adsorption.SMA_NU = [0,] + list(ms['characteristic_charge'])
   model.root.input.model.unit_002.adsorption.SMA_SIGMA = [0,] + list(ms['steric_factor'])
   model.root.input.model.unit_002.adsorption.SMA_LAMBDA = ms['ionic_capacity'] # [mol / L solid phase]  

   ##############################################################################
    
   ## Initial conditions
   model.root.input.model.unit_002.init_c = [ms['load_salt']] + [0 for c in ms['feed_conc']]
   model.root.input.model.unit_002.init_q = []

   ## Discretization
   ### Grid cells
   model.root.input.model.unit_002.discretization.ncol = ms['axial_nodes']
   model.root.input.model.unit_002.discretization.npar = ms['particle_nodes']

   ### Bound states
   model.root.input.model.unit_002.discretization.nbound = [1] + [1 for c in ms['feed_conc']]

   ### Other options
   model.root.input.model.unit_002.discretization.par_disc_type = 'EQUIDISTANT_PAR'    
   model.root.input.model.unit_002.discretization.use_analytic_jacobian = 1
   model.root.input.model.unit_002.discretization.reconstruction = 'WENO'
   model.root.input.model.unit_002.discretization.gs_type = 1
   model.root.input.model.unit_002.discretization.max_krylov = 0
   model.root.input.model.unit_002.discretization.max_restarts = 10
   model.root.input.model.unit_002.discretization.schur_safety = 1e-8

   model.root.input.model.unit_002.discretization.weno.boundary_model = 0
   model.root.input.model.unit_002.discretization.weno.weno_eps = 1e-10
   model.root.input.model.unit_002.discretization.weno.weno_order = 3

   ## Outlet Model ############################################################
   model.root.input.model.unit_003.unit_type = 'OUTLET'
   model.root.input.model.unit_003.ncomp = len(ms['feed_conc_mM']) + 1
 
   ############################################################################
   
   # Sections 
   model.root.input.solver.sections.nsec = 5
   model.root.input.solver.sections.section_times = [0.0, ms['load_time'], ms['wash_time'], ms['step1_time'], ms['step2_time'], ms['step3_time']]   # s
   model.root.input.solver.sections.section_continuity = [0, 0, 0, 0, 0]
   
   # units: inlet > Column > Outlet

   # Switches
   model.root.input.model.connections.nswitches = 1
    
   # mixer is not connected during loading section
   model.root.input.model.connections.switch_000.section = 0
   model.root.input.model.connections.switch_000.connections = [
       0, 1, -1, -1, ms['flow_rate'],  # [unit_000, unit_001, all components, all components, Q/ m^3*s^-1 
       1, 2, -1, -1, ms['flow_rate'],  # [unit_001, unit_002, all components, all components, Q/ m^3*s^-1
       2, 3, -1, -1, ms['flow_rate'],  # [unit_002, unit_003, all components, all components, Q/ m^3*s^-1
       ]  

   # Solver settings
   model.root.input.model.solver.gs_type = 1
   model.root.input.model.solver.max_krylov = 0
   model.root.input.model.solver.max_restarts = 10
   model.root.input.model.solver.schur_safety = 1e-8

   # Number of cores for parallel simulation
   model.root.input.solver.nthreads = 1

   # Tolerances for the time integrator
   model.root.input.solver.time_integrator.abstol = 1e-6       # usually 1e-6
   model.root.input.solver.time_integrator.algtol = 1e-10           # usually 1e-10
   model.root.input.solver.time_integrator.reltol = 1e-6       # usually 1e-6
   model.root.input.solver.time_integrator.init_step_size = 1e-6    # usually 1e-6
   model.root.input.solver.time_integrator.max_steps = 1000000      # usually 1000000

   # Return data
   # dont return data except for outlet unit
   model.root.input['return'].split_components_data = 0
   model.root.input['return'].split_ports_data = 0
   model.root.input['return'].unit_000.write_solution_bulk = 0
   model.root.input['return'].unit_000.write_solution_inlet = 1
   model.root.input['return'].unit_000.write_solution_outlet = 0
   
   model.root.input['return'].unit_001.write_solution_bulk = 0
   model.root.input['return'].unit_001.write_solution_inlet = 0
   model.root.input['return'].unit_001.write_solution_outlet = 0
   
   # return outlet data for outlet unit
   model.root.input['return'].unit_002.write_solution_bulk = 0
   model.root.input['return'].unit_002.write_solution_inlet = 0
   model.root.input['return'].unit_002.write_solution_outlet = 0
   
   # return outlet data for outlet unit
   model.root.input['return'].unit_003.write_solution_bulk = 0
   model.root.input['return'].unit_003.write_solution_inlet = 0
   model.root.input['return'].unit_003.write_solution_outlet = 1
   
   # Solution times
   model.root.input.solver.user_solution_times = np.linspace(0, ms['step3_time'], 1000)

   # Save and run simulation
   model.filename = 'model.h5'
   model.save()

   data = model.run()

   if data.returncode == 0:
       model.load()   
   else:
       print(data)
       raise Exception('Simulation failed')
    
   sim_solution = model.root.output.solution
   return sim_solution

