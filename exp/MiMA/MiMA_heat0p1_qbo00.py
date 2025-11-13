import os

import numpy as np

from isca import IscaCodeBase, DiagTable, Experiment, Namelist, GFDL_BASE

NCORES = 32

# a CodeBase can be a directory on the computer,
# useful for iterative development
cb = IscaCodeBase.from_directory(GFDL_BASE)

# or it can point to a specific git repo and commit id.
# This method should ensure future, independent, reproducibility of results.
# cb = DryCodeBase.from_repo(repo='https://github.com/isca/isca', commit='isca1.1')

# compilation depends on computer specific settings.  The $GFDL_ENV
# environment variable is used to determine which `$GFDL_BASE/src/extra/env` file
# is used to load the correct compilers.  The env file is always loaded from
# $GFDL_BASE and not the checked out git repo.

# create an Experiment object to handle the configuration of model parameters
# and output diagnostics
exp = Experiment('mima_heat0p1_qbo00', codebase=cb)

exp.inputfiles = [os.path.join(GFDL_BASE,'input/rrtm_input_files/ozone_1990.nc')]

#Tell model how to write diagnostics
diag = DiagTable()
diag.add_file('atmos_monthly', 30, 'days', time_units='days')

# === Basic state variables ===
diag.add_field('dynamics', 'ps', time_avg=True)
diag.add_field('dynamics', 'ucomp', time_avg=True)
diag.add_field('dynamics', 'vcomp', time_avg=True)
diag.add_field('dynamics', 'temp', time_avg=True)
diag.add_field('dynamics', 'omega', time_avg=True)
diag.add_field('dynamics', 'height', time_avg=True)

# === Eddy momentum fluxes (CRITICAL for EP flux!) ===
diag.add_field('dynamics', 'ucomp_vcomp', time_avg=True)  # u'v' - horizontal EP flux
diag.add_field('dynamics', 'vcomp_temp', time_avg=True)   # v'T' - vertical EP flux component
diag.add_field('dynamics', 'vcomp_omega', time_avg=True)  # v'ω' - alternative vertical flux

# === Eddy heat fluxes ===
diag.add_field('dynamics', 'ucomp_temp', time_avg=True)   # u'T'
diag.add_field('dynamics', 'omega_temp', time_avg=True)   # ω'T'

# === Variances (for eddy kinetic energy, etc.) ===
diag.add_field('dynamics', 'ucomp_sq', time_avg=True)
diag.add_field('dynamics', 'vcomp_sq', time_avg=True)
diag.add_field('dynamics', 'temp_sq', time_avg=True)
diag.add_field('dynamics', 'omega_sq', time_avg=True)

# === Circulation diagnostics ===
diag.add_field('dynamics', 'vor', time_avg=True)
diag.add_field('dynamics', 'div', time_avg=True)

# === Pressure coordinates ===
diag.add_field('dynamics', 'pres_full', time_avg=True)
diag.add_field('dynamics', 'pres_half', time_avg=True)

# === HS forcing diagnostics ===
diag.add_field('hs_forcing', 'teq', time_avg=True)
diag.add_field('hs_forcing', 'tdt_ndamp', time_avg=True)

# === Spectral truncation info (static fields) ===
diag.add_field('dynamics', 'bk')
diag.add_field('dynamics', 'pk')

exp.diag_table = diag


#Empty the run directory ready to run
exp.clear_rundir()

#Define values for the 'core' namelist
exp.namelist = namelist = Namelist({
    'main_nml': {
        'days'   : 30,
        'hours'  : 0,
        'minutes': 0,
        'seconds': 0,
        'dt_atmos':600,
        'current_date' : [1,1,1,0,0,0],
        'calendar' : 'thirty_day'
    },

    'idealized_moist_phys_nml': {
        'two_stream_gray': False,
        'do_rrtm_radiation': True,    #Use RRTM radiation, not grey
        'convection_scheme': 'SIMPLE_BETTS_MILLER',     #Use the simple Betts Miller convection scheme
        'do_damping': True,
        'turb':True,
        'mixed_layer_bc':True,
        'do_virtual' :False,
        'do_simple': True,
        'roughness_mom':3.21e-05,
        'roughness_heat':3.21e-05,
        'roughness_moist':3.21e-05,                
    },

    'vert_turb_driver_nml': {
        'do_mellor_yamada': False,     # default: True
        'do_diffusivity': True,        # default: False
        'do_simple': True,             # default: False
        'constant_gust': 0.0,          # default: 1.0
        'use_tau': False
    },
    
    'diffusivity_nml': {
        'do_entrain':False,
        'do_simple': True,
    },

    'surface_flux_nml': {
        'use_virtual_temp': False,
        'do_simple': True,
        'old_dtaudv': True    
    },

    'atmosphere_nml': {
        'idealized_moist_model': True
    },

    #Use a large mixed-layer depth, and the Albedo of the CTRL case in Jucker & Gerber, 2017
    'mixed_layer_nml': {
        'depth': 100,
        'albedo_value': 0.205,
        'tconst' : 285.,
        'prescribe_initial_dist':True,
        'evaporation':True,
        'do_qflux': True        
    },

    'qe_moist_convection_nml': {
        'rhbm':0.7,
        'Tmin':160.,
        'Tmax':350.   
    },
    
    'lscale_cond_nml': {
        'do_simple':True,
        'do_evap':True
    },
    
    'sat_vapor_pres_nml': {
        'do_simple':True
    },
    
    'damping_driver_nml': {
        'do_rayleigh': True,
        'trayfric': -0.25,              # quarter-day damping to prevent cold stratosphere crash; neg. value: time in *days*
        'sponge_pbottom':  800., # pws sponge bottom to 8 hPa for consistency
        'do_conserve_energy': True,
        'do_ewa_htg': True,
        'do_sin_qbo': True,
        'h_amp':0.1/86400.,
        'qbo_amp': 0.,
    },

    'qflux_nml': {
        'qflux_amp': 30.0
    },

    'rrtm_radiation_nml': {
        'solr_cnst': 1360,  #s set solar constant to 1360, rather than default of 1368.22
        'dt_rad': 7200, #Use long RRTM timestep
        'do_read_ozone':True,
        'ozone_file':'ozone_1990'
    },

    # FMS Framework configuration
    'diag_manager_nml': {
        'mix_snapshot_average_fields': False  # time avg fields are labelled with time in middle of window
    },

    'fms_nml': {
        'domains_stack_size': 600000                        # default: 0
    },

    'fms_io_nml': {
        'threading_write': 'single',                         # default: multi
        'fileset_write': 'single',                           # default: multi
    },

    'spectral_dynamics_nml': {
        'damping_order': 4,             
        'water_correction_limit': 200.e2,
        'reference_sea_level_press':1.0e5,
        'num_levels':50,               # pws levels from Frierson 2006
        'valid_range_t':[100.,800.],
        'initial_sphum':[2.e-6],
        'vert_coord_option':'input', #Use the vertical levels from Frierson 2006
        'robert_coeff':0.03
    },
    #'vert_coordinate_nml': {
    #    'bk': [0.000000, 0.0117665, 0.0196679, 0.0315244, 0.0485411, 0.0719344, 0.1027829, 0.1418581, 0.1894648, 0.2453219, 0.3085103, 0.3775033, 0.4502789, 0.5244989, 0.5977253, 0.6676441, 0.7322627, 0.7900587, 0.8400683, 0.8819111, 0.9157609, 0.9422770, 0.9625127, 0.9778177, 0.9897489, 1.0000000],
    #    'pk': [0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000],
#       }
    'vert_coordinate_nml': {
        'bk': [0.000000, 0.006777, 0.006896, 0.007098, 0.007389, 0.007778, 0.008278, 0.008904, 0.009676, 0.010620, 0.011767, 0.013154, 0.014829, 0.016849, 0.019284, 0.022216, 0.025745, 0.029990, 0.035091, 0.041212, 0.048541, 0.057295, 0.067716, 0.080072, 0.094651, 0.111759, 0.131704, 0.154786, 0.181278, 0.211406, 0.245322, 0.283080, 0.324608, 0.369689, 0.417937, 0.468794, 0.521533, 0.575275, 0.629024, 0.681712, 0.732263, 0.779655, 0.822998, 0.861593, 0.894993, 0.923042, 0.945905, 0.964070, 0.978349, 0.989857, 1.000000],
        'pk': [0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000],
    }
    
    
})

#Lets do a run!
if __name__=="__main__":
    exp.run(1, use_restart=True, num_cores=NCORES)
    for i in range(2,757): # 84 months spinup + 2*28*12
        exp.run(i, num_cores=NCORES)
