import numpy as np

from isca import DryCodeBase, DiagTable, Experiment, Namelist, GFDL_BASE

NCORES = 4
RESOLUTION = 'T42', 50  # T42 horizontal resolution, 25 levels in pressure

# a CodeBase can be a directory on the computer,
# useful for iterative development
cb = DryCodeBase.from_directory(GFDL_BASE)

# or it can point to a specific git repo and commit id.
# This method should ensure future, independent, reproducibility of results.
# cb = DryCodeBase.from_repo(repo='https://github.com/isca/isca', commit='isca1.1')

# compilation depends on computer specific settings.  The $GFDL_ENV
# environment variable is used to determine which `$GFDL_BASE/src/extra/env` file
# is used to load the correct compilers.  The env file is always loaded from
# $GFDL_BASE and not the checked out git repo.

cb.compile()  # compile the source code to working directory $GFDL_WORK/codebase

# create an Experiment object to handle the configuration of model parameters
# and output diagnostics

exp = Experiment('hs_heat0p1_qbo20', codebase=cb)

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

# define namelist values as python dictionary
# wrapped as a namelist object.
exp.namelist = namelist = Namelist({
    'main_nml':{
     'days'   : 30,
     'hours'  : 0,
     'minutes': 0,
     'seconds': 0,
     'dt_atmos':600, # pws reduced from 720 for stability with 50 levels
     'current_date' : [1,1,1,0,0,0],
     'calendar' : 'thirty_day'
    },

    'atmosphere_nml': {
        'idealized_moist_model': False  # False for Newtonian Cooling.  True for Isca/Frierson
    },

    # configure the relaxation profile
    'hs_forcing_nml': {
        't_zero': 315.,    # temperature at reference pressure at equator (default 315K)
        't_strat': 200.,   # stratosphere temperature (default 200K)
        'delh': 60.,       # equator-pole temp gradient (default 60K)
        'delv': 10.,       # lapse rate (default 10K)
        'eps': 0.,         # stratospheric latitudinal variation (default 0K)
        'sigma_b': 0.7,    # boundary layer friction height (default p/ps = sigma = 0.7)

        # negative sign is a flag indicating that the units are days
        'ka':   -15.,      # pws faster stratospheric relaxation (default 40 days)
        'ks':    -4.,      # Boundary layer dependent cooling timescale (default 4 days)
        'kf':   -1.,       # BL momentum frictional timescale (default 1 days)

        'do_conserve_energy':   True,  # convert dissipated momentum into heat (default True)
        'p_trop': 1.e4, # default
        'stratosphere_t_option': 'hs_like', # pws; colder CPT, classic framework, a bit less numerically stable
        'do_ewa_htg': True # run Ewa's heating
        'do_sin_qbo': True,  # enable QBO nudging
        'h_amp': 0.1,
        'qbo_amp': 20.0,  # QBO amplitude (m/s)
    },

    'diag_manager_nml': {
        'mix_snapshot_average_fields': False
    },

    'fms_nml': {
        'domains_stack_size': 600000                        # default: 0
    },

    'fms_io_nml': {
        'threading_write': 'single',                         # default: multi
        'fileset_write': 'single',                           # default: multi
    },

    'spectral_dynamics_nml': {
        'damping_order'           : 4,                      # default: 2
        'water_correction_limit'  : 200.e2,                 # default: 0
        'reference_sea_level_press': 1.0e5,                  # default: 101325
        'valid_range_t'           : [100., 800.],           # default: (100, 500)
        'initial_sphum'           : 0.0,                  # default: 0
        'vert_coord_option'       : 'input',         # default: 'even_sigma', pws changed from uneven_sigma
    },
    #'vert_coordinate_nml': {
    #    'bk': [0.000000, 0.0117665, 0.0196679, 0.0315244, 0.0485411, 0.0719344, 0.1027829, 0.1418581, 0.1894648, 0.2453219, 0.3085103, 0.3775033, 0.4502789, 0.5244989, 0.5977253, 0.6676441, 0.7322627, 0.7900587, 0.8400683, 0.8819111, 0.9157609, 0.9422770, 0.9625127, 0.9778177, 0.9897489, 1.0000000],
    #    'pk': [0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000],
#       }
    'vert_coordinate_nml': {
        'bk': [0.000000, 0.006777, 0.006896, 0.007098, 0.007389, 0.007778, 0.008278, 0.008904, 0.009676, 0.010620, 0.011767, 0.013154, 0.014829, 0.016849, 0.019284, 0.022216, 0.025745, 0.029990, 0.035091, 0.041212, 0.048541, 0.057295, 0.067716, 0.080072, 0.094651, 0.111759, 0.131704, 0.154786, 0.181278, 0.211406, 0.245322, 0.283080, 0.324608, 0.369689, 0.417937, 0.468794, 0.521533, 0.575275, 0.629024, 0.681712, 0.732263, 0.779655, 0.822998, 0.861593, 0.894993, 0.923042, 0.945905, 0.964070, 0.978349, 0.989857, 1.000000],
        'pk': [0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000],
    },

})

exp.set_resolution(*RESOLUTION)

#Lets do a run!
if __name__ == '__main__':
    exp.run(1, num_cores=NCORES, use_restart=True)
    for i in range(2, 757):
        exp.run(i, num_cores=NCORES)  # use the restart i-1 by default
