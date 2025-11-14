import numpy as np

from isca import DryCodeBase, DiagTable, Experiment, Namelist, GFDL_BASE

NCORES = 16
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
        'h_amp': 0.1/86400.,
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
        'damping_order': 4,             
        'water_correction_limit': 200.e2,
        'reference_sea_level_press':1.0e5,
        'num_levels':50,               # pws levels from Frierson 2006
        'valid_range_t':[100.,800.],
        'initial_sphum':[2.e-6],
        'vert_coord_option':'uneven_sigma', # trying to get a better stratosphere
        'robert_coeff':0.03,
        'scale_heights': 6.0,      # Model top at ~e^(-6) ≈ 0.25 hPa (much higher!)
        'surf_res': 0.1,           # 10% of range concentrated near surface
        'exponent': 2.5,           # Moderate surface concentration
    }
})
exp.set_resolution(*RESOLUTION)

#Lets do a run!
if __name__ == '__main__':
    exp.run(1, num_cores=NCORES, use_restart=True)
    for i in range(2, 757):
        exp.run(i, num_cores=NCORES)  # use the restart i-1 by default
