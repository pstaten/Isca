#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Plot the imposed QBO-like wind and SAI-like heating STRUCTURE for one run config.

Reads a run's namelist (&damping_driver_nml) and plots, in several panels, the structure of
that config's imposed forcings: the fake_qbo descending-wind pattern and the ewa_heating
(Bednarz-profile) SAI-like heating. Use this to verify what a single experiment actually forces.
For the sweep-wide MAP of heating locations, see plot_heating_location_grid.py instead.
"""

import sys
import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

def parse_namelist(namelist_path):
    """
    Parse the namelist file and extract damping_driver_nml parameters.

    Parameters
    ----------
    namelist_path : str
        Path to the namelist Python file

    Returns
    -------
    dict
        Dictionary containing the damping_driver_nml parameters
    """
    # Read the file and extract the damping_driver_nml section
    with open(namelist_path, 'r') as f:
        content = f.read()

    # Find the damping_driver_nml section
    # Look for 'damping_driver_nml': {
    import re

    # Pattern to match the damping_driver_nml dictionary
    pattern = r"'damping_driver_nml'\s*:\s*\{([^}]+)\}"
    match = re.search(pattern, content)

    if not match:
        raise ValueError("Could not find 'damping_driver_nml' in namelist file")

    # Extract the dictionary content
    dict_content = match.group(1)

    # Parse key-value pairs
    params = {}

    # Pattern to match key: value pairs
    # Handles various types: booleans, integers, floats
    param_pattern = r"'(\w+)'\s*:\s*([^,\n]+)"

    for match in re.finditer(param_pattern, dict_content):
        key = match.group(1)
        value_str = match.group(2).strip()

        # Parse the value
        if value_str == 'True':
            value = True
        elif value_str == 'False':
            value = False
        elif '/' in value_str:
            # Handle expressions like 1.0/86400.
            value = eval(value_str)
        else:
            # Try to parse as number
            try:
                if '.' in value_str or 'e' in value_str.lower():
                    value = float(value_str)
                else:
                    value = int(value_str)
            except ValueError:
                value = value_str

        params[key] = value

    return params

def compute_ewa_heating(lat, p_full, h_amp, p_s, p_center, p_width, lat_width):
    """
    Compute the EWA heating forcing.

    This implements the heating profile from Ewa Bednarz' work:
    Q(lat, p) = Q_max * exp[-(lat^2/(2*0.4^2) + (p/p_s - p_center/p_s)^2/(2*p_width^2))]

    Parameters
    ----------
    lat : ndarray (nlat,)
        Latitude in radians
    p_full : ndarray (nlev,)
        Pressure levels in Pa
    h_amp : float
        Heating amplitude (K/s)
    p_s : float
        Surface pressure (Pa)
    p_center : float
        Center pressure level (Pa)
    p_width : float
        Pressure width parameter
    lat_width : float
        Latitude width parameter (radians)

    Returns
    -------
    ndarray (nlat, nlev)
        Heating rate (K/day) as a function of latitude and pressure
    """
    nlat = len(lat)
    nlev = len(p_full)

    ttnd = np.zeros((nlat, nlev))

    for k in range(nlev):
        for j in range(nlat):
            # Vertical component
            p_factor = np.exp(-(p_full[k]/p_s - p_center/p_s)**2 / (2*p_width**2))

            # Latitudinal component
            lat_factor = np.exp(-lat[j]**2 / (2*lat_width**2))

            # Combine
            ttnd[j, k] = h_amp * p_factor * lat_factor

    # Convert from K/s to K/day
    ttnd_per_day = ttnd * 86400.0

    return ttnd_per_day

def compute_fake_qbo(lat, p_full, z_full, qbo_amp, time_days=0):
    """
    Compute the fake QBO forcing.

    Parameters
    ----------
    lat : ndarray (nlat,)
        Latitude in radians
    p_full : ndarray (nlev,)
        Pressure levels in Pa
    z_full : ndarray (nlev,)
        Height levels in m
    qbo_amp : float
        QBO amplitude (m/s)
    time_days : float
        Time in days (default: 0)

    Returns
    -------
    qbo : ndarray (nlat, nlev)
        QBO wind profile (m/s) at time t
    qbo_time : ndarray (ntimes, nlev_qbo)
        QBO wind profile at equator vs time
    times : ndarray (ntimes,)
        Time array in days
    p_qbo : ndarray (nlev_qbo,)
        Pressure levels where QBO is applied
    """
    nlat = len(lat)
    nlev = len(p_full)

    # QBO parameters (from fake_qbo subroutine)
    daysperyear = 360  # thirty_day calendar
    daypsec = 1.0 / 86400.0

    # QBO frequency (28-month period)
    nu = 2.0 * np.pi / ((28.0/12.0) * daysperyear / daypsec)  # rad/s

    # QBO vertical wavenumber (40 km vertical wavelength)
    m = -2.0 * np.pi / 40000.0  # rad/m

    # Altitude profile parameters
    p_b = 8000.0  # Pa (80 hPa)
    del_p_b = 500.0  # Pa (5 hPa)
    p_t = 800.0  # Pa (8 hPa)
    del_p_t = 50.0  # Pa (0.5 hPa)

    # Latitude profile parameters
    lat_n = np.pi * 12.0 / 180.0  # 12 deg N in radians
    lat_s = np.pi * (-12.0) / 180.0  # 12 deg S in radians
    del_lat = np.pi * 2.0 / 180.0  # 2 deg in radians

    # Compute QBO profile at time t
    seconds = 0
    days = time_days
    qbo_phase = np.zeros((nlat, nlev))

    for k in range(nlev):
        for j in range(nlat):
            qbo_phase[j, k] = qbo_amp * np.cos(m * z_full[k] - nu * (days*24*3600 + seconds))

    # Altitude profile
    alt_profile = np.zeros(nlev)
    for k in range(nlev):
        alt_profile[k] = (1 - np.tanh((p_full[k] - p_b) / del_p_b)) * \
                        (1 + np.tanh((p_full[k] - p_t) / del_p_t))

    # Latitude profile
    lat_profile = np.zeros(nlat)
    for j in range(nlat):
        lat_profile[j] = (1 + np.tanh((lat[j] - lat_s) / del_lat)) * \
                        (1 - np.tanh((lat[j] - lat_n) / del_lat))

    # Apply profiles to QBO
    qbo = np.zeros((nlat, nlev))
    for k in range(nlev):
        for j in range(nlat):
            qbo[j, k] = qbo_phase[j, k] * alt_profile[k] * lat_profile[j] / 16.0

    # Compute QBO at equator vs time
    # Find equator index
    eq_idx = np.argmin(np.abs(lat))

    # Time array (0 to 84 months = 7 years with 28-month period)
    times = np.linspace(0, 84*30, 500)  # 30-day months

    # Compute QBO time series for all levels at equator
    qbo_time = np.zeros((len(times), nlev))
    for i, t in enumerate(times):
        days_t = t
        for k in range(nlev):
            qbo_time[i, k] = qbo_amp * np.cos(m * z_full[k] - nu * days_t*24*3600) * \
                             alt_profile[k] * lat_profile[eq_idx] / 16.0

    return qbo, qbo_time, times, p_full

def compute_canonical_qbo_time(p_full, times_days, table_path, qbo_tab_scale=1.0):
    """
    Equatorial time series of the TABULATED canonical QBO nudging target (do_tab_qbo),
    mirroring the Fortran: periodic-linear in time (30-day months), linear in ln(p)
    clamped outside the table range, then the same tanh altitude window and equatorial
    latitude factor as the sinusoidal fake_qbo.

    Returns ndarray (ntimes, nlev) in m/s.
    """
    with open(table_path) as f:
        lines = [l for l in f if not l.startswith('#')]
    nlev_tab, ntime = (int(x) for x in lines[0].split())
    plev = np.array([float(x) for x in lines[1].split()])  # hPa
    utab = np.array([[float(x) for x in lines[2 + j].split()] for j in range(nlev_tab)])
    order = np.argsort(plev)  # ascending p for interpolation in ln(p)
    lnp_tab = np.log(plev[order] * 100.0)
    utab = utab[order]

    # same windows as compute_fake_qbo
    p_b, del_p_b, p_t, del_p_t = 8000.0, 500.0, 800.0, 50.0
    alt_profile = (1 - np.tanh((p_full - p_b) / del_p_b)) * (1 + np.tanh((p_full - p_t) / del_p_t))
    lat_eq_factor = (1 + np.tanh((0 - np.pi*(-12.0)/180.0) / (np.pi*2.0/180.0))) * \
                    (1 - np.tanh((0 - np.pi*12.0/180.0) / (np.pi*2.0/180.0)))

    lnp = np.log(np.maximum(p_full, 1.0))
    out = np.zeros((len(times_days), len(p_full)))
    for i, t in enumerate(times_days):
        tmon = np.mod(t / 30.0, ntime)
        it1 = int(tmon) % ntime
        it2 = (it1 + 1) % ntime
        wt = tmon - int(tmon)
        ut = (1 - wt) * utab[:, it1] + wt * utab[:, it2]
        uz = np.interp(lnp, lnp_tab, ut)  # np.interp clamps outside the table range
        out[i] = qbo_tab_scale * uz * alt_profile * lat_eq_factor / 16.0
    return out


def create_pressure_height_grid(num_levels=50, scale_heights=6.0, surf_res=0.1, exponent=2.5):
    """
    Create pressure and height grids matching the model configuration.

    Parameters
    ----------
    num_levels : int
        Number of vertical levels
    scale_heights : float
        Model top height in scale heights
    surf_res : float
        Fraction of range concentrated near surface
    exponent : float
        Surface concentration exponent

    Returns
    -------
    p_full : ndarray
        Pressure at full levels (Pa)
    z_full : ndarray
        Height at full levels (m)
    """
    # Create sigma levels (uneven_sigma coordinate)
    k = np.arange(1, num_levels + 1)

    # Uneven sigma coordinate (approximation)
    # This creates enhanced resolution near the surface and at the top
    sigma = ((k - 0.5) / num_levels) ** exponent

    # Reference surface pressure
    p_surf = 1.0e5  # Pa

    # Pressure at full levels
    # For simplicity, use a hybrid coordinate approximation
    # p approx p_surf * sigma for lower levels
    # Scale to get stratospheric levels
    p_full = p_surf * np.exp(-scale_heights * (1 - sigma[::-1]))

    # Height using hypsometric equation (H = scale_height * ln(p_surf/p))
    H_scale = 7000.0  # meters (scale height)
    z_full = H_scale * np.log(p_surf / p_full)

    return p_full, z_full

def main():
    if len(sys.argv) == 1:
        namelist_path = "exp/frierson50lev/frierson50lev_heat0p1_qbo20.py"
        print(f"No namelist file provided. Using default: {namelist_path}")
    elif len(sys.argv) == 2:
        namelist_path = sys.argv[1]
    else:
        print("Usage: python plot_imposed_forcings.py <namelist_file.py>")
        print("Example: python plot_imposed_forcings.py exp/frierson50lev/frierson50lev_heat1p0_qbo20.py")
        sys.exit(1)

    if not os.path.exists(namelist_path):
        print(f"Error: Namelist file '{namelist_path}' not found")
        sys.exit(1)

    # Parse namelist
    params = parse_namelist(namelist_path)

    # Extract parameters with defaults from damping_driver.f90
    h_amp = params.get('h_amp', 0.1/86400.0)  # K/s
    p_s = 100000.0  # Pa (hardcoded in Fortran, not in namelist)
    p_center = 5000.0  # Pa (hardcoded in Fortran, not in namelist)
    p_width = 0.0175  # (hardcoded in Fortran, not in namelist)
    lat_width = 0.4  # radians (hardcoded in Fortran, not in namelist)
    qbo_amp = params.get('qbo_amp', 20.0)  # m/s

    do_ewa_htg = params.get('do_ewa_htg', False)
    do_sin_qbo = params.get('do_sin_qbo', False)

    # Create grids
    num_levels = 50  # From namelist: spectral_dynamics_nml['num_levels']
    nlat = 64  # Typical for T42 resolution

    # Latitude grid
    lat = np.linspace(-np.pi/2, np.pi/2, nlat)
    lat_deg = np.degrees(lat)

    # Pressure and height grids
    p_full, z_full = create_pressure_height_grid(num_levels=num_levels)
    p_hPa = p_full / 100.0  # Convert to hPa

    # Compute forcings
    if do_ewa_htg:
        ewa_heating = compute_ewa_heating(lat, p_full, h_amp, p_s, p_center, p_width, lat_width)
    else:
        ewa_heating = np.zeros((nlat, num_levels))

    if do_sin_qbo:
        qbo, qbo_time, times, p_qbo = compute_fake_qbo(lat, p_full, z_full, qbo_amp, time_days=0)
        p_qbo_hPa = p_qbo / 100.0
    else:
        qbo = np.zeros((nlat, num_levels))
        qbo_time = None
        times = None
        p_qbo_hPa = None

    # Pad pressure grid to extend from 1000 hPa to 1 hPa with zeros
    p_extended_hPa = np.concatenate([
        np.array([1000.0]),  # Add surface level
        p_hPa,
        np.array([1.0])      # Add top level
    ])

    # Pad the forcing arrays with zeros at top and bottom
    ewa_heating_padded = np.zeros((nlat, num_levels + 2))
    ewa_heating_padded[:, 1:-1] = ewa_heating

    qbo_padded = np.zeros((nlat, num_levels + 2))
    qbo_padded[:, 1:-1] = qbo

    # Pad QBO time series with zeros at top and bottom
    if qbo_time is not None:
        qbo_time_padded = np.zeros((len(times), num_levels + 2))
        qbo_time_padded[:, 1:-1] = qbo_time
    else:
        qbo_time_padded = None

    # Create figure with three panels
    # Increase font sizes for readability
    plt.rcParams.update({
        'font.size': 18,
        'axes.titlesize': 24,
        'axes.labelsize': 18,
        'xtick.labelsize': 16,
        'ytick.labelsize': 16,
        'legend.fontsize': 16,
    })
    # Optional 4th panel: the canonical (tabulated, observed-composite) QBO cycle
    canonical_table = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                   'input', 'qbo_canonical_cycle.txt')
    have_canonical = os.path.exists(canonical_table) and do_sin_qbo
    npanels = 4 if have_canonical else 3

    fig = plt.figure(figsize=(4.43 * npanels, 5.5))
    gs = GridSpec(1, npanels, figure=fig, wspace=0.3)

    ax1 = fig.add_subplot(gs[0, 0])
    ax2 = fig.add_subplot(gs[0, 1])
    ax3 = fig.add_subplot(gs[0, 2])
    ax4 = fig.add_subplot(gs[0, 3]) if have_canonical else None

    cb_pad = 0.09
    # Panel 1: EWA heating vs pressure and latitude
    # n_levels_ewa = 11
    # levels_ewa = np.linspace(0, h_amp*86400.*n_levels_ewa/(n_levels_ewa-1), n_levels_ewa)
    levels_ewa = np.arange(0.1, 1.1, 0.1)*(h_amp*86400.0)

    if do_ewa_htg:
        c1 = ax1.contourf(lat_deg, p_extended_hPa, ewa_heating_padded.T, levels=levels_ewa, cmap='magma_r')
        ax1.contour(lat_deg, p_extended_hPa, ewa_heating_padded.T, levels=levels_ewa, colors='k', linewidths=0.5, alpha=0.3)
        cb = plt.colorbar(c1, ax=ax1, orientation='horizontal', label=r'(K day$^{-1}$)', pad=cb_pad)
        ticks = cb.get_ticks()
        if len(ticks) >= 2:
            new_ticks = [ticks[0], ticks[-1]]
            cb.set_ticks(new_ticks)
            cb.set_ticklabels([f"{t:.3g}" for t in new_ticks])
    ax1.set_yscale('log')
    ax1.invert_yaxis()
    ax1.set_xticks([-90, -60, -30, 0, 30, 60, 90])
    ax1.set_xticklabels(['90°S', '', '', '0°', '', '', '90°N'])
    ax1.set_ylabel('Pressure (hPa)')
    ax1.set_title('(a) $T$ tendency')
    ax1.set_ylim(200, 2)
    ax1.set_yticks([100, 10])
    ax1.set_yticklabels(['100', '10'])
    ax1.grid(True, alpha=0.3)

    # Panel 2: QBO forcing at t=0 vs pressure and latitude
    levels_qbo = np.linspace(-qbo_amp, qbo_amp, 21)
    levels_qbo = np.arange(-18, 22, 4)
    levels_qbo = np.array([-20, -16, -12, -8, -4, 4, 8, 12, 16, 20])
    if do_sin_qbo:
        c2 = ax2.contourf(lat_deg, p_extended_hPa, qbo_padded.T, levels=levels_qbo, cmap='RdBu_r', extend='both')
        ax2.contour(lat_deg, p_extended_hPa, qbo_padded.T, levels=levels_qbo, colors='k', linewidths=0.5, alpha=0.3)
        cb = plt.colorbar(c2, ax=ax2, orientation='horizontal', label=r'(m s$^{-1}$)', pad=cb_pad)
        ticks = cb.get_ticks()
        if len(ticks) >= 2:
            new_ticks = [ticks[0], ticks[-1]]
            cb.set_ticks(new_ticks)
            cb.set_ticklabels([f"{t:.3g}" for t in new_ticks])
    ax2.set_yscale('log')
    ax2.invert_yaxis()
    ax2.set_xticks([-90, -60, -30, 0, 30, 60, 90])
    ax2.set_xticklabels(['90°S', '', '', '0°', '', '', '90°N'])
    # ax2.set_ylabel('Pressure (hPa)')
    ax2.set_title(r'(b) $\overline{u}_\mathrm{QBO}$, $t$=0')
    ax2.set_ylim(200, 2)
    # show plain decade labels instead of 10^2, 10^1
    ax2.set_yticks([100, 10])
    ax2.set_yticklabels(['100', '10'])
    ax2.grid(True, alpha=0.3)

    # Panel 3: QBO forcing at equator vs pressure and time
    if do_sin_qbo and qbo_time_padded is not None:
        times_months = times / 30.0  # Convert to months
        c3 = ax3.contourf(times_months, p_extended_hPa, qbo_time_padded.T, levels=levels_qbo, cmap='RdBu_r', extend='both')
        ax3.contour(times_months, p_extended_hPa, qbo_time_padded.T, levels=levels_qbo, colors='k', linewidths=0.5, alpha=0.3)
        cb = plt.colorbar(c3, ax=ax3, orientation='horizontal', label=r'(m s$^{-1}$)', pad=cb_pad)
        ticks = cb.get_ticks()
        if len(ticks) >= 2:
            new_ticks = [ticks[0], ticks[-1]]
            cb.set_ticks(new_ticks)
            cb.set_ticklabels([f"{t:.3g}" for t in new_ticks])
    ax3.set_yscale('log')
    ax3.invert_yaxis()
    # ax3.set_ylabel('Pressure (hPa)')
    ax3.set_xticks([0, 12, 24, 36, 48, 60, 72, 84])
    ax3.set_xticklabels(['yr 0', '1', '2', '3', '4', '5', '6', '7'])
    ax3.set_title(r'(c) $\overline{u}_\mathrm{QBO}$, $\phi$=0')
    ax3.set_ylim(200, 2)
    ax3.set_yticks([100, 10])
    ax3.set_yticklabels(['100', '10'])
    ax3.grid(True, alpha=0.3)

    # Panel 4: canonical (observed-composite) QBO at equator vs pressure and time
    if have_canonical:
        qbo_can = compute_canonical_qbo_time(p_full, times, canonical_table)
        qbo_can_padded = np.zeros((len(times), num_levels + 2))
        qbo_can_padded[:, 1:-1] = qbo_can
        times_months = times / 30.0
        c4 = ax4.contourf(times_months, p_extended_hPa, qbo_can_padded.T, levels=levels_qbo,
                          cmap='RdBu_r', extend='both')
        ax4.contour(times_months, p_extended_hPa, qbo_can_padded.T, levels=levels_qbo,
                    colors='k', linewidths=0.5, alpha=0.3)
        cb = plt.colorbar(c4, ax=ax4, orientation='horizontal', label=r'(m s$^{-1}$)', pad=cb_pad)
        ticks = cb.get_ticks()
        if len(ticks) >= 2:
            cb.set_ticks([ticks[0], ticks[-1]])
            cb.set_ticklabels([f"{t:.3g}" for t in [ticks[0], ticks[-1]]])
        ax4.set_yscale('log')
        ax4.invert_yaxis()
        ax4.set_xticks([0, 12, 24, 36, 48, 60, 72, 84])
        ax4.set_xticklabels(['yr 0', '1', '2', '3', '4', '5', '6', '7'])
        ax4.set_title(r'(d) canonical $\overline{u}_\mathrm{QBO}$ (obs)')
        ax4.set_ylim(200, 2)
        ax4.set_yticks([100, 10])
        ax4.set_yticklabels(['100', '10'])
        ax4.grid(True, alpha=0.3)

    plt.tight_layout()

    # Save figure
    output_name = os.path.splitext(os.path.basename(namelist_path))[0] + '_forcings'
    plt.savefig(output_name + '.png', dpi=300, bbox_inches='tight')
    plt.savefig(output_name + '.pdf', bbox_inches='tight')
    print(f"Saved plot to: {output_name}")
    plt.close()

if __name__ == "__main__":
    main()
