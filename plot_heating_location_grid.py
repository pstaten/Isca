#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Plot the GRID of SAI-like heating locations across the experiment sweep.

Draws one contour (0.75x max) of the ewa_heating (Bednarz-profile) forcing for each
(lat_center, p_center) combination in the sweep, all overlaid on a single latitude-pressure
plane -- i.e. a map of WHERE heating is applied across the experiments. For the QBO + heating
STRUCTURE of a single run's namelist, see plot_imposed_qbo_and_heating.py instead.
"""

import numpy as np
import matplotlib.pyplot as plt

def compute_ewa_heating(lat, p_full, h_amp, p_s, p_center, p_width, lat_center, lat_width):
    """
    Compute the EWA heating forcing.

    This implements the heating profile from Ewa Bednarz' work:
    Q(lat, p) = Q_max * exp[-(lat-lat_center)^2/(2*lat_width^2) + (p/p_s - p_center/p_s)^2/(2*p_width^2))]

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
    lat_center : float
        Center latitude (radians)
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

            # Latitudinal component (centered at lat_center)
            lat_factor = np.exp(-(lat[j] - lat_center)**2 / (2*lat_width**2))

            # Combine
            ttnd[j, k] = h_amp * p_factor * lat_factor

    # Convert from K/s to K/day
    ttnd_per_day = ttnd * 86400.0

    return ttnd_per_day

def analytical_integral(h_amp, lat_width, p_width, p_s):
    """
    Compute the analytical integral of the heating forcing.

    The definite integral (from -infinity to infinity in both lat and p dimensions) is:
    2π Q_max |σ_φ σ_p p_s|

    Parameters
    ----------
    h_amp : float
        Heating amplitude (K/s)
    lat_width : float
        Latitude width parameter (radians)
    p_width : float
        Pressure width parameter (dimensionless)
    p_s : float
        Surface pressure (Pa)

    Returns
    -------
    float
        Analytical integral (K Pa rad / s)
    """
    return 2.0 * np.pi * h_amp * np.abs(lat_width * p_width * p_s)

def numerical_integral(heating_per_sec, lat, p_full):
    """
    Compute the numerical integral of the heating on the discretized grid.

    Parameters
    ----------
    heating_per_sec : ndarray (nlat, nlev)
        Heating rate in K/s
    lat : ndarray (nlat,)
        Latitude in radians
    p_full : ndarray (nlev,)
        Pressure levels in Pa

    Returns
    -------
    float
        Numerical integral (K Pa rad / s)
    """
    nlat = len(lat)
    nlev = len(p_full)

    # Compute grid cell centers and widths
    # For latitude: use cell-centered approach
    if nlat > 1:
        dlat = np.diff(lat)
        # Extend to get width for all cells
        dlat_all = np.zeros(nlat)
        dlat_all[0] = dlat[0]
        dlat_all[1:-1] = (dlat[:-1] + dlat[1:]) / 2.0
        dlat_all[-1] = dlat[-1]
    else:
        dlat_all = np.array([np.pi])  # Full sphere

    # For pressure: use cell-centered approach
    if nlev > 1:
        dp = np.abs(np.diff(p_full))  # Absolute difference
        # Extend to get width for all cells
        dp_all = np.zeros(nlev)
        dp_all[0] = dp[0]
        dp_all[1:-1] = (dp[:-1] + dp[1:]) / 2.0
        dp_all[-1] = dp[-1]
    else:
        dp_all = np.array([p_full[0]])

    # Sum heating * dlat * dp
    total = 0.0
    for j in range(nlat):
        for k in range(nlev):
            total += heating_per_sec[j, k] * dlat_all[j] * dp_all[k]

    return total

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
    sigma = ((k - 0.5) / num_levels) ** exponent

    # Reference surface pressure
    p_surf = 1.0e5  # Pa

    # Pressure at full levels
    p_full = p_surf * np.exp(-scale_heights * (1 - sigma[::-1]))

    # Height using hypsometric equation
    H_scale = 7000.0  # meters (scale height)
    z_full = H_scale * np.log(p_surf / p_full)

    return p_full, z_full

def main():
    # Parameters
    h_amp = 0.1 / 86400.0  # K/s (0.1 K/day)
    p_s = 100000.0  # Pa
    p_width = 0.0175  # dimensionless
    lat_width = 0.4  # radians (~23 degrees)

    # Define all center positions (latitude in degrees, pressure in hPa)
    positions = [
        (0, 50),
        (0, 30),
        (15, 30),
        (-15, 30),
        (15, 50),
        (-15, 50),
        (30, 30),
        (-30, 30),
        (30, 50),
        (-30, 50),
        (45, 30),
        (-45, 30),
        (45, 50),
        (-45, 50),
        (45, 70),
        (-45, 70),
        (60, 30),
        (-60, 30),
        (60, 50),
        (-60, 50),
        (60, 70),
        (-60, 70),
        (75, 30),
        (-75, 30),
        (75, 50),
        (-75, 50),
        (75, 70),
        (-75, 70),
        (90, 30),
        (-90, 30),
        (90, 50),
        (-90, 50),
        (90, 70),
        (-90, 70),
    ]

    # Create grids
    num_levels = 50
    nlat = 128  # Higher resolution for smoother contours

    # Latitude grid
    lat = np.linspace(-np.pi/2, np.pi/2, nlat)
    lat_deg = np.degrees(lat)

    # Pressure grid
    p_full, z_full = create_pressure_height_grid(num_levels=num_levels)
    p_hPa = p_full / 100.0

    # Pad pressure grid
    p_extended_hPa = np.concatenate([
        np.array([1000.0]),
        p_hPa,
        np.array([1.0])
    ])

    # Native model-resolution latitude grid, for honest pixel rendering (pcolormesh) of the
    # example footprints instead of smooth contours that overstate the grid resolution.
    from matplotlib.colors import LinearSegmentedColormap
    nlat_model = 64
    lat_model_rad = np.linspace(-np.pi / 2, np.pi / 2, nlat_model)
    lat_model_deg = np.degrees(lat_model_rad)

    # Create figure
    plt.rcParams.update({
        'font.size': 12,
        'axes.titlesize': 14,
        'axes.labelsize': 12,
        'xtick.labelsize': 10,
        'ytick.labelsize': 10,
    })

    fig, ax = plt.subplots(figsize=(10, 6))

    # Alternating black/gray contours per pressure row (instead of distinct colors).
    # Rows are the pressure centers, ordered top->bottom (lowest pressure first). Within a
    # row, both outermost (highest |lat|) contours take the row's base color and the color
    # alternates moving inward toward the equator, so the pattern is symmetric about the
    # equator. The base (outermost) color alternates by row: black, gray, black, ...
    GREEN = (41/255, 94/255, 17/255)    # rgb(41,94,17)
    PURPLE = (88/255, 9/255, 79/255)    # rgb(88,9,79)
    GRAY = (0.22, 0.22, 0.22)           # ~same lightness as GREEN, to start

    # Only these example cases get a shaded footprint and a colored dot; every other center is
    # just a black dot. Maps position (lat_deg, p_hPa) -> shade color.
    HIGHLIGHT_COLOR = {(-90, 30): GREEN, (0, 50): PURPLE, (60, 70): GRAY}

    # Compute analytical integral (same for all positions)
    analytical_int = analytical_integral(h_amp, lat_width, p_width, p_s)
    print(f"Analytical integral: {analytical_int:.6e} K Pa rad / s")
    print("\nNormalization factors for each position:")

    # Plot contours for each position
    for i, (lat_center_deg, p_center_hPa) in enumerate(positions):
        lat_center_rad = np.radians(lat_center_deg)
        p_center = p_center_hPa * 100.0  # Convert to Pa

        # Compute heating (in K/day)
        heating_per_day = compute_ewa_heating(lat, p_full, h_amp, p_s, p_center,
                                              p_width, lat_center_rad, lat_width)

        # Convert to K/s for integration
        heating_per_sec = heating_per_day / 86400.0

        # Compute numerical integral
        numerical_int = numerical_integral(heating_per_sec, lat, p_full)

        # Compute normalization factor
        normalization = analytical_int / numerical_int

        # Print normalization info
        print(f"  Position ({lat_center_deg:+4.0f}°, {p_center_hPa:3.0f} hPa): "
              f"numerical_int = {numerical_int:.6e}, normalization = {normalization:.4f}")

        # Apply normalization to heating (in K/day)
        heating = heating_per_day * normalization

        # Pad heating array
        heating_padded = np.zeros((nlat, num_levels + 2))
        heating_padded[:, 1:-1] = heating

        # Five contour levels per patch, from 0.75x max (widest) to 0.90x max (smallest,
        # near the Gaussian center). Each level is drawn as a translucent filled region
        # (value >= level) plus a thin outline. The five fills nest and stack, so the middle
        # of each Gaussian is darker where all five overlap.
        pos = (lat_center_deg, p_center_hPa)

        # Example cases get a footprint at the model's NATIVE grid resolution (pcolormesh
        # pixels) so the coarseness is honest -- no smoothing implied -- and a colored dot.
        if pos in HIGHLIGHT_COLOR:
            color = HIGHLIGHT_COLOR[pos]
            heat_native = compute_ewa_heating(lat_model_rad, p_full, h_amp, p_s, p_center,
                                              p_width, lat_center_rad, lat_width)
            cmap = LinearSegmentedColormap.from_list(
                'shade', [(color[0], color[1], color[2], 0.0),
                          (color[0], color[1], color[2], 0.85)])
            ax.pcolormesh(lat_model_deg, p_hPa, heat_native.T,
                          cmap=cmap, vmin=0.0, vmax=np.max(heat_native),
                          shading='nearest', zorder=1)
            dot_color = color
        else:
            dot_color = 'black'

        # Center dot (drawn on top of any footprint).
        ax.plot(lat_center_deg, p_center_hPa, marker='o', color=dot_color,
                markersize=4, alpha=0.9, zorder=3)

        # Add label for this position
        label = f'{lat_center_deg:+.0f}°, {p_center_hPa:.0f} hPa'
        #cs.collections[0].set_label(label)

    # Set plot properties
    ax.set_yscale('log')
    ax.invert_yaxis()
    ax.set_xlim(-90, 90)
    ax.set_ylim(200, 2)
    ax.set_xticks([-90, -60, -30, 0, 30, 60, 90])
    ax.set_xticklabels(['90°S', '60°S', '30°S', '0°', '30°N', '60°N', '90°N'])
    ax.set_yticks([100, 10])
    ax.set_yticklabels(['100', '10'])
    ax.set_xlabel('Latitude')
    ax.set_ylabel('Pressure (hPa)')
    ax.set_title('SAI Heating Locations (centers; example footprints at model resolution)')
    ax.grid(True, alpha=0.3, which='both')

    plt.tight_layout()

    # Save figure
    output_name = 'plot_heating_location_grid'
    plt.savefig(output_name + '.png', dpi=300, bbox_inches='tight')
    plt.savefig(output_name + '.pdf', bbox_inches='tight')
    print(f"Saved plot to: {output_name}.png and {output_name}.pdf")

    # Also show the plot
    plt.show()

if __name__ == "__main__":
    main()
