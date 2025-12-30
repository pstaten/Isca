#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Plot proposed heating forcings at multiple latitudes and pressure levels.

This script plots contours of the ewa_heating forcing centered at various
combinations of latitude and pressure level.
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

    # Create figure
    plt.rcParams.update({
        'font.size': 12,
        'axes.titlesize': 14,
        'axes.labelsize': 12,
        'xtick.labelsize': 10,
        'ytick.labelsize': 10,
    })

    fig, ax = plt.subplots(figsize=(10, 6))

    # Define colormap for different positions
    # Use a set of distinct colors
    colors = plt.cm.tab20(np.linspace(0, 1, 20))

    # Plot contours for each position
    for i, (lat_center_deg, p_center_hPa) in enumerate(positions):
        lat_center_rad = np.radians(lat_center_deg)
        p_center = p_center_hPa * 100.0  # Convert to Pa

        # Compute heating
        heating = compute_ewa_heating(lat, p_full, h_amp, p_s, p_center,
                                      p_width, lat_center_rad, lat_width)

        # Pad heating array
        heating_padded = np.zeros((nlat, num_levels + 2))
        heating_padded[:, 1:-1] = heating

        # Choose contour level (e.g., 50% of maximum)
        max_heating = np.max(heating)
        contour_level = 0.75 * max_heating

        # Plot single contour
        color_idx = i % len(colors)
        cs = ax.contour(lat_deg, p_extended_hPa, heating_padded.T,
                       levels=[contour_level],
                       colors=[colors[color_idx]],
                       linewidths=1.5,
                       alpha=0.8)

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
    ax.set_title('Proposed Heating Forcing Positions (0.075 K/day contour)')
    ax.grid(True, alpha=0.3, which='both')

    # Add legend
    # Split legend into multiple columns to fit
    ax.legend(loc='upper left', ncol=2, fontsize=8, framealpha=0.9)

    plt.tight_layout()

    # Save figure
    output_name = 'plot_proposed_forcings'
    plt.savefig(output_name + '.png', dpi=300, bbox_inches='tight')
    plt.savefig(output_name + '.pdf', bbox_inches='tight')
    print(f"Saved plot to: {output_name}.png and {output_name}.pdf")

    # Also show the plot
    plt.show()

if __name__ == "__main__":
    main()
