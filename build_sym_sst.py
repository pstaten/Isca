#!/usr/bin/env python3
"""Symmetrize the prescribed SST climatology: zonal mean + aquaplanet R-fold.

The harvested clim (input/sst_ctrl_clim.nc) carries the slab control's bistable-state
signature (antisymmetric dipole peaking +/-1.2 K at +/-4 deg; NH-SH +0.59 K) and
~0.12 K RMS of frozen zonal weather noise. This script writes
input/sst_ctrl_clim_sym.nc with

    SST_sym(lat, lon, m) = zonal_mean( 0.5*(Z(lat, m) + Z(-lat, (m+6) mod 12)) )

i.e. zonally uniform and exactly invariant under the model's (lat -> -lat,
month -> month+6) symmetry -- the average of the two mirror bistable states, the
slab's unstable symmetric fixed point. File structure is preserved byte-for-byte
except the SST values (copy + in-place value replacement via netCDF4).
"""
import shutil

import netCDF4
import numpy as np

SRC = "input/sst_ctrl_clim.nc"
DST = "input/sst_ctrl_clim_sym.nc"

shutil.copy(SRC, DST)
ds = netCDF4.Dataset(DST, "r+")
name = next(v for v in ds.variables
            if ds.variables[v].ndim == 3 and v not in ("latb", "lonb"))
sst = ds.variables[name][:]                      # (time=12, lat, lon)
lat = ds.variables["lat"][:]
assert sst.shape[0] == 12, sst.shape

Z = sst.mean(axis=2)                             # zonal mean (12, lat)
# R-fold: average with the latitude mirror shifted by 6 months
Zm = Z[:, ::-1]                                  # lat -> -lat (grid symmetric)
Zs = 0.5 * (Z + np.roll(Zm, 6, axis=0))          # (m) with mirror at (m+6)
resid_asym = 0.5 * (Z - np.roll(Zm, 6, axis=0))
zonal_rms = float(np.sqrt(((sst - Z[:, :, None]) ** 2).mean()))
print(f"removed R-asymmetric component: max |A| = {np.abs(resid_asym).max():.3f} K, "
      f"rms = {np.sqrt((resid_asym**2).mean()):.3f} K")
print(f"removed zonal deviations: rms = {zonal_rms:.3f} K")

ds.variables[name][:] = np.repeat(Zs[:, :, None], sst.shape[2], axis=2)
ds.setncattr("history", "symmetrized (zonal mean + lat/month+6 R-fold) from sst_ctrl_clim.nc "
                        "by build_sym_sst.py")
ds.close()

# verification
ds = netCDF4.Dataset(DST)
s = ds.variables[name][:]
Z2 = s.mean(axis=2)
a = 0.5 * (Z2 - np.roll(Z2[:, ::-1], 6, axis=0))
print(f"verify: residual asymmetry {np.abs(a).max():.2e} K; "
      f"zonal std {float((s - Z2[:, :, None]).std()):.2e} K")
print(f"wrote {DST}")
