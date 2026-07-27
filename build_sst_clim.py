#!/usr/bin/env python3
"""Build a prescribed-SST climatology file from the control (mima_heat0p0_noqbo) restarts, for the
fixed-SST sensitivity experiments. Extracts t_surf from mixed_layer.res.nc across ~10 yr of monthly
restarts, averages by calendar month -> 12-month THIRTY_DAY_MONTHS climatology on the T42 grid, and
writes it via Isca's create_timeseries (FMS interpolator-compatible). Output: sst_ctrl_clim.nc with
variable 'sst_ctrl_clim'. Run in the qbo_sai HPC env on BigRed (needs GFDL_BASE, GFDL_DATA)."""
import os, sys, tarfile, tempfile
import numpy as np, xarray as xr

GFDL_BASE = os.environ["GFDL_BASE"]; GFDL_DATA = os.environ["GFDL_DATA"]
sys.path.insert(0, os.path.join(GFDL_BASE, "src/extra/python/scripts"))
import create_timeseries as cts

lons, lats, lonbs, latbs, nlon, nlat, nlonb, nlatb = cts.create_grid(False)  # T42 model grid
lats = np.asarray(lats)
print(f"grid: nlon={nlon} nlat={nlat}; lat[0]={lats[0]:.2f} lat[-1]={lats[-1]:.2f}")

restart_dir = os.path.join(GFDL_DATA, "mima_heat0p0_noqbo", "restarts")
N0, N1 = 1309, 1428  # last 10 years (120 monthly restarts)
acc = np.zeros((12, nlat, nlon)); cnt = np.zeros(12)
with tempfile.TemporaryDirectory() as td:
    for N in range(N0, N1 + 1):
        rp = os.path.join(restart_dir, f"res{N:04d}.tar.gz")
        if not os.path.isfile(rp):
            continue
        with tarfile.open(rp) as tf:
            m = [x for x in tf.getnames() if x.endswith("mixed_layer.res.nc")][0]
            tf.extract(m, path=td)
        ds = xr.open_dataset(os.path.join(td, m))
        ts = np.asarray(ds["t_surf"].values).squeeze()  # (nlat, nlon)
        ds.close()
        mon = (N - 1) % 12
        acc[mon] += ts; cnt[mon] += 1
sst = acc / cnt[:, None, None]
print(f"months filled: {int((cnt>0).sum())}/12; counts/month: {cnt.astype(int)}")

# orientation sanity check: annual+zonal mean SST must peak near the equator, not the poles
zm = sst.mean(axis=(0, 2))  # (nlat,)
kmax = int(np.argmax(zm))
print(f"warmest zonal-mean SST at lat={lats[kmax]:.1f} (should be ~0); tropics~{zm[np.argmin(np.abs(lats))]:.1f}K poles~{zm[0]:.1f}/{zm[-1]:.1f}K")
if abs(lats[kmax]) > 45:
    print("  -> SST peaks off-equator; latitude axis is flipped relative to the T42 grid -> reversing")
    sst = sst[:, ::-1, :]

time_arr, day_number, ntime, time_units, time_bounds = cts.create_time_arr(1, True, 12)
number_dict = {"nlat": nlat, "nlon": nlon, "nlatb": nlatb, "nlonb": nlonb}
# p_full=None -> writes a 2D (time,lat,lon) field (is_thd keyed off p_full being None)
cts.output_to_file(sst.astype("f4"), lats, lons, latbs, lonbs, None, None,
                   time_arr, time_units, "sst_ctrl_clim", "sst_ctrl_clim", number_dict,
                   time_bounds=time_bounds)
print("wrote sst_ctrl_clim.nc")
