#!/usr/bin/env python3
"""Generate fixed-SST versions of every exp/MiMA/MiMA_*.py.  Each _fixsst copy:
 - renames the Experiment (adds '_fixsst' -> new data dir)
 - adds the prescribed-SST climatology to exp.inputfiles
 - switches mixed_layer_nml to prescribed SST (do_read_sst/do_sc_sst/sst_file='sst_ctrl_clim')
 - outputs t_surf (to verify the SST is held)
Everything else (heating config, QBO config, diagnostics, 1428-month run loop) is identical, so the
only difference vs the interactive runs is that the surface is prescribed to the control climatology.
Run from the Isca repo root."""
import glob, os, re, sys

files = sorted(f for f in glob.glob("exp/MiMA/MiMA_*.py") if "_fixsst" not in f)
n = 0
for f in files:
    txt = open(f).read()
    m = re.search(r"Experiment\('(mima_[^']+)'", txt)
    if not m:
        print(f"SKIP (no Experiment name): {f}"); continue
    name = m.group(1); new = name + "_fixsst"
    orig = txt

    txt = txt.replace(f"Experiment('{name}'", f"Experiment('{new}'")

    txt = txt.replace(
        "exp.inputfiles = [os.path.join(GFDL_BASE,'input/rrtm_input_files/ozone_1990.nc')]",
        "exp.inputfiles = [os.path.join(GFDL_BASE,'input/rrtm_input_files/ozone_1990.nc'), "
        "os.path.join(GFDL_BASE,'input/sst_ctrl_clim.nc')]")

    txt = re.sub(
        r"'do_qflux': True\s*\n(\s*)\},",
        "'do_qflux': True,\n"
        "        'do_read_sst': True,   # fixed SST: prescribe the control's own climatology\n"
        "        'do_sc_sst': True,\n"
        "        'sst_file': 'sst_ctrl_clim',\n"
        r"\1},", txt, count=1)

    txt = txt.replace(
        "exp.diag_table = diag",
        "diag.add_field('mixed_layer', 't_surf', time_avg=True)  # verify prescribed SST is held\n"
        "exp.diag_table = diag", 1)

    # sanity: all four edits must have fired
    checks = {"name": new in txt, "inputfile": "sst_ctrl_clim.nc" in txt,
              "namelist": "'sst_file': 'sst_ctrl_clim'" in txt, "t_surf": "'t_surf'" in txt}
    if not all(checks.values()) or txt == orig:
        print(f"FAILED edits on {f}: {checks}"); sys.exit(1)

    out = f[:-3] + "_fixsst.py"
    open(out, "w").write(txt); n += 1
print(f"generated {n} fixed-SST scripts")
