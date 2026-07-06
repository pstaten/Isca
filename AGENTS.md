# AGENTS.md — QBO × SAI-heating fork of Isca

Working notes for Peter Staten (IU) and AI agents. This is a fork of **Isca** (the
GFDL-based idealized GCM framework) modified to study **how stratospheric-aerosol-injection
(SAI)-like heating modulates the impacts of the QBO**. Collaboration with **Ben Kravitz**
and **Ewa Bednarz**.

> **Agents: read this whole file at the start of every session** to get up to speed before
> doing anything. It documents the *scientific plan*, *the edits in this fork*, and the
> *run / debug / analysis workflow* — not upstream Isca internals. For stock Isca, see
> `ReadMe.md`, `docs/`, and https://execlim.github.io/Isca.

---

## Scientific goal & experimental design

Run idealized GCM experiments **with and without an imposed QBO-like stratospheric zonal-wind
signal**, while **systematically applying SAI-like heating at different locations throughout
the lower and middle stratosphere**, and ask how the heating changes the QBO's influence on
climate.

**Three tiers of complexity** (not every experiment is run in all three):
- **Held–Suarez (`hs`)** — simplest; dry dynamical core with Newtonian relaxation.
- **Frierson (`frierson25lev` / `frierson50lev`)** — grey-radiation moist aquaplanet.
- **MiMA (`MiMA`)** — most complete (RRTM radiation). The main production sweep lives here.

**The two imposed forcings** (both in `damping_driver.f90`; in Held–Suarez they live in the
HS forcing code instead — see below):
- **QBO-like nudging** — relaxes equatorial `u` toward a downward-propagating sinusoid.
- **SAI-like heating** — a localized Gaussian stratospheric heating (the "EWA" heating, after
  Bednarz's profile), swept across latitude and pressure.

**Core comparisons the experiments are designed to support:**
1. **QBO impact.** `qbo20 − qbo00` = (20 m/s QBO-like forcing) minus (0 m/s nudging, i.e. `u`
   nudged toward *zero* stratospheric wind). Holding the nudging machinery identical in both
   isolates the QBO's impact. ⚠️ **`qbo00` ≠ `noqbo`**: `qbo00` keeps the nudging ON but at
   amplitude 0 (actively pins `u→0`); `noqbo` turns nudging OFF (free dynamics). The clean
   "with vs without QBO" contrast is **`qbo20 − qbo00`**.
2. **SAI modulation of the QBO impact.** Repeat that difference under different SAI-heating
   configurations (heating amplitude, and especially heating *location* in lat × pressure) to
   see how SAI-like heating changes the QBO impact.
3. **Nonlinear (rectified) QBO impact.** Over an integer number of complete QBO cycles the
   *time-mean forcing/nudging* is ~zero, but the *time-mean impact* of the QBO need not be —
   it may be small but nonzero. Quantifying this rectified response requires **running long
   enough to cover complete QBO cycles** (period = 28 months; production runs do ~84-month
   spinup + several 28-month cycles). This is a primary science target, so protect run length.

---

## The scientific edits — all in `damping_driver.f90`

Everything custom lives in **`src/atmos_param/damping_driver/damping_driver.f90`**
(Held–Suarez configs carry the equivalent forcing in the HS code). Two forcings, both driven
from `&damping_driver_nml`, both with diagnostics. Namelist decl + defaults:
`damping_driver.f90:57–80`.

### 1. `fake_qbo` — nudged equatorial QBO  (`do_sin_qbo = .true.`)
- Relaxes zonal wind `u` toward a downward-propagating sinusoid
  `qbo = qbo_amp · cos(m·z − ν·t)` (`:861–864`): 28-month period, 40 km vertical wavelength.
- Confined to the equator (±12°, 2° tanh edges) and to ~8–80 hPa (tanh windows in `pfull`).
- 5-day relaxation; tendency `utnd = (qbo − u)·(1 − exp(−relax_rate))`, added to `udt`.
- Knob: **`qbo_amp`** (m/s; experiments use 0 or 20). Diagnostic: **`udt_qbo`** (group `damping`).

### 2. `ewa_heating` — localized SAI-like stratospheric heating  (`do_ewa_htg = .true.`)
- Gaussian in latitude × pressure (Bednarz profile, `:744–756`):
  `Q = h_amp · exp(−(p−p_center)²/2p_width²) · exp(−(lat−lat_center)²/2lat_width²)`.
- **Normalized** (`:762`): the discrete grid integral is rescaled to the analytic integral
  `2π·h_amp·|lat_width·p_width·p_s|`, so total injected heat is grid-independent.
  ⚠️ Normalization divides by the numerical integral — guard against degenerate grids.
- Knobs (all `&damping_driver_nml`): `h_amp` (K/s; scripts pass `0.1/86400.` etc. — the
  `heatXpY` in names is K/day), `p_center`, `p_width` (Pa), `p_s` (=100000), `lat_center`,
  `lat_width` (radians), `both_hemispheres` (.true. → symmetric via `abs(lat)`; `:747–752`).
- Diagnostic: **`tdt_ewa`** (K/s, group `damping`).

Apply sites (added to `udt`/`tdt`): `:343–361`. Diagnostics registered in
`damping_driver_init` (`id_udt_qbo`, `id_tdt_ewa` declared at `:108`).

---

## Experiment naming convention

Scripts live in `exp/{hs, frierson25lev, frierson50lev, MiMA}/`. Filenames encode params:

| Token        | Meaning                                                                    |
|--------------|-----------------------------------------------------------------------------|
| `heat0p1`    | `h_amp = 0.1` K/day (`heat0p0` = heating off); `hiheat` = larger `h_amp`     |
| `qbo20`      | `qbo_amp = 20` m/s; `qbo00` = nudging ON but amp 0 (u→0); `noqbo` = `do_sin_qbo=False` |
| `latcent30S` | `lat_center = 30°S`; `30B` = 30° **B**oth hemispheres; `00B` = equatorial    |
| `pcent050`   | `p_center = 50` hPa                                                          |

To add an experiment: copy the closest existing script, rename, set the `Experiment(...)` name
to match the filename (lowercase), edit `&damping_driver_nml`, keep the diag_table block intact.
The big no-QBO heating sweep is `MiMA_heat0p1_latcent*_pcent*_noqbo.py`, batched via
`submit_mima_noqbo_batches.sh`. The QBO/hiheat runs (`*_qbo00`, `*_qbo20`, `hiheat*`) go via
`run_mima_jobs.sh`.

---

## Running (IU HPC — BigRed200)

Runs execute on **BigRed200** (Cray EX; conda OpenMPI + `srun`). Env vars (user's Big Red
`~/.bashrc`, set before running/submitting):
- `GFDL_BASE=/N/slate/pwstaten/Projects/Isca` — source code (on **slate**; = repo root)
- `GFDL_ENV=bigred_conda` — selects compiler env file under `src/extra/env`
- `GFDL_WORK=/N/scratch/pwstaten/isca_work_bigred` — build + transient run scratch.
  ⚠️ **scratch AUTO-PURGES unaccessed files** — the build dir (incl. the `mppnccombine.x`
  symlink) vanishes while away. This is the root of the recurring failure (see below).
- `GFDL_DATA=/N/project/pfec_staten/isca_data` — model output (on **project**; does NOT purge)
- conda env `isca_env`. SLURM: account `r00132`, partitions `general` (48 h) / `debug` (1 h),
  1 node × 4 tasks × 32 cpus.

**Each exp script** does `cb = IscaCodeBase.from_directory(GFDL_BASE)` then loops
`exp.run(1, use_restart=True)` … `for i in range(2,1429): exp.run(i)` (30-day segments;
~84-month spinup + 4×28-month QBO cycles). **The scripts do NOT call `cb.compile()`** — they
assume a pre-built model in `$GFDL_WORK`. With `overwrite_data=False`, `run(i)` skips segments
whose `run{i}` dir already exists, then runs the first missing one, loading restart `res{i-1}`.

**SLURM launchers** (4 experiments/node in parallel via `srun`):
`run_mima_jobs.sh`, `run_frierson50lev_jobs.sh`, `run_hs_jobs.sh` (general, 48 h);
`debug_*_jobs.sh` (debug, 1 h); `submit_mima_noqbo_batches.sh` (auto-detects unfinished sweep
runs by checking for `run1428`, batches of 4). `cleanup_last_runs.sh` deletes the highest
`run*` per exp — **use with care** (see below), it does not check whether that run is complete.

**SSH access for agents** (BigRed200 needs Duo, so no headless login): user opens a
multiplexed master in a real terminal —
`ssh -M -S ~/.ssh/cm-bigred -o ControlPersist=8h -o IdentitiesOnly=yes bigred200` — then the
agent runs read-only commands via `ssh -S ~/.ssh/cm-bigred bigred200 '…'`. SLURM
`*_batch*_<jobid>.{err,txt}` logs live in `$GFDL_BASE` (slate); the model's own stdout
(`logfile.*.out`, `fms.out`) is in the per-run `$GFDL_WORK` dir (transient, on scratch).

---

## Current status & the recurring failure mode (diagnosed 2026-07-06)

The MiMA no-QBO sweep got **stuck at a different segment for each experiment** (1029, 1040,
1017, …). Fully diagnosed — **not physics, not blow-ups**. Causal chain:

1. Scripts never `compile()`; they rely on the `$GFDL_WORK` build, which **scratch purged**
   (`$GFDL_WORK` down to ~17 MB; `build/isca/` and its `mppnccombine.x` gone).
2. This fork's `compile.sh` doesn't build mppnccombine (the `cc` lines are commented) — it
   just **symlinks `$GFDL_BASE/postprocessing/mppnccombine.x`** into the build dir. Nothing
   recreates that symlink once purged.
3. `experiment.py:run()` order: model runs → log "Run N complete" → `mkdir(run{N})` (line 302)
   → combine via `mppnccombine.x` → write `restarts/res{N}.tar.gz` (line 337). Missing combiner
   → **`ErrorReturnCode_127` after the empty `run{N}` dir exists but before `res{N}`**.
4. **Feb 16 2026 10:00**: every running exp died this way, each at a different segment → leaves
   **empty `run{N}`, no `res{N}`**.
5. **Every resubmit since** skips the empty `run{N}` (overwrite_data=False) and asks for `res{N}`
   → `IOError: Restart file not found` → dies in seconds, zero progress. That is the "stuck" loop.

(An *earlier*, since-fixed failure ~Jan 30 was `Disk quota exceeded` at run ~534 — the old
"~530–584" cluster. Quota was raised to 120 TB. **But `/N/project` is now ~96% full (~5 TB
free), `pfec_staten` ~95%** — watch space when finishing the sweep.)

**Fix (as of 2026-07-06, plan; confirm before destructive steps):**
1. **Recompile once** → rebuilds `isca.x` and restores the `mppnccombine.x` symlink.
2. **Delete only the empty poison dirs** (top `run{N}` with 0 files and no `res{N}`). Audit
   2026-07-06: 38/41 exps qualify. **Do NOT use blanket `cleanup_last_runs.sh`** — 3 exps have a
   *complete* top run (`latcent45s_pcent050/run0724`, `hiheat0p1_qbo00/run1343`,
   `hiheat0p1_qbo20/run0338`) that it would wrongly delete. Those 3 resume fine as-is.
3. **Resubmit** (`submit_mima_noqbo_batches.sh` / `run_mima_jobs.sh`).
4. **Durable fix** so the purge can't re-trap us: **add `cb.compile()` to the exp scripts** (each
   job self-heals a purged build), or move `$GFDL_WORK` off auto-purge scratch.

Crash-mitigation knobs already in the MiMA config (scars of an old, since-fixed cold-strat
crash — revisit only if a genuine NaN reappears): `valid_range_t:[100,800]`, `num_levels:50`,
`scale_heights:6.0`, `damping_order:4`, `robert_coeff:0.03`; `do_rayleigh:True`,
`trayfric:-0.25` (¼-day sponge), `sponge_pbottom:1200.` (Pa).

---

## Analyzing output

Output is CF-style NetCDF: `$GFDL_DATA/<exp>/run<NNN>/atmos_monthly.nc` (30-day means);
restarts are `$GFDL_DATA/<exp>/restarts/res<NNN>.tar.gz`. Standard diag set (see any MiMA
script) by group:
- `dynamics`: `ps, ucomp, vcomp, temp, omega, height`, eddy fluxes (`ucomp_vcomp`,
  `vcomp_temp`, `vcomp_omega`, `ucomp_omega`, `ucomp_temp`, `omega_temp` → EP flux), variances
  (`*_sq`), `vor, div`, `pres_full/half`, `sphum` + moisture fluxes, `bk/pk`.
- `atmosphere`: moisture/temperature tendencies (`dt_qg_*`, `dt_tg_*`), `rh`.
- `hs_forcing`: `teq`, `tdt_ndamp`.  ·  `damping`: **`udt_qbo`, `tdt_ewa`** (the custom forcings).
- `rrtm_radiation`: `tdt_rad, tdt_sw, tdt_lw`.

Typical analysis: `qbo20 − qbo00` for the QBO impact (per SAI-heating config); time-mean over
complete 28-month QBO cycles for the nonlinear/rectified impact. Forcing-visualization scripts
(run locally): `plot_proposed_forcings.py`, `plot_imposed_forcings.py` → PNG/PDF (gitignored).

---

## Git sync — keep both checkouts + GitHub in step

Same remote/branch on both machines:

```
MacBook  /Users/pwstaten/Projects/Isca  ⟷  github.com:pstaten/Isca (master)  ⟷  BigRed  $GFDL_BASE
```

**Standing rule: any edits made in a session get committed and pushed, and both checkouts pulled
up to date before the session ends.** After editing (edits usually happen on the MacBook):
1. On the machine that was edited: `git add … && git commit && git push origin master`.
2. On the other machine, fast-forward it:
   - BigRed: `ssh -S ~/.ssh/cm-bigred bigred200 'bash -lc "cd \$GFDL_BASE && git pull --ff-only"'`
   - MacBook: `git pull --ff-only`
3. Confirm `git status` on both shows `up to date with origin/master`.

Note: BigRed's `$GFDL_BASE` carries some untracked local-only files (test-case scripts, `.swp`,
a nested `Isca/`) — harmless, they don't block a fast-forward pull. `$GFDL_ENV=bigred_conda` and
the env are set from `~/.bashrc` there, not the repo.

---

## Conventions for agents editing this repo

- **Read this file at the start of every session.**
- Keep all custom physics in `damping_driver.f90` (or the HS forcing code for Held–Suarez);
  don't scatter it. Match the existing Fortran style (4-space indent, namelist defaults near
  the top, diag fields registered in `damping_driver_init`).
- New namelist params must be added to **both** the declaration block **and** the
  `namelist /damping_driver_nml/` list (typo here = silent default).
- **Sync every session's edits (standing rule).** Any time we make edits in a session, commit
  them and keep **both** checkouts up to date with GitHub. Don't commit run output, logs, or
  `$GFDL_DATA` contents.
- This is a fork — keep upstream Isca files untouched unless the change is intentional and noted
  in the commit message.
- Destructive or outward-facing cluster actions (deleting run dirs, submitting jobs, recompiling)
  — confirm with Peter first.
