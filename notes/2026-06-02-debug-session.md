# Session notes — debugging the stuck QBO runs (2026-06-02)

First session back on the project after a break. Goal: get the stuck simulations un-stuck and
set up to analyze output. Outcome: **the failures are environment/workflow, not a physics
blow-up** — and we found the current blocker. No code changes made this session.

## What we set up
- Wrote/updated `AGENTS.md` (focused on the custom edits + run/debug/analysis workflow).
- Established agent access to the cluster via an SSH ControlMaster socket:
  - User runs once in a real terminal (Duo prompt):
    `ssh -M -S ~/.ssh/cm-quartz -o ControlPersist=8h -o IdentitiesOnly=yes quartz`
    then `unset TMOUT` (and optionally `tmux`) to keep it alive.
  - Agent reuses it read-only: `ssh -S ~/.ssh/cm-quartz quartz '…'`.
- Memory file updated: `current-debug-goal.md`.

## Cluster facts (settled)
- Production runs executed on **BigRed200** (Cray EX; nodes `nid####` / `x1002c6s5b1n1`;
  `GFDL_WORK=/N/scratch/pwstaten/isca_work_bigred`).
- Current login/`.bashrc` is configured for **Quartz** (`h1.quartz.uits.iu.edu`,
  `GFDL_ENV=quartz_conda`, `isca_work_quartz`). → a BigRed200→Quartz migration was underway.
- Paths: code on slate (`/N/slate/pwstaten/Projects/Isca`), output on scratch
  (`/N/scratch/pwstaten/isca_data`), post on `/N/project/pfec_staten/isca_data`.

## Key finding: NOT a NaN/blow-up
Grepped every SLURM `.err` for `NaN|valid range|FATAL|forrtl|backtrace` → **0 hits**.
Three real failure modes, none physics:
1. **`mppnccombine.x: No such file or directory` → ErrorReturnCode_127** — the post-processing
   tile-combiner is missing from the build dir, so every segment fails at the combine step.
   This is the **current primary blocker**. Compiled copy exists in repo `postprocessing/`.
   Leading cause: scratch auto-purge of the build dir. (newest `mima_runner` 6350239 ran 20 h
   then this; newest `hiheat_runner` 6336791 same.)
2. **OpenMPI `orte_init`/`orte_session_dir` mkdir `/tmp` failure** — older hiheat 6310285
   (Feb 4), died in 1 s. conda-OpenMPI-on-Cray TMPDIR issue.
3. **"Data for run N already exists, overwrite_data False. Stopping."** — benign resume skip
   (the `mima_noqbo_batch*` jobs).

Run counts reached ~530–584 of 1428 requested, **independent of forcing strength** → looks
like **48 h wall-clock truncation**, not instability.

## Open questions / decisions for next time
- Does the user actually remember seeing NaN/"valid range" text, or just "jobs kept dying"?
  (Leaning: misread `FailedRunError` tracebacks as crashes; the crash-mitigation namelist
  knobs are scars from an earlier, already-fixed crash.)
- **Resume on BigRed200** (rebuild there, keep existing `isca_data`) **or restart on Quartz**
  (fresh build + work dir)? Affects everything downstream.

## Next steps (when we resume)
1. Decide BigRed200-vs-Quartz (above).
2. Fix `mppnccombine.x`: rebuild it into the active `$GFDL_WORK/codebase/.../build/isca/`
   (or confirm Isca rebuilds it on a clean codebase compile). `postprocessing/compile_mppn.sh`
   is the standalone compile recipe.
3. Sanity-check one existing `atmos_monthly.nc` (temps/winds) to confirm the completed data is
   physically fine — rules out a real blow-up definitively.
4. (Optional) confirm a run's `fms.out`/`logfile.*.out` shows clean integration.
5. Then resume the run loop and move to analysis (EP fluxes, QBO descent, EWA response).
