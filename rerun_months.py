#!/usr/bin/env python3
"""Re-run specific bad/missing months for Isca experiments — the 'repair' half of the
validate-and-repair loop (detector: idealized_qbo_sai_analysis/validate_model_output.py).

For each flagged month it regenerates run<NNN>/atmos_monthly.nc via exp.run(N,
overwrite_data=True) off res(N-1), and removes the stale zonal-mean file so the analysis
pipeline recomputes it. Because the model is deterministic, re-running a month reproduces an
identical res(N), so downstream months stay valid — only the bad month is redone.

Each experiment's fully-configured Experiment object is rebuilt by exec'ing its exp/*/*.py with
Experiment.run temporarily neutralized, so the module constructs `exp` WITHOUT running its
1428-month loop. No edits to the exp scripts are needed.

Requirements / safety:
  * MUST run inside a SLURM allocation — the model launches via srun (like the normal launchers).
  * Default is DRY-RUN (finds the exp script, rebuilds `exp`, prints the plan; runs no model).
    Pass --apply to actually re-run.
  * Do NOT --apply months of an experiment whose run job is currently active (dir conflicts).

Usage:
  rerun_months.py --worklist bad_months.txt [--apply]
  rerun_months.py --exp mima_heat0p1_latcent00b_pcent030_noqbo --runs 843,844 [--apply]
"""
import argparse, glob, os, sys

GFDL_BASE = os.environ.get("GFDL_BASE")
GFDL_DATA = os.environ.get("GFDL_DATA")


def find_exp_script(dataname):
    """Map a data-dir name (lowercase) to its exp/*/*.py (basename lowercased == data name)."""
    for py in glob.glob(os.path.join(GFDL_BASE, "exp", "*", "*.py")):
        if os.path.basename(py)[:-3].lower() == dataname.lower():
            return py
    return None


def load_experiment(script):
    """exec the exp script to rebuild its `exp` object, with all destructive Experiment
    side effects neutralized, so this has NO filesystem impact — critical because the target
    experiment may be running and shares its work rundir by name. We stub out run() (the
    module-level 1428-month loop) AND clear_rundir() (which the constructor calls to EMPTY the
    live work directory). Both are restored before any real re-run."""
    import isca
    saved = {name: getattr(isca.Experiment, name) for name in ("run", "clear_rundir")}
    isca.Experiment.run = lambda self, *a, **k: None
    isca.Experiment.clear_rundir = lambda self, *a, **k: None
    ns = {"__name__": "__rerun__", "__file__": script}
    try:
        with open(script) as fh:
            exec(compile(fh.read(), script, "exec"), ns)
    finally:
        for name, fn in saved.items():
            setattr(isca.Experiment, name, fn)
    return ns.get("exp"), int(ns.get("NCORES", 32))


def experiment_busy(dataname, minutes=15):
    """Heuristic guard: True if the experiment's shared work rundir was touched recently,
    i.e. a run job is probably active — re-running it now would corrupt the live segment."""
    import time
    rd = os.path.join(os.environ.get("GFDL_WORK", ""), "experiment", dataname, "run")
    try:
        return (time.time() - os.path.getmtime(rd)) < minutes * 60
    except OSError:
        return False


def parse_worklist(path):
    jobs = {}
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t") if "\t" in line else line.split()
            if len(parts) < 2:
                continue
            jobs.setdefault(parts[0], set()).add(int(parts[1]))
    return {k: sorted(v) for k, v in jobs.items()}


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--worklist", help="worklist file from validate_model_output.py")
    ap.add_argument("--exp", help="single experiment data-dir name")
    ap.add_argument("--runs", help="comma-separated run numbers (with --exp)")
    ap.add_argument("--apply", action="store_true", help="actually re-run (default: dry-run)")
    ap.add_argument("--force", action="store_true", help="re-run even if the experiment looks active")
    ap.add_argument("--keep-zm", action="store_true", help="don't delete stale atmos_monthly_zm.nc")
    args = ap.parse_args()

    if not GFDL_BASE or not GFDL_DATA:
        sys.exit("ERROR: GFDL_BASE and GFDL_DATA must be set")
    if args.worklist:
        jobs = parse_worklist(args.worklist)
    elif args.exp and args.runs:
        jobs = {args.exp: sorted(int(x) for x in args.runs.split(","))}
    else:
        sys.exit("ERROR: provide --worklist FILE, or --exp NAME --runs N,N,...")

    mode = "APPLY" if args.apply else "DRY-RUN"
    n_months = sum(len(v) for v in jobs.values())
    print(f"[{mode}] {n_months} month(s) across {len(jobs)} experiment(s)\n")

    failures = 0
    for dataname, runs in jobs.items():
        script = find_exp_script(dataname)
        if not script:
            print(f"  {dataname}: NO exp script found — SKIP"); failures += 1; continue
        # Rebuild the Experiment (also validates the exec mechanism in dry-run).
        try:
            exp, ncores = load_experiment(script)
        except Exception as e:
            print(f"  {dataname}: exec failed ({type(e).__name__}: {e}) — SKIP"); failures += 1; continue
        if exp is None:
            print(f"  {dataname}: no `exp` object after exec — SKIP"); failures += 1; continue
        print(f"  {dataname}: {len(runs)} month(s) {runs}")
        print(f"      script={os.path.relpath(script, GFDL_BASE)}  exp={exp.name!r}  ncores={ncores}")
        if not args.apply:
            continue
        if experiment_busy(dataname) and not args.force:
            print(f"      SKIP: work rundir touched <15 min ago — experiment looks ACTIVE. "
                  f"Wait for its job to finish, or pass --force."); failures += 1; continue
        for n in runs:
            if not args.keep_zm:
                zm = os.path.join(GFDL_DATA, dataname, f"run{n:04d}", "atmos_monthly_zm.nc")
                if os.path.isfile(zm):
                    os.remove(zm); print(f"      run{n:04d}: removed stale atmos_monthly_zm.nc")
            print(f"      run{n:04d}: re-running (overwrite_data=True) ...")
            exp.run(n, overwrite_data=True, num_cores=ncores)

    print(f"\ndone ({failures} experiment(s) skipped)")
    sys.exit(1 if failures else 0)


if __name__ == "__main__":
    main()
