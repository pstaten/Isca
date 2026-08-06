#!/bin/bash
#SBATCH --job-name=mima_qbocan
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --cpus-per-task=32
#SBATCH --time=48:00:00
#SBATCH -p general
#SBATCH -A r00132
#SBATCH -o mima_qbocan_%j.txt
#SBATCH -e mima_qbocan_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=pwstaten@iu.edu

# Canonical observed-QBO fixed-SST pair (do_tab_qbo): heat0p0 and heat0p1 (SAI-like
# heating), mirroring the qbo20 pair; differenced against the existing qbo00 controls.
# REQUIRES the build compiled after the do_tab_qbo commit (the gate job handles that).

source ~/.bashrc
conda activate isca_env

# Self-heal a scratch-purged build (see run_mima_jobs.sh).
BUILD=$(ls -d $GFDL_WORK/codebase/*/build/isca 2>/dev/null | head -1)
if [ -z "$BUILD" ] || [ ! -e "$BUILD/isca.x" ] || [ ! -e "$BUILD/mppnccombine.x" ]; then
  echo "codebase build missing - compiling..."
  python -c "from isca import IscaCodeBase, GFDL_BASE; cb=IscaCodeBase.from_directory(GFDL_BASE); cb.compile()" || { echo "COMPILE FAILED - aborting"; exit 1; }
else
  echo "codebase build present - skipping compile"
fi

# Clear incomplete (restart-less) run dirs so resume re-runs them (see run_mima_jobs.sh).
for e in mima_heat0p0_qbocan_fixsst mima_heat0p1_qbocan_fixsst; do
  for rd in "$GFDL_DATA/$e"/run[0-9]*; do
    [ -d "$rd" ] || continue
    n=$(basename "$rd"); n=${n#run}
    [ -f "$GFDL_DATA/$e/restarts/res${n}.tar.gz" ] || rm -rf "$rd"
  done
done

srun -n1 --cpus-per-task=32 python /N/slate/pwstaten/Projects/Isca/exp/MiMA/MiMA_heat0p0_qbocan_fixsst.py &
srun -n1 --cpus-per-task=32 python /N/slate/pwstaten/Projects/Isca/exp/MiMA/MiMA_heat0p1_qbocan_fixsst.py &

wait
