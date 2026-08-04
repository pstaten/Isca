#!/bin/bash
#SBATCH --job-name=mima_fixh2o
#SBATCH --nodes=1
#SBATCH --ntasks=3
#SBATCH --cpus-per-task=32
#SBATCH --time=48:00:00
#SBATCH -p general
#SBATCH -A r00132
#SBATCH -o mima_fixh2o_%j.txt
#SBATCH -e mima_fixh2o_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=pwstaten@iu.edu

# Mechanism-denial runs for Ewa's strat-H2O question: three slab experiments with
# RRTM fed FIXED stratospheric water vapor (do_fixed_water, p < 100 hPa), so the
# heating-induced H2O increase cannot act radiatively. Perturbation - control with
# radiation blind to H2O, vs the standard pairs, isolates the H2O greenhouse
# contribution to the slab surface warming.

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
for e in mima_heat0p0_noqbo_fixh2o mima_heat0p1_latcent00b_pcent050_noqbo_fixh2o mima_heat0p1_latcent90b_pcent050_noqbo_fixh2o; do
  for rd in "$GFDL_DATA/$e"/run[0-9]*; do
    [ -d "$rd" ] || continue
    n=$(basename "$rd"); n=${n#run}
    [ -f "$GFDL_DATA/$e/restarts/res${n}.tar.gz" ] || rm -rf "$rd"
  done
done

srun -n1 --cpus-per-task=32 python /N/slate/pwstaten/Projects/Isca/exp/MiMA/MiMA_heat0p0_noqbo_fixh2o.py &
srun -n1 --cpus-per-task=32 python /N/slate/pwstaten/Projects/Isca/exp/MiMA/MiMA_heat0p1_latcent00B_pcent050_noqbo_fixh2o.py &
srun -n1 --cpus-per-task=32 python /N/slate/pwstaten/Projects/Isca/exp/MiMA/MiMA_heat0p1_latcent90B_pcent050_noqbo_fixh2o.py &

wait
