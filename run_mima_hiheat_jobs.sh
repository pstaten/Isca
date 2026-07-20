#!/bin/bash
#SBATCH --job-name=mima_hiheat_runner
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --cpus-per-task=32
#SBATCH --time=48:00:00
#SBATCH -p general
#SBATCH -A r00132
#SBATCH -o hiheat_runner_%j.txt
#SBATCH -e hiheat_runner_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=pwstaten@iu.edu

# Load environment
source ~/.bashrc
conda activate isca_env

# Compile the codebase before launching, but ONLY if the build is missing. This
# self-heals a scratch-purged build (rebuilds isca.x + restores the mppnccombine.x
# symlink). Guarding on presence means jobs skip compiling when a good build already
# exists, so parallel jobs don't race on the shared build dir.
BUILD=$(ls -d $GFDL_WORK/codebase/*/build/isca 2>/dev/null | head -1)
if [ -z "$BUILD" ] || [ ! -e "$BUILD/isca.x" ] || [ ! -e "$BUILD/mppnccombine.x" ]; then
  echo "codebase build missing - compiling..."
  python -c "from isca import IscaCodeBase, GFDL_BASE; cb=IscaCodeBase.from_directory(GFDL_BASE); cb.compile()" || { echo "COMPILE FAILED - aborting"; exit 1; }
else
  echo "codebase build present - skipping compile"
fi

# Clear any INCOMPLETE run dirs (run{N} with no matching res{N}.tar.gz -- empty, 0-byte output,
# or full output but no restart; all poison that would make resume skip run{N} and die on the
# missing res{N}). Deleting lets resume re-run month N off the previous good restart. A run with
# a matching restart is never touched. Scoped to THIS job's experiments only.
for e in mima_hiheat0p1_qbo00 mima_hiheat0p1_qbo20; do
  for rd in "$GFDL_DATA/$e"/run[0-9]*; do
    [ -d "$rd" ] || continue
    n=$(basename "$rd"); n=${n#run}
    [ -f "$GFDL_DATA/$e/restarts/res${n}.tar.gz" ] || rm -rf "$rd"
  done
done

# Launch scripts in parallel
srun -n1 --cpus-per-task=32 python /N/slate/pwstaten/Projects/Isca/exp/MiMA/MiMA_hiheat0p1_qbo00.py &
srun -n1 --cpus-per-task=32 python /N/slate/pwstaten/Projects/Isca/exp/MiMA/MiMA_hiheat0p1_qbo20.py &

wait
