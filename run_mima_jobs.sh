#!/bin/bash
#SBATCH --job-name=mima_runner
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=32
#SBATCH --time=48:00:00
#SBATCH -p general
#SBATCH -A r00132
#SBATCH -o mima_runner_%j.txt
#SBATCH -e mima_runner_%j.err
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

# Launch scripts in parallel
srun -n1 --cpus-per-task=32 python /N/slate/pwstaten/Projects/Isca/exp/MiMA/MiMA_heat0p0_qbo00.py &
srun -n1 --cpus-per-task=32 python /N/slate/pwstaten/Projects/Isca/exp/MiMA/MiMA_heat0p0_qbo20.py &
srun -n1 --cpus-per-task=32 python /N/slate/pwstaten/Projects/Isca/exp/MiMA/MiMA_heat0p1_qbo00.py &
srun -n1 --cpus-per-task=32 python /N/slate/pwstaten/Projects/Isca/exp/MiMA/MiMA_heat0p1_qbo20.py &

wait
