#!/bin/bash
#SBATCH --job-name=frierson50lev_debugger
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=32
#SBATCH --time=01:00:00
#SBATCH -p debug
#SBATCH -A r00132
#SBATCH -o frierson50lev_debugger_%j.txt
#SBATCH -e frierson50lev_debugger_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=pwstaten@iu.edu

# Load environment
source ~/.bashrc
conda activate isca_env

# Launch scripts in parallel
srun -n1 --cpus-per-task=32 python /N/slate/pwstaten/Projects/Isca/exp/frierson50lev/frierson50lev_heat0p0_qbo00.py &
srun -n1 --cpus-per-task=32 python /N/slate/pwstaten/Projects/Isca/exp/frierson50lev/frierson50lev_heat0p0_qbo20.py &
srun -n1 --cpus-per-task=32 python /N/slate/pwstaten/Projects/Isca/exp/frierson50lev/frierson50lev_heat0p1_qbo00.py &
srun -n1 --cpus-per-task=32 python /N/slate/pwstaten/Projects/Isca/exp/frierson50lev/frierson50lev_heat0p1_qbo20.py &

wait
