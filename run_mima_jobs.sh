
#!/bin/bash
#SBATCH --job-name=mima_runner
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=32
#SBATCH --time=06:00:00
#SBATCH -p general
#SBATCH -A r00132
#SBATCH -o mima_runner_%j.txt
#SBATCH -e mima_runner_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=pwstaten@iu.edu

# Load any required modules
module load isca_env

# Launch scripts in parallel
srun -n1 /N/slate/pwstaten/Projects/Isca/exp/MiMA/pwdMiMA_heat0p0_qbo00.py &
srun -n1 /N/slate/pwstaten/Projects/Isca/exp/MiMA/pwdMiMA_heat0p0_qbo20.py &
srun -n1 /N/slate/pwstaten/Projects/Isca/exp/MiMA/pwdMiMA_heat0p1_qbo00.py &
srun -n1 /N/slate/pwstaten/Projects/Isca/exp/MiMA/pwdMiMA_heat0p1_qbo20.py &

wait
