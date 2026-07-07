#!/bin/bash
#
# Submit MiMA no-QBO experiments in batches of 4 to SLURM
# This script automatically detects unfinished experiments and creates batch job scripts
#

# Configuration
BATCH_SIZE=4
EXP_DIR="/N/slate/pwstaten/Projects/Isca/exp/MiMA"
SCRIPT_DIR="/N/slate/pwstaten/Projects/Isca"

# Check if GFDL_DATA is set
if [ -z "$GFDL_DATA" ]; then
    echo "ERROR: GFDL_DATA environment variable is not set"
    echo "Please set GFDL_DATA to your data directory"
    exit 1
fi

echo "GFDL_DATA is set to: $GFDL_DATA"
echo ""

# Define all experiment names (without .py extension)
declare -a all_experiments=(
    "MiMA_heat0p1_latcent00B_pcent030_noqbo"
    "MiMA_heat0p1_latcent00B_pcent050_noqbo"
    "MiMA_heat0p1_latcent15B_pcent030_noqbo"
    "MiMA_heat0p1_latcent15B_pcent050_noqbo"
    "MiMA_heat0p1_latcent15S_pcent030_noqbo"
    "MiMA_heat0p1_latcent15S_pcent050_noqbo"
    "MiMA_heat0p1_latcent30B_pcent030_noqbo"
    "MiMA_heat0p1_latcent30B_pcent050_noqbo"
    "MiMA_heat0p1_latcent30S_pcent030_noqbo"
    "MiMA_heat0p1_latcent30S_pcent050_noqbo"
    "MiMA_heat0p1_latcent45B_pcent030_noqbo"
    "MiMA_heat0p1_latcent45B_pcent050_noqbo"
    "MiMA_heat0p1_latcent45B_pcent070_noqbo"
    "MiMA_heat0p1_latcent45S_pcent030_noqbo"
    "MiMA_heat0p1_latcent45S_pcent050_noqbo"
    "MiMA_heat0p1_latcent45S_pcent070_noqbo"
    "MiMA_heat0p1_latcent60B_pcent030_noqbo"
    "MiMA_heat0p1_latcent60B_pcent050_noqbo"
    "MiMA_heat0p1_latcent60B_pcent070_noqbo"
    "MiMA_heat0p1_latcent60S_pcent030_noqbo"
    "MiMA_heat0p1_latcent60S_pcent050_noqbo"
    "MiMA_heat0p1_latcent60S_pcent070_noqbo"
    "MiMA_heat0p1_latcent75B_pcent030_noqbo"
    "MiMA_heat0p1_latcent75B_pcent050_noqbo"
    "MiMA_heat0p1_latcent75B_pcent070_noqbo"
    "MiMA_heat0p1_latcent75S_pcent030_noqbo"
    "MiMA_heat0p1_latcent75S_pcent050_noqbo"
    "MiMA_heat0p1_latcent75S_pcent070_noqbo"
    "MiMA_heat0p1_latcent90B_pcent030_noqbo"
    "MiMA_heat0p1_latcent90B_pcent050_noqbo"
    "MiMA_heat0p1_latcent90B_pcent070_noqbo"
    "MiMA_heat0p1_latcent90S_pcent030_noqbo"
    "MiMA_heat0p1_latcent90S_pcent050_noqbo"
    "MiMA_heat0p1_latcent90S_pcent070_noqbo"
    "MiMA_heat0p0_noqbo"
)

# Find unfinished experiments
declare -a unfinished_experiments=()

echo "Checking for unfinished experiments..."
for exp_name in "${all_experiments[@]}"; do
    # Convert to lowercase for directory name
    exp_dir_name=$(echo "$exp_name" | tr '[:upper:]' '[:lower:]')
    run1428_dir="$GFDL_DATA/$exp_dir_name/run1428"

    if [ ! -d "$run1428_dir" ]; then
        unfinished_experiments+=("$exp_name")
        echo "  Need to run: $exp_name"
    fi
done

num_unfinished=${#unfinished_experiments[@]}
echo ""
echo "Found $num_unfinished unfinished experiments"

if [ $num_unfinished -eq 0 ]; then
    echo "All experiments are complete!"
    exit 0
fi

# Calculate number of batches needed
num_batches=$(( (num_unfinished + BATCH_SIZE - 1) / BATCH_SIZE ))
echo "Will create $num_batches batch job(s) of up to $BATCH_SIZE experiments each"
echo ""

# Create batch scripts and submit them
for ((batch=0; batch<num_batches; batch++)); do
    batch_num=$((batch + 1))
    batch_script="run_mima_noqbo_batch${batch_num}.sh"

    # Start batch script
    cat > "$batch_script" << 'EOF_HEADER'
#!/bin/bash
#SBATCH --job-name=mima_noqbo
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=32
#SBATCH --time=48:00:00
#SBATCH -p general
#SBATCH -A r00132
#SBATCH --mail-type=ALL
#SBATCH --mail-user=pwstaten@iu.edu
EOF_HEADER

    # Add output/error file names with batch number
    echo "#SBATCH -o mima_noqbo_batch${batch_num}_%j.txt" >> "$batch_script"
    echo "#SBATCH -e mima_noqbo_batch${batch_num}_%j.err" >> "$batch_script"
    echo "" >> "$batch_script"

    # Add environment setup
    cat >> "$batch_script" << 'EOF_ENV'
# Load environment
source ~/.bashrc
conda activate isca_env

# Compile the codebase before launching, but ONLY if the build is missing. This
# self-heals a scratch-purged build (rebuilds isca.x + restores the mppnccombine.x
# symlink). Guarding on presence means jobs skip compiling when a good build already
# exists, so concurrent batch jobs don't race on the shared build dir.
BUILD=$(ls -d $GFDL_WORK/codebase/*/build/isca 2>/dev/null | head -1)
if [ -z "$BUILD" ] || [ ! -e "$BUILD/isca.x" ] || [ ! -e "$BUILD/mppnccombine.x" ]; then
  echo "codebase build missing - compiling..."
  python -c "from isca import IscaCodeBase, GFDL_BASE; cb=IscaCodeBase.from_directory(GFDL_BASE); cb.compile()" || { echo "COMPILE FAILED - aborting"; exit 1; }
else
  echo "codebase build present - skipping compile"
fi

# Launch scripts in parallel
EOF_ENV

    # Add experiments for this batch
    start_idx=$((batch * BATCH_SIZE))
    end_idx=$((start_idx + BATCH_SIZE))
    if [ $end_idx -gt $num_unfinished ]; then
        end_idx=$num_unfinished
    fi

    echo "Batch $batch_num will run experiments $((start_idx + 1)) to $end_idx:"

    for ((i=start_idx; i<end_idx; i++)); do
        exp_name="${unfinished_experiments[$i]}"
        echo "  - $exp_name"
        echo "srun -n1 --cpus-per-task=32 python $EXP_DIR/${exp_name}.py &" >> "$batch_script"
    done

    # Add wait command
    echo "" >> "$batch_script"
    echo "wait" >> "$batch_script"

    # Make executable
    chmod +x "$batch_script"

    # Submit to SLURM
    echo "Submitting $batch_script..."
    sbatch "$batch_script"

    echo ""
done

echo "========================================="
echo "Summary:"
echo "  Total experiments: ${#all_experiments[@]}"
echo "  Unfinished: $num_unfinished"
echo "  Batches submitted: $num_batches"
echo "========================================="
echo ""
echo "Generated batch scripts:"
for ((batch=0; batch<num_batches; batch++)); do
    batch_num=$((batch + 1))
    echo "  run_mima_noqbo_batch${batch_num}.sh"
done
echo ""
echo "You can check job status with: squeue -u $USER"
