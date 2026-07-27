#!/bin/bash
# Launch (or resume) all fixed-SST MiMA experiments, batched 4 per 128-core SLURM job.
# Idempotent: exp.run skips already-complete run dirs, so re-running this resumes each experiment
# from its last restart. Each batch job self-heals poison run dirs (run{N} lacking res{N}) and
# compiles only if the build is missing. Run on BigRed.
set -u
ISCA=/N/slate/pwstaten/Projects/Isca
cd "$ISCA"
mapfile -t SCRIPTS < <(ls exp/MiMA/MiMA_*_fixsst.py 2>/dev/null | grep -v 'MiMA_test_fixsst.py' | sort)
n=${#SCRIPTS[@]}
echo "launching $n fixed-SST experiments in batches of 4"
b=0
for ((i=0; i<n; i+=4)); do
  b=$((b+1))
  chunk=("${SCRIPTS[@]:i:4}")
  nt=${#chunk[@]}
  jf=~/fixsst_batch${b}.slurm
  {
    echo "#!/bin/bash"
    echo "#SBATCH --job-name=fixsst_b${b}"
    echo "#SBATCH --nodes=1"
    echo "#SBATCH --ntasks=${nt}"
    echo "#SBATCH --cpus-per-task=32"
    echo "#SBATCH --time=48:00:00"
    echo "#SBATCH -p general"
    echo "#SBATCH -A r00132"
    echo "#SBATCH -o fixsst_b${b}_%j.txt"
    echo "#SBATCH -e fixsst_b${b}_%j.err"
    echo "#SBATCH --mail-type=FAIL"
    echo "#SBATCH --mail-user=pwstaten@iu.edu"
    echo 'source ~/.bashrc; conda activate isca_env'
    echo 'BUILD=$(ls -d $GFDL_WORK/codebase/*/build/isca 2>/dev/null | head -1)'
    echo 'if [ -z "$BUILD" ] || [ ! -e "$BUILD/isca.x" ] || [ ! -e "$BUILD/mppnccombine.x" ]; then echo COMPILING; python -c "from isca import IscaCodeBase, GFDL_BASE; cb=IscaCodeBase.from_directory(GFDL_BASE); cb.compile()" || { echo COMPILE FAILED; exit 1; }; fi'
    for s in "${chunk[@]}"; do
      e=$(grep -oE "Experiment\(.[a-z0-9_]+_fixsst" "$ISCA/$s" | sed "s/Experiment(.//")
      echo "for rd in \"\$GFDL_DATA/$e\"/run[0-9]*; do [ -d \"\$rd\" ] || continue; nn=\$(basename \"\$rd\"); nn=\${nn#run}; [ -f \"\$GFDL_DATA/$e/restarts/res\${nn}.tar.gz\" ] || rm -rf \"\$rd\"; done"
    done
    for s in "${chunk[@]}"; do
      echo "srun -n1 --cpus-per-task=32 python $ISCA/$s &"
    done
    echo "wait"
  } > "$jf"
  jid=$(sbatch "$jf" | awk '{print $NF}')
  echo "  batch $b ($nt exps): job $jid"
done
echo "done."
