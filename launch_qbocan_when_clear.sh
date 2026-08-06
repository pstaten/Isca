#!/bin/bash
#SBATCH --job-name=qbocan_gate
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --time=01:00:00
#SBATCH -p general
#SBATCH -A r00132
#SBATCH -o /dev/null
#SBATCH -e /dev/null

# Gate for the canonical-QBO pair: waits (hourly) until the running heat0p4 and fixh2o
# campaigns are COMPLETELY finished, then recompiles the shared build ONCE (picking up
# do_tab_qbo + the QBO phase int-overflow fix -- a recompile mid-campaign would change
# those campaigns' remaining segments), launches run_qbocan_jobs.sh + its babysitter,
# and stops. Log: ~/isca_qbocan_gate.log

source ~/.bashrc
LOG=$HOME/isca_qbocan_gate.log

if squeue -u pwstaten -h -n mima_heat0p4,mima_fixh2o 2>/dev/null | grep -q .; then
  echo "$(date '+%F %T') campaign jobs still queued - waiting." >> "$LOG"
  sbatch --begin=now+60minutes "$GFDL_BASE/launch_qbocan_when_clear.sh" >/dev/null 2>&1
  exit 0
fi

alldone=1
for e in $(cd "$GFDL_BASE/exp/MiMA" && ls MiMA_*heat0p4*.py | sed 's/\.py$//' | tr '[:upper:]' '[:lower:]') \
         mima_heat0p0_noqbo_fixh2o mima_heat0p1_latcent00b_pcent050_noqbo_fixh2o mima_heat0p1_latcent90b_pcent050_noqbo_fixh2o; do
  last=$(ls -d "$GFDL_DATA/$e"/run???? 2>/dev/null | tail -1)
  last=${last##*/run}; last=${last:-0000}
  [ "$((10#$last))" -ge 1428 ] || { alldone=0; echo "$(date '+%F %T') waiting on $e (at $last)" >> "$LOG"; break; }
done

if [ "$alldone" -eq 0 ]; then
  sbatch --begin=now+60minutes "$GFDL_BASE/launch_qbocan_when_clear.sh" >/dev/null 2>&1
  exit 0
fi

echo "$(date '+%F %T') campaigns complete - recompiling with do_tab_qbo + phase fix." >> "$LOG"
conda activate isca_env
cd "$GFDL_BASE" || exit 1
python -c "from isca import IscaCodeBase, GFDL_BASE; cb=IscaCodeBase.from_directory(GFDL_BASE); cb.compile()" >> "$LOG" 2>&1 \
  || { echo "$(date '+%F %T') COMPILE FAILED - stopping (fix manually, then sbatch run_qbocan_jobs.sh)." >> "$LOG"; exit 1; }
echo "$(date '+%F %T') compile OK - launching canonical-QBO pair + babysitter." >> "$LOG"
sbatch run_qbocan_jobs.sh >> "$LOG" 2>&1
sbatch --begin=now+60minutes auto_resubmit_qbocan.sh >> "$LOG" 2>&1
echo "$(date '+%F %T') gate done." >> "$LOG"
