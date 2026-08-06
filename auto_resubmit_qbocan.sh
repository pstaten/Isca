#!/bin/bash
#SBATCH --job-name=auto_resub_qbocan
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=00:10:00
#SBATCH -p general
#SBATCH -A r00132
#SBATCH -o /dev/null
#SBATCH -e /dev/null

# Hourly babysitter for the two canonical-QBO fixed-SST runs (same pattern as
# auto_resubmit_controller.sh): while any experiment is below run1428 and no
# runner is queued, resubmit run_qbocan_jobs.sh; stop when done, or when nothing
# advanced since the last resubmit (anti-storm guard for a permanently stuck run).
# Log: ~/isca_autoresub_qbocan.log  State: ~/.isca_autoresub_qbocan_state

source ~/.bashrc
EXPS="mima_heat0p0_qbocan_fixsst mima_heat0p1_qbocan_fixsst"
LOG=$HOME/isca_autoresub_qbocan.log
STATE=$HOME/.isca_autoresub_qbocan_state

progress=""; done_all=1
for e in $EXPS; do
  last=$(ls -d "$GFDL_DATA/$e"/run???? 2>/dev/null | tail -1)
  last=${last##*/run}; last=${last:-0000}
  progress="$progress$e:$((10#$last)) "
  [ "$((10#$last))" -ge 1428 ] || done_all=0
done

if [ "$done_all" -eq 1 ]; then
  echo "$(date '+%F %T') ALL qbocan experiments reached run1428 - complete. Stopping." >> "$LOG"
  exit 0
fi

if squeue -u pwstaten -h -n mima_qbocan 2>/dev/null | grep -q .; then
  echo "$(date '+%F %T') qbocan runner still in queue - waiting. [$progress]" >> "$LOG"
else
  if [ -f "$STATE" ] && [ "$(cat "$STATE")" = "$progress" ]; then
    echo "$(date '+%F %T') NO PROGRESS since last resubmit [$progress] - stopping (anti-storm)." >> "$LOG"
    exit 0
  fi
  echo "$progress" > "$STATE"
  cd "$GFDL_BASE" && sbatch run_qbocan_jobs.sh >> "$LOG" 2>&1
  echo "$(date '+%F %T') resubmitted runner. [$progress]" >> "$LOG"
fi

sbatch --begin=now+60minutes "$GFDL_BASE/auto_resubmit_qbocan.sh" >> "$LOG" 2>&1
echo "$(date '+%F %T') rescheduled +1h" >> "$LOG"
