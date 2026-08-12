#!/bin/bash
#SBATCH --job-name=auto_resub_symsst
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=00:10:00
#SBATCH -p general
#SBATCH -A r00132
#SBATCH -o /dev/null
#SBATCH -e /dev/null

# Hourly babysitter for the 38-experiment 4x-heating (0.4 K/day) slab campaign
# (same pattern as auto_resubmit_controller.sh / auto_resubmit_fixh2o.sh): while any
# experiment is below run0600 and no mima_symsst batch job is queued, rerun the
# self-detecting batch submitter; stop when all done, or when nothing advanced since
# the last resubmit (anti-storm guard).
# Log: ~/isca_autoresub_symsst.log  State: ~/.isca_autoresub_symsst_state

source ~/.bashrc
EXPS="mima_heat0p0_noqbo_symsst mima_heat0p1_latcent00b_pcent030_noqbo_symsst mima_heat0p1_latcent00b_pcent050_noqbo_symsst mima_heat0p1_latcent30b_pcent050_noqbo_symsst mima_heat0p1_latcent45s_pcent050_noqbo_symsst mima_heat0p1_latcent45s_pcent125_noqbo_symsst mima_heat0p1_latcent60b_pcent050_noqbo_symsst mima_heat0p1_latcent60b_pcent070_noqbo_symsst mima_heat0p1_latcent60s_pcent125_noqbo_symsst mima_heat0p1_latcent75s_pcent125_noqbo_symsst mima_heat0p1_latcent90b_pcent050_noqbo_symsst mima_heat0p1_latcent90s_pcent050_noqbo_symsst"
LOG=$HOME/isca_autoresub_symsst.log
STATE=$HOME/.isca_autoresub_symsst_state

progress=""; done_all=1
for e in $EXPS; do
  last=$(ls -d "$GFDL_DATA/$e"/run???? 2>/dev/null | tail -1)
  last=${last##*/run}; last=${last:-0000}
  progress="$progress$((10#$last)),"
  [ "$((10#$last))" -ge 600 ] || done_all=0
done

if [ "$done_all" -eq 1 ]; then
  echo "$(date '+%F %T') ALL symsst experiments reached run0600 - complete. Stopping." >> "$LOG"
  exit 0
fi

if squeue -u pwstaten -h -n mima_symsst 2>/dev/null | grep -q .; then
  echo "$(date '+%F %T') symsst batch jobs still in queue - waiting." >> "$LOG"
else
  if [ -f "$STATE" ] && [ "$(cat "$STATE")" = "$progress" ]; then
    echo "$(date '+%F %T') NO PROGRESS since last resubmit - stopping (anti-storm). [$progress]" >> "$LOG"
    exit 0
  fi
  echo "$progress" > "$STATE"
  cd "$GFDL_BASE" && bash submit_mima_symsst_batches.sh >> "$LOG" 2>&1
  echo "$(date '+%F %T') ran batch submitter. [$progress]" >> "$LOG"
fi

sbatch --begin=now+60minutes "$GFDL_BASE/auto_resubmit_symsst.sh" >/dev/null 2>&1
echo "$(date '+%F %T') rescheduled +1h" >> "$LOG"
