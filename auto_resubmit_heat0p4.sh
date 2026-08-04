#!/bin/bash
#SBATCH --job-name=auto_resub_heat0p4
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
# experiment is below run1428 and no mima_heat0p4 batch job is queued, rerun the
# self-detecting batch submitter; stop when all done, or when nothing advanced since
# the last resubmit (anti-storm guard).
# Log: ~/isca_autoresub_heat0p4.log  State: ~/.isca_autoresub_heat0p4_state

source ~/.bashrc
EXPS="mima_heat0p4_latcent00b_pcent030_noqbo mima_heat0p4_latcent00b_pcent050_noqbo mima_heat0p4_latcent15b_pcent030_noqbo mima_heat0p4_latcent15b_pcent050_noqbo mima_heat0p4_latcent15s_pcent030_noqbo mima_heat0p4_latcent15s_pcent050_noqbo mima_heat0p4_latcent30b_pcent030_noqbo mima_heat0p4_latcent30b_pcent050_noqbo mima_heat0p4_latcent30s_pcent030_noqbo mima_heat0p4_latcent30s_pcent050_noqbo mima_heat0p4_latcent45b_pcent030_noqbo mima_heat0p4_latcent45b_pcent050_noqbo mima_heat0p4_latcent45b_pcent070_noqbo mima_heat0p4_latcent45s_pcent030_noqbo mima_heat0p4_latcent45s_pcent050_noqbo mima_heat0p4_latcent45s_pcent070_noqbo mima_heat0p4_latcent60b_pcent030_noqbo mima_heat0p4_latcent60b_pcent050_noqbo mima_heat0p4_latcent60b_pcent070_noqbo mima_heat0p4_latcent60s_pcent030_noqbo mima_heat0p4_latcent60s_pcent050_noqbo mima_heat0p4_latcent60s_pcent070_noqbo mima_heat0p4_latcent75b_pcent030_noqbo mima_heat0p4_latcent75b_pcent050_noqbo mima_heat0p4_latcent75b_pcent070_noqbo mima_heat0p4_latcent75s_pcent030_noqbo mima_heat0p4_latcent75s_pcent050_noqbo mima_heat0p4_latcent75s_pcent070_noqbo mima_heat0p4_latcent90b_pcent030_noqbo mima_heat0p4_latcent90b_pcent050_noqbo mima_heat0p4_latcent90b_pcent070_noqbo mima_heat0p4_latcent90s_pcent030_noqbo mima_heat0p4_latcent90s_pcent050_noqbo mima_heat0p4_latcent90s_pcent070_noqbo mima_heat0p4_qbo00 mima_heat0p4_qbo20 mima_hiheat0p4_qbo00 mima_hiheat0p4_qbo20"
LOG=$HOME/isca_autoresub_heat0p4.log
STATE=$HOME/.isca_autoresub_heat0p4_state

progress=""; done_all=1
for e in $EXPS; do
  last=$(ls -d "$GFDL_DATA/$e"/run???? 2>/dev/null | tail -1)
  last=${last##*/run}; last=${last:-0000}
  progress="$progress$((10#$last)),"
  [ "$((10#$last))" -ge 1428 ] || done_all=0
done

if [ "$done_all" -eq 1 ]; then
  echo "$(date '+%F %T') ALL heat0p4 experiments reached run1428 - complete. Stopping." >> "$LOG"
  exit 0
fi

if squeue -u pwstaten -h -n mima_heat0p4 2>/dev/null | grep -q .; then
  echo "$(date '+%F %T') heat0p4 batch jobs still in queue - waiting." >> "$LOG"
else
  if [ -f "$STATE" ] && [ "$(cat "$STATE")" = "$progress" ]; then
    echo "$(date '+%F %T') NO PROGRESS since last resubmit - stopping (anti-storm). [$progress]" >> "$LOG"
    exit 0
  fi
  echo "$progress" > "$STATE"
  cd "$GFDL_BASE" && bash submit_mima_heat0p4_batches.sh >> "$LOG" 2>&1
  echo "$(date '+%F %T') ran batch submitter. [$progress]" >> "$LOG"
fi

sbatch --begin=now+60minutes "$GFDL_BASE/auto_resubmit_heat0p4.sh" >/dev/null 2>&1
echo "$(date '+%F %T') rescheduled +1h" >> "$LOG"
