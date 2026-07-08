#!/bin/bash
#SBATCH --job-name=auto_resubmit
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=00:10:00
#SBATCH -p general
#SBATCH -A r00132
#SBATCH -o auto_resubmit_%j.txt
#SBATCH -e auto_resubmit_%j.err
#
# Self-rescheduling babysitter for the MiMA runs. Wakes ~hourly; when no run-jobs are
# in the queue, resubmits whatever hasn't reached run TARGET via the existing launchers,
# then reschedules itself. Stops when everything is done OR nothing advanced since the
# last resubmit (guards against an endless loop on a permanently-stuck experiment).
#
# Kick off ONCE:  sbatch auto_resubmit_controller.sh
# Stop early:     scancel the pending auto_resubmit job (and it won't reschedule).
#
# Test without side effects:
#   DRY_RUN=1 bash auto_resubmit_controller.sh                  # exercise the wait path
#   DRY_RUN=1 TEST_FORCE_EMPTY=1 bash auto_resubmit_controller.sh   # exercise resubmit decision

source ~/.bashrc 2>/dev/null

SELF="$GFDL_BASE/auto_resubmit_controller.sh"
TARGET=1428
INTERVAL_H=1
RUNJOBS_RE="mima_noqbo|mima_runner|mima_hiheat_runner"
STATE="$HOME/.isca_autoresubmit_state"
LOG="$HOME/isca_autoresubmit.log"

QBO_EXPS="mima_heat0p0_qbo00 mima_heat0p0_qbo20 mima_heat0p1_qbo00 mima_heat0p1_qbo20"
HIHEAT_EXPS="mima_hiheat0p1_qbo00 mima_hiheat0p1_qbo20"
SWEEP_EXPS=$(cd "$GFDL_DATA" 2>/dev/null && ls -d mima_*_noqbo 2>/dev/null | tr '\n' ' ')

log() { echo "$(date '+%F %T') $*" | tee -a "$LOG"; }
max_run() { ls -d "$GFDL_DATA/$1"/run[0-9]* 2>/dev/null | sed 's:.*/run::' | sort -n | tail -1; }
any_unfinished() { local e c; for e in $1; do c=$(max_run "$e"); [ "$((10#${c:-0}))" -lt "$TARGET" ] && return 0; done; return 1; }

do_sbatch() { if [ -n "$DRY_RUN" ]; then log "[dry-run] would: $*"; else "$@" >>"$LOG" 2>&1; fi; }
do_submit_sweep() { if [ -n "$DRY_RUN" ]; then log "[dry-run] would: bash submit_mima_noqbo_batches.sh"; else bash "$GFDL_BASE/submit_mima_noqbo_batches.sh" >>"$LOG" 2>&1; fi; }
reschedule() { if [ -n "$DRY_RUN" ]; then log "[dry-run] would reschedule +${INTERVAL_H}h"; else sbatch --begin=now+${INTERVAL_H}hour "$SELF" >>"$LOG" 2>&1 && log "rescheduled +${INTERVAL_H}h"; fi; }

# 1) Run-jobs still active? Wait and re-check later.
if [ -z "$TEST_FORCE_EMPTY" ] && squeue -u "$USER" -h -o '%j' 2>/dev/null | grep -qE "$RUNJOBS_RE"; then
  log "run-jobs still in queue — waiting."
  reschedule; exit 0
fi

# 2) Queue empty. Measure current progress and compare to the state at last resubmit.
progressed=0; anyunfinished=0
newstate="$(mktemp)"
for e in $SWEEP_EXPS $QBO_EXPS $HIHEAT_EXPS; do
  raw=$(max_run "$e" 2>/dev/null); c=$((10#${raw:-0}))
  echo "$e $c" >> "$newstate"
  [ "$c" -lt "$TARGET" ] && anyunfinished=1
  prev=$(awk -v x="$e" '$1==x{print $2}' "$STATE" 2>/dev/null); prev=${prev:--1}
  [ "$c" -gt "$prev" ] && progressed=1
done

if [ "$anyunfinished" = 0 ]; then
  log "ALL experiments reached run$TARGET — auto-resubmit complete. Stopping (no reschedule)."
  [ -z "$DRY_RUN" ] && mv "$newstate" "$STATE" || rm -f "$newstate"; exit 0
fi
if [ "$progressed" = 0 ]; then
  log "Nothing advanced since last resubmit — likely stuck. NOT resubmitting; stopping. Investigate."
  [ -z "$DRY_RUN" ] && mv "$newstate" "$STATE" || rm -f "$newstate"; exit 0
fi

# 3) Unfinished and progressing -> resubmit the unfinished launchers.
log "queue empty, unfinished + progressing -> resubmitting"
any_unfinished "$SWEEP_EXPS"  && do_submit_sweep                                    && log "-> noqbo sweep"
any_unfinished "$QBO_EXPS"    && do_sbatch sbatch "$GFDL_BASE/run_mima_jobs.sh"     && log "-> qbo runner"
any_unfinished "$HIHEAT_EXPS" && do_sbatch sbatch "$GFDL_BASE/run_mima_hiheat_jobs.sh" && log "-> hiheat runner"
[ -z "$DRY_RUN" ] && mv "$newstate" "$STATE" || rm -f "$newstate"
reschedule
