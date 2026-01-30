#!/bin/bash
#
# Delete the last run directory for each MiMA experiment.
# The last run may have incomplete output due to job termination.
#

if [ -z "$GFDL_DATA" ]; then
    echo "ERROR: GFDL_DATA environment variable is not set"
    exit 1
fi

echo "GFDL_DATA: $GFDL_DATA"
echo ""

deleted=0
skipped=0

for exp_dir in "$GFDL_DATA"/mima_*/; do
    # Skip if no matching directories
    [ -d "$exp_dir" ] || continue

    exp_name=$(basename "$exp_dir")

    # Find the highest-numbered run directory
    last_run=$(ls -d "$exp_dir"run[0-9]* 2>/dev/null | sort -t 'n' -k1 -V | tail -1)

    if [ -z "$last_run" ]; then
        echo "SKIP: $exp_name (no run directories found)"
        ((skipped++))
        continue
    fi

    last_run_name=$(basename "$last_run")
    echo "DELETE: $exp_name/$last_run_name"
    rm -rf "$last_run"
    ((deleted++))
done

echo ""
echo "========================================="
echo "Summary:"
echo "  Deleted: $deleted"
echo "  Skipped: $skipped"
echo "========================================="
