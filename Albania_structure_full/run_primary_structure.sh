#!/usr/bin/env bash

set -u

K="$1"
REP="$2"

SEED=$((20260805 + K * 1000 + REP))

PREFIX=$(printf "Albania_primary_K%02d_rep%02d" "$K" "$REP")

OUTPUT="primary_output/${PREFIX}"
LOG="primary_logs/${PREFIX}.console.log"

echo "Starting K=${K}, replicate=${REP}, seed=${SEED}"

./structure \
  -m mainparams_primary \
  -e extraparams_primary \
  -i Albania_2000_STRUCTURE.txt \
  -o "$OUTPUT" \
  -K "$K" \
  -N 2000 \
  -L 16 \
  -D "$SEED" \
  > "$LOG" 2>&1

STATUS=$?

if [ "$STATUS" -eq 0 ] && [ -s "${OUTPUT}_f" ]; then
    echo "Completed K=${K}, replicate=${REP}"
else
    echo "FAILED K=${K}, replicate=${REP}, status=${STATUS}"
fi

exit "$STATUS"
