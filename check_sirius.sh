#!/bin/bash
# Watch a running SIRIUS job and report when it finishes.
#
# The output directory defaults to $METABO_INPUT_DIR/outputs, matching
# constants.py, and can be given as the first argument:
#
#   ./check_sirius.sh                    # $METABO_INPUT_DIR/outputs, or data/outputs
#   ./check_sirius.sh path/to/outputs
#   INTERVAL=30 ./check_sirius.sh        # poll every 30s instead of 120s
set -uo pipefail

OUT_DIR="${1:-${METABO_INPUT_DIR:-data}/outputs}"
INTERVAL="${INTERVAL:-120}"

if [ ! -d "$OUT_DIR" ]; then
    echo "Output directory not found: $OUT_DIR" >&2
    echo "Pass it as the first argument or set METABO_INPUT_DIR." >&2
    exit 2
fi

echo "Watching $OUT_DIR every ${INTERVAL}s..."
while true; do
    if ps aux | grep "java.*sirius.*-i.*sirius_unknown" | grep -v grep > /dev/null; then
        echo "[$(date '+%H:%M:%S')] SIRIUS still running..."
        ls -lh "$OUT_DIR"/sirius*.sirius 2>/dev/null | tail -2
    else
        echo "[$(date '+%H:%M:%S')] SIRIUS COMPLETED!"
        ls -lh "$OUT_DIR"/sirius*.sirius 2>/dev/null
        echo ""
        echo "Now run: metabo final"
        break
    fi
    sleep "$INTERVAL"
done
