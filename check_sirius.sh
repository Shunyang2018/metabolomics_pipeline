#!/bin/bash
# Check SIRIUS progress every 2 minutes
while true; do
    if ps aux | grep "java.*sirius.*-i.*sirius_unknown" | grep -v grep > /dev/null; then
        echo "[$(date '+%H:%M:%S')] SIRIUS still running..."
        ls -lh /Users/wangs261/Documents/project/excel_merge/test_csv/outputs/sirius*.sirius 2>/dev/null | tail -2
    else
        echo "[$(date '+%H:%M:%S')] SIRIUS COMPLETED!"
        ls -lh /Users/wangs261/Documents/project/excel_merge/test_csv/outputs/sirius*.sirius 2>/dev/null
        echo ""
        echo "Now run: metabo final"
        break
    fi
    sleep 120
done
