#!/usr/bin/env bash
set -euo pipefail

PROJECT="/home/zhuweiyu/codex-r/GBC_STK31_public_scRNA"
COMPONENT="${1:-}"

case "$COMPONENT" in
  malignant)
    SCRIPT="$PROJECT/04_scripts/stage3_malignant_detection.R"
    LOG="$PROJECT/logs/stage3_malignant_detection.log"
    TIME_FILE="$PROJECT/06_feasibility/stage3/stage3_malignant_detection.time.txt"
    ;;
  nk)
    SCRIPT="$PROJECT/04_scripts/stage3_nk_detection.R"
    LOG="$PROJECT/logs/stage3_nk_detection.log"
    TIME_FILE="$PROJECT/06_feasibility/stage3/stage3_nk_detection.time.txt"
    ;;
  primary)
    SCRIPT="$PROJECT/04_scripts/stage3_primary_matching.R"
    LOG="$PROJECT/logs/stage3_primary_matching.log"
    TIME_FILE="$PROJECT/06_feasibility/stage3/stage3_primary_matching.time.txt"
    ;;
  *)
    echo "Usage: $0 {malignant|nk|primary}" >&2
    exit 64
    ;;
esac

mkdir -p "$PROJECT/logs" "$PROJECT/06_feasibility/stage3"
/usr/bin/time -v -o "$TIME_FILE" Rscript "$SCRIPT" >"$LOG" 2>&1
echo "STAGE3_${COMPONENT^^}_DETECTION_OK"
