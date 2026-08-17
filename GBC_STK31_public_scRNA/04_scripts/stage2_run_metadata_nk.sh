#!/usr/bin/env bash
set -euo pipefail

PROJECT="/home/zhuweiyu/codex-r/GBC_STK31_public_scRNA"
cd "${PROJECT}"

/usr/bin/time -v Rscript 04_scripts/stage2_inspect_metadata_nk.R \
  > logs/stage2_inspect_metadata_nk.stdout.log \
  2> logs/stage2_inspect_metadata_nk.stderr.log

tail -40 logs/stage2_inspect_metadata_nk.stdout.log
tail -40 logs/stage2_inspect_metadata_nk.stderr.log
