#!/usr/bin/env bash
set -euo pipefail

PROJECT="/home/zhuweiyu/codex-r/GBC_STK31_public_scRNA"
RAR="${PROJECT}/02_processed_data/All_single_cells/GBC_counts.rar"
DEST="${PROJECT}/02_processed_data/All_single_cells/counts"
UNRAR="${PROJECT}/00_tools/bin/unrar"
LOG="${PROJECT}/logs/stage2_extract_counts.log"

mkdir -p "${DEST}"
"${UNRAR}" x -o+ -idq "${RAR}" "${DEST}/" |& tee "${LOG}"

find "${DEST}" -type f -printf 'EXTRACTED_COUNT_FILE\t%p\t%s bytes\n' | sort
printf 'STAGE2_COUNTS_EXTRACTION_COMPLETE\n'
