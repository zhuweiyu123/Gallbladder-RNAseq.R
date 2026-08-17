#!/usr/bin/env bash
set -euo pipefail

TOOL_ROOT="/home/zhuweiyu/codex-r/GBC_STK31_public_scRNA/00_tools/aria2"
export LD_LIBRARY_PATH="${TOOL_ROOT}/usr/lib/x86_64-linux-gnu${LD_LIBRARY_PATH:+:${LD_LIBRARY_PATH}}"
exec "${TOOL_ROOT}/usr/bin/aria2c" "$@"
