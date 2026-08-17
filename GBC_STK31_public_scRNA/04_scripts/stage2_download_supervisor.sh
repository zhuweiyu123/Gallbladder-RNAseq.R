#!/usr/bin/env bash
set -u

PROJECT="/home/zhuweiyu/codex-r/GBC_STK31_public_scRNA"
MAIN="${PROJECT}/04_scripts/stage2_download_unpack.sh"

attempt=0
while true; do
  attempt=$((attempt + 1))
  printf 'SUPERVISOR_ATTEMPT\t%s\t%s\n' "${attempt}" "$(date --iso-8601=seconds)"
  "${MAIN}"
  rc=$?
  if [[ ${rc} -eq 0 ]]; then
    printf 'SUPERVISOR_COMPLETE\t%s\n' "$(date --iso-8601=seconds)"
    exit 0
  fi
  if [[ ${rc} -eq 2 || ${rc} -eq 3 ]]; then
    printf 'SUPERVISOR_HARD_STOP\trc=%s\t%s\n' "${rc}" "$(date --iso-8601=seconds)" >&2
    exit "${rc}"
  fi
  printf 'SUPERVISOR_RETRY\trc=%s\twait=60s\n' "${rc}" >&2
  sleep 60
done
