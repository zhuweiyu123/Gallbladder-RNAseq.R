#!/usr/bin/env bash
set -u

PROJECT="/home/zhuweiyu/codex-r/GBC_STK31_public_scRNA"
RAW="${PROJECT}/01_raw_data"
ARIA2C="${PROJECT}/04_scripts/run_aria2c.sh"
URL="https://zenodo.org/api/records/15400138/files/NK.RDS/content"
EXPECTED="2282d943cac04844c9580275945672af"
FILE="NK.RDS"

mkdir -p "${RAW}"
cd "${RAW}"
attempt=0
while true; do
  attempt=$((attempt + 1))
  printf 'NK_DOWNLOAD_ATTEMPT\t%s\t%s\n' "${attempt}" "$(date --iso-8601=seconds)"
  "${ARIA2C}" --continue=true --split=16 --max-connection-per-server=16 \
    --max-tries=0 --retry-wait=30 --timeout=60 --min-split-size=16M \
    --file-allocation=none --auto-file-renaming=false --summary-interval=10 \
    --console-log-level=notice --out="${FILE}" "${URL}"
  rc=$?
  if [[ ${rc} -eq 0 && -f "${FILE}" ]]; then
    actual="$(md5sum "${FILE}" | awk '{print $1}')"
    if [[ "${actual}" == "${EXPECTED}" ]]; then
      printf 'NK_MD5_OK\t%s\n' "${actual}"
      exit 0
    fi
    printf 'NK_MD5_MISMATCH\texpected=%s\tactual=%s\n' "${EXPECTED}" "${actual}" >&2
    exit 2
  fi
  printf 'NK_DOWNLOAD_RETRY\trc=%s\twait=60s\n' "${rc}" >&2
  sleep 60
done
