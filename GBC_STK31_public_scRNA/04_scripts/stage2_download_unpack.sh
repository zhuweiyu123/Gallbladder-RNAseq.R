#!/usr/bin/env bash
set -euo pipefail

PROJECT="/home/zhuweiyu/codex-r/GBC_STK31_public_scRNA"
RAW="${PROJECT}/01_raw_data"
PROCESSED="${PROJECT}/02_processed_data"
FEASIBILITY="${PROJECT}/06_feasibility"
UNRAR="${PROJECT}/00_tools/bin/unrar"
ARIA2C="${PROJECT}/04_scripts/run_aria2c.sh"

ALL_FILE="scRNA_processed_data_all.zip"
NK_FILE="NK.RDS"
ALL_MD5="922631b810f83653b85bd479084b7a77"
NK_MD5="2282d943cac04844c9580275945672af"
ALL_URL="https://zenodo.org/api/records/15400138/files/scRNA_processed_data_all.zip/content"
NK_URL="https://zenodo.org/api/records/15400138/files/NK.RDS/content"

mkdir -p "${RAW}" "${PROCESSED}" "${FEASIBILITY}"

download_file() {
  local url="$1"
  local output="$2"
  local expected_md5="$3"
  if [[ -f "${output}" && ! -f "${output}.aria2" ]]; then
    local existing_md5
    existing_md5="$(md5sum "${output}" | awk '{print $1}')"
    if [[ "${existing_md5}" == "${expected_md5}" ]]; then
      printf 'DOWNLOAD_ALREADY_COMPLETE\t%s\t%s\n' "${output}" "${existing_md5}"
      return 0
    fi
  fi
  if [[ -x "${ARIA2C}" ]]; then
    "${ARIA2C}" \
      --continue=true \
      --split=16 \
      --max-connection-per-server=16 \
      --max-tries=0 \
      --retry-wait=30 \
      --timeout=60 \
      --min-split-size=16M \
      --file-allocation=none \
      --auto-file-renaming=false \
      --summary-interval=10 \
      --console-log-level=notice \
      --out="${output}" \
      "${url}"
  else
    curl -fL \
      --retry 10 \
      --retry-all-errors \
      --retry-delay 30 \
      --connect-timeout 30 \
      --continue-at - \
      --output "${output}" \
      "${url}"
  fi
}

verify_md5() {
  local expected="$1"
  local file="$2"
  local actual
  actual="$(md5sum "${file}" | awk '{print $1}')"
  if [[ "${actual}" != "${expected}" ]]; then
    printf 'MD5_MISMATCH\t%s\texpected=%s\tactual=%s\n' "${file}" "${expected}" "${actual}" >&2
    exit 2
  fi
  printf 'MD5_OK\t%s\t%s\n' "${file}" "${actual}"
}

cd "${RAW}"
download_file "${ALL_URL}" "${ALL_FILE}" "${ALL_MD5}"
verify_md5 "${ALL_MD5}" "${ALL_FILE}"

download_file "${NK_URL}" "${NK_FILE}" "${NK_MD5}"
verify_md5 "${NK_MD5}" "${NK_FILE}"

unzip -o "${ALL_FILE}" -d "${PROCESSED}"

COUNTS_RAR="${PROCESSED}/All_single_cells/GBC_counts.rar"
METADATA="${PROCESSED}/All_single_cells/GBC_Metadata.txt"
SAMPLE_INFO="${PROCESSED}/All_single_cells/sample_info.xlsx"

for required in "${COUNTS_RAR}" "${METADATA}" "${SAMPLE_INFO}"; do
  if [[ ! -f "${required}" ]]; then
    printf 'REQUIRED_FILE_MISSING\t%s\n' "${required}" >&2
    exit 3
  fi
  stat --printf='EXTRACTED_FILE\t%n\t%s bytes\n' "${required}"
done

"${UNRAR}" l "${COUNTS_RAR}" | tee "${FEASIBILITY}/counts_structure_unrar_l.txt"

printf 'STAGE2_DOWNLOAD_UNPACK_COMPLETE\n'
