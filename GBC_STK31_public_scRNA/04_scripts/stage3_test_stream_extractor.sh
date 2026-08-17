#!/usr/bin/env bash
set -euo pipefail

PROJECT="/home/zhuweiyu/codex-r/GBC_STK31_public_scRNA"
WORK="${PROJECT}/logs/stage3_stream_unit_test"
SOURCE="${PROJECT}/04_scripts/stage3_stream_extract_malignant.cpp"
BINARY="${WORK}/extractor"

mkdir -p "${WORK}"
g++ -O3 -std=c++17 -Wall -Wextra -Wpedantic "${SOURCE}" -lz -o "${BINARY}"

printf '%%%%MatrixMarket matrix coordinate integer general\n4 5 6\n1 1 2\n4 1 3\n2 2 4\n3 3 5\n1 4 6\n2 5 7\n' | gzip -c > "${WORK}/input.mtx.gz"

Rscript -e '
  writeBin(as.integer(c(0, 1, 0, 2, 3)), "'"${WORK}"'/col.bin", size = 4L, endian = "little")
  writeBin(as.integer(c(1, 0, 2, 0)), "'"${WORK}"'/row.bin", size = 4L, endian = "little")
'

"${BINARY}" \
  "${WORK}/input.mtx.gz" "${WORK}/col.bin" "${WORK}/row.bin" \
  4 5 3 2 \
  "${WORK}/full.mtx.gz" "${WORK}/target.mtx.gz" \
  "${WORK}/library.tsv" "${WORK}/audit.txt" "${WORK}/tmp"

gzip -cd "${WORK}/full.mtx.gz" > "${WORK}/full.actual"
gzip -cd "${WORK}/target.mtx.gz" > "${WORK}/target.actual"

printf '%%%%MatrixMarket matrix coordinate integer general\n4 3 4\n2 1 4\n1 2 6\n2 3 7\n' > "${WORK}/full.expected"
# Full subset also contains row 1/col 2 only and row 2/col 1+3; expected nnz is 3.
printf '%%%%MatrixMarket matrix coordinate integer general\n4 3 3\n2 1 4\n1 2 6\n2 3 7\n' > "${WORK}/full.expected"
printf '%%%%MatrixMarket matrix coordinate integer general\n2 3 1\n1 2 6\n' > "${WORK}/target.expected"
printf 'malignant_col\tall_gene_raw_umi\n1\t4\n2\t6\n3\t7\n' > "${WORK}/library.expected"

cmp "${WORK}/full.expected" "${WORK}/full.actual"
cmp "${WORK}/target.expected" "${WORK}/target.actual"
cmp "${WORK}/library.expected" "${WORK}/library.tsv"
grep -q '^actual_nnz=6$' "${WORK}/audit.txt"
grep -q '^malignant_output_nnz=3$' "${WORK}/audit.txt"
grep -q '^target_output_nnz=1$' "${WORK}/audit.txt"

printf 'STAGE3_STREAM_UNIT_TEST_OK\n'
