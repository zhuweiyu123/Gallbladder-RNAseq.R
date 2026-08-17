#include <zlib.h>

#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

int main(int argc, char* argv[]) {
  if (argc != 8) {
    std::cerr << "Usage: " << argv[0]
              << " <matrix.mtx.gz> <cell_to_sample.int32.bin> <n_rows> <n_cols>"
              << " <n_samples> <output.mtx.gz> <audit.tsv>\n";
    return 64;
  }

  const std::string matrix_path = argv[1];
  const std::string map_path = argv[2];
  const uint64_t expected_rows = std::strtoull(argv[3], nullptr, 10);
  const uint64_t expected_cols = std::strtoull(argv[4], nullptr, 10);
  const int32_t n_samples = static_cast<int32_t>(std::strtol(argv[5], nullptr, 10));
  const std::string output_path = argv[6];
  const std::string audit_path = argv[7];

  std::ifstream map_in(map_path, std::ios::binary);
  if (!map_in) {
    std::cerr << "Cannot open map: " << map_path << "\n";
    return 65;
  }
  std::vector<int32_t> cell_to_sample(expected_cols, 0);
  map_in.read(reinterpret_cast<char*>(cell_to_sample.data()),
              static_cast<std::streamsize>(expected_cols * sizeof(int32_t)));
  if (map_in.gcount() != static_cast<std::streamsize>(expected_cols * sizeof(int32_t))) {
    std::cerr << "Cell-to-sample map has unexpected byte length\n";
    return 66;
  }

  gzFile input = gzopen(matrix_path.c_str(), "rb");
  if (!input) {
    std::cerr << "Cannot open matrix: " << matrix_path << "\n";
    return 67;
  }
  constexpr int buffer_size = 65536;
  std::vector<char> buffer(buffer_size);
  if (!gzgets(input, buffer.data(), buffer_size)) {
    std::cerr << "Cannot read Matrix Market banner\n";
    gzclose(input);
    return 68;
  }
  const std::string banner(buffer.data());
  if (banner.rfind("%%MatrixMarket matrix coordinate integer general", 0) != 0) {
    std::cerr << "Unexpected Matrix Market banner: " << banner;
    gzclose(input);
    return 69;
  }
  std::string dimension_line;
  while (gzgets(input, buffer.data(), buffer_size)) {
    std::string line(buffer.data());
    if (!line.empty() && line[0] != '%') {
      dimension_line = line;
      break;
    }
  }
  uint64_t rows = 0, cols = 0, declared_nnz = 0;
  {
    std::istringstream parser(dimension_line);
    if (!(parser >> rows >> cols >> declared_nnz) || rows != expected_rows || cols != expected_cols) {
      std::cerr << "Unexpected matrix dimensions: " << dimension_line;
      gzclose(input);
      return 70;
    }
  }

  std::vector<uint64_t> sums(static_cast<size_t>(rows) * static_cast<size_t>(n_samples), 0);
  uint64_t actual_nnz = 0;
  while (gzgets(input, buffer.data(), buffer_size)) {
    uint64_t row = 0, col = 0, value = 0;
    std::istringstream parser(buffer.data());
    if (!(parser >> row >> col >> value)) {
      std::cerr << "Malformed coordinate at nonzero " << (actual_nnz + 1U) << "\n";
      gzclose(input);
      return 71;
    }
    ++actual_nnz;
    if (row == 0 || row > rows || col == 0 || col > cols) {
      std::cerr << "Out-of-range coordinate at nonzero " << actual_nnz << "\n";
      gzclose(input);
      return 72;
    }
    const int32_t sample = cell_to_sample[static_cast<size_t>(col - 1U)];
    if (sample < 0 || sample > n_samples) {
      std::cerr << "Invalid sample index in map\n";
      gzclose(input);
      return 73;
    }
    if (sample > 0) {
      sums[(static_cast<size_t>(row - 1U) * static_cast<size_t>(n_samples)) + static_cast<size_t>(sample - 1)] += value;
    }
  }
  const int close_status = gzclose(input);
  if (close_status != Z_OK || actual_nnz != declared_nnz) {
    std::cerr << "Matrix scan failed integrity check: declared=" << declared_nnz
              << " actual=" << actual_nnz << " gzclose=" << close_status << "\n";
    return 74;
  }

  uint64_t output_nnz = 0;
  for (const uint64_t value : sums) if (value > 0) ++output_nnz;
  gzFile output = gzopen(output_path.c_str(), "wb");
  if (!output) {
    std::cerr << "Cannot open output matrix: " << output_path << "\n";
    return 75;
  }
  gzprintf(output, "%%%%MatrixMarket matrix coordinate integer general\n");
  gzprintf(output, "%llu %d %llu\n", static_cast<unsigned long long>(rows), n_samples,
           static_cast<unsigned long long>(output_nnz));
  for (uint64_t row = 0; row < rows; ++row) {
    for (int32_t sample = 0; sample < n_samples; ++sample) {
      const uint64_t value = sums[(static_cast<size_t>(row) * static_cast<size_t>(n_samples)) + static_cast<size_t>(sample)];
      if (value > 0) {
        gzprintf(output, "%llu %d %llu\n", static_cast<unsigned long long>(row + 1U), sample + 1,
                 static_cast<unsigned long long>(value));
      }
    }
  }
  if (gzclose(output) != Z_OK) {
    std::cerr << "Cannot finalize output matrix\n";
    return 76;
  }

  std::ofstream audit(audit_path);
  if (!audit) {
    std::cerr << "Cannot open audit: " << audit_path << "\n";
    return 77;
  }
  audit << "input_rows\tinput_cols\tdeclared_nnz\tactual_nnz\toutput_rows\toutput_cols\toutput_nnz\n";
  audit << rows << '\t' << cols << '\t' << declared_nnz << '\t' << actual_nnz << '\t'
        << rows << '\t' << n_samples << '\t' << output_nnz << '\n';
  std::cout << "STAGE4A_MALIGNANT_PSEUDOBULK_STREAM_OK declared_nnz=" << declared_nnz
            << " actual_nnz=" << actual_nnz << " output_nnz=" << output_nnz << "\n";
  return 0;
}
