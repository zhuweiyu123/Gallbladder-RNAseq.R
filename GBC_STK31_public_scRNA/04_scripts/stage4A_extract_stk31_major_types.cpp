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
              << " <matrix.mtx.gz> <cell_to_group.int32.bin> <n_rows> <n_cols>"
              << " <stk31_row_1based> <n_groups> <output.tsv>\n";
    return 64;
  }

  const std::string matrix_path = argv[1];
  const std::string map_path = argv[2];
  const uint64_t expected_rows = std::strtoull(argv[3], nullptr, 10);
  const uint64_t expected_cols = std::strtoull(argv[4], nullptr, 10);
  const uint64_t stk31_row = std::strtoull(argv[5], nullptr, 10);
  const int32_t n_groups = static_cast<int32_t>(std::strtol(argv[6], nullptr, 10));
  const std::string output_path = argv[7];

  std::ifstream map_in(map_path, std::ios::binary);
  if (!map_in) {
    std::cerr << "Cannot open map: " << map_path << "\n";
    return 65;
  }
  std::vector<int32_t> cell_to_group(expected_cols, 0);
  map_in.read(reinterpret_cast<char*>(cell_to_group.data()),
              static_cast<std::streamsize>(expected_cols * sizeof(int32_t)));
  if (map_in.gcount() != static_cast<std::streamsize>(expected_cols * sizeof(int32_t))) {
    std::cerr << "Cell-to-group map has unexpected byte length\n";
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
  if (dimension_line.empty()) {
    std::cerr << "Missing Matrix Market dimensions\n";
    gzclose(input);
    return 70;
  }

  uint64_t rows = 0, cols = 0, declared_nnz = 0;
  {
    std::istringstream parser(dimension_line);
    if (!(parser >> rows >> cols >> declared_nnz) || rows != expected_rows || cols != expected_cols) {
      std::cerr << "Unexpected matrix dimensions: " << dimension_line;
      gzclose(input);
      return 71;
    }
  }

  std::vector<uint64_t> library_umi(static_cast<size_t>(n_groups) + 1U, 0);
  std::vector<uint64_t> stk31_umi(static_cast<size_t>(n_groups) + 1U, 0);
  std::vector<uint64_t> stk31_detected(static_cast<size_t>(n_groups) + 1U, 0);
  uint64_t actual_nnz = 0;

  while (gzgets(input, buffer.data(), buffer_size)) {
    uint64_t row = 0, col = 0, value = 0;
    std::istringstream parser(buffer.data());
    if (!(parser >> row >> col >> value)) {
      std::cerr << "Malformed coordinate at nonzero " << (actual_nnz + 1U) << "\n";
      gzclose(input);
      return 72;
    }
    ++actual_nnz;
    if (row == 0 || row > rows || col == 0 || col > cols) {
      std::cerr << "Out-of-range coordinate at nonzero " << actual_nnz << "\n";
      gzclose(input);
      return 73;
    }
    const int32_t group = cell_to_group[static_cast<size_t>(col - 1U)];
    if (group < 0 || group > n_groups) {
      std::cerr << "Invalid group index in map\n";
      gzclose(input);
      return 74;
    }
    if (group > 0) {
      library_umi[static_cast<size_t>(group)] += value;
      if (row == stk31_row && value > 0) {
        stk31_umi[static_cast<size_t>(group)] += value;
        ++stk31_detected[static_cast<size_t>(group)];
      }
    }
  }
  const int close_status = gzclose(input);
  if (close_status != Z_OK || actual_nnz != declared_nnz) {
    std::cerr << "Matrix scan failed integrity check: declared=" << declared_nnz
              << " actual=" << actual_nnz << " gzclose=" << close_status << "\n";
    return 75;
  }

  std::ofstream output(output_path);
  if (!output) {
    std::cerr << "Cannot open output: " << output_path << "\n";
    return 76;
  }
  output << "group_index\tgroup_library_UMI\tSTK31_raw_UMI\tSTK31_detected_cell_n\tdeclared_nnz\tactual_nnz\n";
  for (int32_t group = 1; group <= n_groups; ++group) {
    output << group << '\t' << library_umi[static_cast<size_t>(group)] << '\t'
           << stk31_umi[static_cast<size_t>(group)] << '\t'
           << stk31_detected[static_cast<size_t>(group)] << '\t'
           << declared_nnz << '\t' << actual_nnz << '\n';
  }
  std::cout << "STAGE4A_STK31_MAJOR_TYPE_STREAM_OK declared_nnz=" << declared_nnz
            << " actual_nnz=" << actual_nnz << "\n";
  return 0;
}
