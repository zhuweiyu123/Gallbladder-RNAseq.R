#include <zlib.h>

#include <algorithm>
#include <cerrno>
#include <cinttypes>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <sys/stat.h>
#include <vector>

struct Triplet {
  uint32_t row;
  uint32_t col;
  uint64_t value;
};
static_assert(sizeof(Triplet) == 16, "Unexpected Triplet layout");

static void fail(const std::string &message) { throw std::runtime_error(message); }

static uint64_t file_size(const std::string &path) {
  struct stat st {};
  if (stat(path.c_str(), &st) != 0) fail("stat failed for " + path + ": " + std::strerror(errno));
  return static_cast<uint64_t>(st.st_size);
}

static std::vector<int32_t> read_map(const std::string &path, size_t expected) {
  std::ifstream in(path, std::ios::binary);
  if (!in) fail("Cannot open map: " + path);
  std::vector<int32_t> result(expected);
  in.read(reinterpret_cast<char *>(result.data()), static_cast<std::streamsize>(expected * sizeof(int32_t)));
  if (in.gcount() != static_cast<std::streamsize>(expected * sizeof(int32_t))) fail("Short map: " + path);
  char extra;
  if (in.read(&extra, 1)) fail("Map has trailing bytes: " + path);
  return result;
}

static void gz_write_all(gzFile output, const char *data, size_t length) {
  while (length > 0) {
    const unsigned int chunk = static_cast<unsigned int>(std::min<size_t>(length, 1U << 30));
    const int written = gzwrite(output, data, chunk);
    if (written <= 0) {
      int err = 0;
      const char *message = gzerror(output, &err);
      fail(std::string("gzwrite failed: ") + (message ? message : "unknown"));
    }
    data += written;
    length -= static_cast<size_t>(written);
  }
}

static void append_uint(std::vector<char> &buffer, size_t &position, uint64_t value) {
  char digits[32];
  size_t n = 0;
  do {
    digits[n++] = static_cast<char>('0' + value % 10);
    value /= 10;
  } while (value > 0);
  while (n > 0) buffer[position++] = digits[--n];
}

static void binary_to_mtx_gz(const std::string &body_path, const std::string &output_path,
                             uint32_t rows, uint32_t cols, uint64_t nnz) {
  const std::string partial = output_path + ".partial";
  std::ifstream input(body_path, std::ios::binary);
  if (!input) fail("Cannot open temporary body: " + body_path);
  gzFile output = gzopen(partial.c_str(), "wb6");
  if (!output) fail("Cannot open gzip output: " + partial);
  gzbuffer(output, 16U * 1024U * 1024U);

  const std::string header = "%%MatrixMarket matrix coordinate integer general\n" +
                             std::to_string(rows) + " " + std::to_string(cols) + " " +
                             std::to_string(nnz) + "\n";
  gz_write_all(output, header.data(), header.size());

  const size_t records_per_block = 1U << 20;
  std::vector<Triplet> records(records_per_block);
  std::vector<char> text(16U * 1024U * 1024U);
  size_t text_position = 0;
  uint64_t records_written = 0;

  while (input) {
    input.read(reinterpret_cast<char *>(records.data()),
               static_cast<std::streamsize>(records.size() * sizeof(Triplet)));
    const std::streamsize bytes = input.gcount();
    if (bytes % static_cast<std::streamsize>(sizeof(Triplet)) != 0) fail("Truncated temporary body");
    const size_t count = static_cast<size_t>(bytes) / sizeof(Triplet);
    for (size_t index = 0; index < count; ++index) {
      if (text.size() - text_position < 96) {
        gz_write_all(output, text.data(), text_position);
        text_position = 0;
      }
      append_uint(text, text_position, records[index].row);
      text[text_position++] = ' ';
      append_uint(text, text_position, records[index].col);
      text[text_position++] = ' ';
      append_uint(text, text_position, records[index].value);
      text[text_position++] = '\n';
      ++records_written;
    }
  }
  if (input.bad()) fail("Read error in temporary body");
  if (text_position > 0) gz_write_all(output, text.data(), text_position);
  if (records_written != nnz) fail("Temporary body record count differs from expected nnz");
  if (gzclose(output) != Z_OK) fail("Failed to finalize gzip output: " + partial);
  if (std::rename(partial.c_str(), output_path.c_str()) != 0)
    fail("Cannot rename completed output: " + std::string(std::strerror(errno)));
}

int main(int argc, char **argv) {
  try {
    if (argc != 13) {
      std::cerr << "Usage: extractor input.mtx.gz col_map.bin row_map.bin full_rows full_cols "
                   "malignant_cols target_rows full_output.mtx.gz target_output.mtx.gz "
                   "cell_library.tsv audit.txt\n";
      return 64;
    }

    const std::string input_path = argv[1];
    const std::string col_map_path = argv[2];
    const std::string row_map_path = argv[3];
    const uint32_t expected_rows = static_cast<uint32_t>(std::stoul(argv[4]));
    const uint32_t expected_cols = static_cast<uint32_t>(std::stoul(argv[5]));
    const uint32_t malignant_cols = static_cast<uint32_t>(std::stoul(argv[6]));
    const uint32_t target_rows = static_cast<uint32_t>(std::stoul(argv[7]));
    const std::string full_output_path = argv[8];
    const std::string target_output_path = argv[9];
    const std::string library_path = argv[10];
    const std::string audit_path = argv[11];
    const std::string temp_prefix = argv[12];

    if (file_size(col_map_path) != static_cast<uint64_t>(expected_cols) * 4U) fail("Column map size mismatch");
    if (file_size(row_map_path) != static_cast<uint64_t>(expected_rows) * 4U) fail("Row map size mismatch");
    const auto col_map = read_map(col_map_path, expected_cols);
    const auto row_map = read_map(row_map_path, expected_rows);

    std::vector<uint8_t> seen_cols(static_cast<size_t>(malignant_cols) + 1U, 0U);
    uint64_t mapped_cols = 0;
    for (int32_t value : col_map) {
      if (value < 0 || static_cast<uint32_t>(value) > malignant_cols) fail("Invalid malignant column map value");
      if (value > 0) {
        if (seen_cols[static_cast<size_t>(value)]) fail("Duplicated malignant column map value");
        seen_cols[static_cast<size_t>(value)] = 1U;
        ++mapped_cols;
      }
    }
    if (mapped_cols != malignant_cols) fail("Malignant column map is not complete") ;

    std::vector<uint8_t> seen_rows(static_cast<size_t>(target_rows) + 1U, 0U);
    uint64_t mapped_rows = 0;
    for (int32_t value : row_map) {
      if (value < 0 || static_cast<uint32_t>(value) > target_rows) fail("Invalid target row map value");
      if (value > 0) {
        if (seen_rows[static_cast<size_t>(value)]) fail("Duplicated target row map value");
        seen_rows[static_cast<size_t>(value)] = 1U;
        ++mapped_rows;
      }
    }
    if (mapped_rows != target_rows) fail("Target row map is not complete");

    gzFile input = gzopen(input_path.c_str(), "rb");
    if (!input) fail("Cannot open input gzip: " + input_path);
    gzbuffer(input, 32U * 1024U * 1024U);
    char header_line[4096];
    if (!gzgets(input, header_line, sizeof(header_line))) fail("Cannot read Matrix Market banner");
    std::string banner(header_line);
    while (!banner.empty() && (banner.back() == '\n' || banner.back() == '\r')) banner.pop_back();
    if (banner != "%%MatrixMarket matrix coordinate integer general") fail("Unexpected Matrix Market banner: " + banner);

    uint64_t declared_nnz = 0;
    uint32_t header_rows = 0, header_cols = 0;
    while (true) {
      if (!gzgets(input, header_line, sizeof(header_line))) fail("Cannot read Matrix Market dimensions");
      if (header_line[0] == '%' || header_line[0] == '\n' || header_line[0] == '\r') continue;
      unsigned long long parsed_nnz = 0;
      unsigned int parsed_rows = 0, parsed_cols = 0;
      if (std::sscanf(header_line, "%u %u %llu", &parsed_rows, &parsed_cols, &parsed_nnz) != 3)
        fail("Malformed Matrix Market dimension line");
      header_rows = parsed_rows;
      header_cols = parsed_cols;
      declared_nnz = parsed_nnz;
      break;
    }
    if (header_rows != expected_rows || header_cols != expected_cols)
      fail("Input Matrix Market dimensions differ from manifest");

    const std::string full_body_path = temp_prefix + ".full.body.bin";
    const std::string target_body_path = temp_prefix + ".target.body.bin";
    FILE *full_body = std::fopen(full_body_path.c_str(), "wb");
    FILE *target_body = std::fopen(target_body_path.c_str(), "wb");
    if (!full_body || !target_body) fail("Cannot open temporary binary body files");
    std::vector<char> full_file_buffer(16U * 1024U * 1024U);
    std::vector<char> target_file_buffer(4U * 1024U * 1024U);
    setvbuf(full_body, full_file_buffer.data(), _IOFBF, full_file_buffer.size());
    setvbuf(target_body, target_file_buffer.data(), _IOFBF, target_file_buffer.size());

    std::vector<uint64_t> cell_library(malignant_cols, 0U);
    std::vector<uint64_t> target_sums(target_rows, 0U);
    uint64_t actual_nnz = 0, full_nnz = 0, target_nnz = 0;
    uint64_t adjacent_duplicates = 0;
    bool sorted_col_row = true;
    uint32_t previous_row = 0, previous_col = 0;

    auto process_triplet = [&](uint64_t row64, uint64_t col64, uint64_t value) {
      if (row64 < 1 || row64 > expected_rows || col64 < 1 || col64 > expected_cols)
        fail("Matrix triplet index out of bounds at entry " + std::to_string(actual_nnz + 1));
      if (value == 0) fail("Explicit zero encountered in input Matrix Market file");
      const uint32_t row = static_cast<uint32_t>(row64);
      const uint32_t col = static_cast<uint32_t>(col64);
      if (actual_nnz > 0) {
        if (col < previous_col || (col == previous_col && row < previous_row)) sorted_col_row = false;
        if (col == previous_col && row == previous_row) ++adjacent_duplicates;
      }
      previous_row = row;
      previous_col = col;
      ++actual_nnz;

      const int32_t malignant_col = col_map[col - 1U];
      if (malignant_col > 0) {
        const Triplet full_record{row, static_cast<uint32_t>(malignant_col), value};
        if (std::fwrite(&full_record, sizeof(full_record), 1, full_body) != 1) fail("Failed writing full temporary body");
        ++full_nnz;
        uint64_t &library = cell_library[static_cast<size_t>(malignant_col - 1)];
        if (std::numeric_limits<uint64_t>::max() - library < value) fail("Cell library-size overflow");
        library += value;

        const int32_t target_row = row_map[row - 1U];
        if (target_row > 0) {
          const Triplet target_record{static_cast<uint32_t>(target_row), static_cast<uint32_t>(malignant_col), value};
          if (std::fwrite(&target_record, sizeof(target_record), 1, target_body) != 1) fail("Failed writing target temporary body");
          ++target_nnz;
          target_sums[static_cast<size_t>(target_row - 1)] += value;
        }
      }
    };

    std::vector<unsigned char> input_buffer(32U * 1024U * 1024U);
    uint64_t token_value = 0;
    uint64_t fields[3] = {0, 0, 0};
    int field_index = 0;
    bool in_number = false;
    while (true) {
      const int bytes = gzread(input, input_buffer.data(), static_cast<unsigned int>(input_buffer.size()));
      if (bytes < 0) {
        int err = 0;
        const char *message = gzerror(input, &err);
        fail(std::string("gzread failed: ") + (message ? message : "unknown"));
      }
      if (bytes == 0) break;
      for (int index = 0; index < bytes; ++index) {
        const unsigned char character = input_buffer[static_cast<size_t>(index)];
        if (character >= '0' && character <= '9') {
          const uint64_t digit = character - '0';
          if (token_value > (std::numeric_limits<uint64_t>::max() - digit) / 10U) fail("Integer overflow while parsing input");
          token_value = token_value * 10U + digit;
          in_number = true;
        } else if (character == ' ' || character == '\t' || character == '\r' || character == '\n') {
          if (in_number) {
            fields[field_index++] = token_value;
            token_value = 0;
            in_number = false;
            if (field_index == 3) {
              process_triplet(fields[0], fields[1], fields[2]);
              field_index = 0;
            }
          }
        } else {
          fail("Unexpected non-numeric byte in Matrix Market body");
        }
      }
    }
    if (in_number) {
      fields[field_index++] = token_value;
      if (field_index == 3) {
        process_triplet(fields[0], fields[1], fields[2]);
        field_index = 0;
      }
    }
    if (field_index != 0) fail("Incomplete Matrix Market triplet at EOF");
    int gzip_error = Z_OK;
    const char *gzip_message = gzerror(input, &gzip_error);
    if (gzip_error != Z_OK && gzip_error != Z_STREAM_END)
      fail(std::string("Input gzip did not end cleanly: ") + (gzip_message ? gzip_message : "unknown"));
    if (gzclose(input) != Z_OK) fail("Input gzip CRC/finalization check failed");
    if (std::fclose(full_body) != 0 || std::fclose(target_body) != 0) fail("Failed closing temporary body files");

    if (actual_nnz != declared_nnz) fail("Actual triplet count differs from Matrix Market declaration");
    if (adjacent_duplicates != 0) fail("Duplicate adjacent coordinates detected in input");
    if (!sorted_col_row) fail("Input is not sorted by column then row; duplicate safety cannot be guaranteed");

    binary_to_mtx_gz(full_body_path, full_output_path, expected_rows, malignant_cols, full_nnz);
    binary_to_mtx_gz(target_body_path, target_output_path, target_rows, malignant_cols, target_nnz);

    const std::string library_partial = library_path + ".partial";
    {
      std::ofstream library_output(library_partial);
      if (!library_output) fail("Cannot write cell library sizes");
      library_output << "malignant_col\tall_gene_raw_umi\n";
      for (uint32_t col = 1; col <= malignant_cols; ++col)
        library_output << col << '\t' << cell_library[col - 1U] << '\n';
    }
    if (std::rename(library_partial.c_str(), library_path.c_str()) != 0) fail("Cannot finalize cell library-size output");

    {
      std::ofstream audit(audit_path);
      if (!audit) fail("Cannot write extractor audit");
      audit << "input=" << input_path << '\n'
            << "input_rows=" << header_rows << '\n'
            << "input_cols=" << header_cols << '\n'
            << "declared_nnz=" << declared_nnz << '\n'
            << "actual_nnz=" << actual_nnz << '\n'
            << "input_sorted_col_then_row=" << (sorted_col_row ? "TRUE" : "FALSE") << '\n'
            << "adjacent_duplicate_coordinates=" << adjacent_duplicates << '\n'
            << "malignant_cols=" << malignant_cols << '\n'
            << "malignant_output_nnz=" << full_nnz << '\n'
            << "target_rows=" << target_rows << '\n'
            << "target_output_nnz=" << target_nnz << '\n'
            << "target_row_sums=";
      for (size_t index = 0; index < target_sums.size(); ++index) {
        if (index) audit << '|';
        audit << (index + 1) << ':' << target_sums[index];
      }
      audit << '\n'
            << "malignant_matrix_gz_bytes=" << file_size(full_output_path) << '\n'
            << "target_matrix_gz_bytes=" << file_size(target_output_path) << '\n';
    }

    std::remove(full_body_path.c_str());
    std::remove(target_body_path.c_str());
    std::cout << "STAGE3_STREAM_EXTRACTION_OK\n"
              << "actual_nnz=" << actual_nnz << '\n'
              << "malignant_output_nnz=" << full_nnz << '\n'
              << "target_output_nnz=" << target_nnz << '\n';
    return 0;
  } catch (const std::exception &error) {
    std::cerr << "FATAL: " << error.what() << '\n';
    return 1;
  }
}
