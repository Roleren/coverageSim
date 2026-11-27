#include <Rcpp.h>
#include "raw_read_id.h"
#include "raw_read_buffer.h"
#include <charconv>
#include <algorithm>
#include <array>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <future>
#include <memory>
#include <stdexcept>
#include <vector>
#if defined(COVERAGESIM_ISAL)
#include <isa-l/igzip_lib.h>
#elif defined(COVERAGESIM_LIBDEFLATE)
#include <libdeflate.h>
#else
#include <zlib.h>
#endif

namespace {

class GzipCompressor {
public:
#if defined(COVERAGESIM_ISAL)
  void compress(const unsigned char* input, size_t length, RawReadBuffer& output) {
    // Reserve space for a gzip header/trailer and any stored DEFLATE blocks.
    output.resize(length + 5 * (length / 65535 + 1) + 18);
    isal_deflate_stateless_init(&stream);
    stream.gzip_flag = IGZIP_GZIP;
    stream.end_of_stream = 1;
    stream.level = 1;
    stream.level_buf = workspace.data();
    stream.level_buf_size = workspace.size();
    stream.next_in = const_cast<unsigned char*>(input);
    stream.avail_in = length;
    stream.next_out = output.data();
    stream.avail_out = output.size();
    if (isal_deflate_stateless(&stream) != COMP_OK) {
      throw std::runtime_error("read output gzip compression failed");
    }
    output.resize(stream.total_out);
  }

private:
  isal_zstream stream;
  std::array<unsigned char, ISAL_DEF_LVL1_DEFAULT> workspace;
#elif defined(COVERAGESIM_LIBDEFLATE)
  GzipCompressor() : compressor(libdeflate_alloc_compressor(1), libdeflate_free_compressor) {
    if (!compressor) throw std::runtime_error("cannot initialize read output gzip compressor");
  }

  void compress(const unsigned char* input, size_t length, RawReadBuffer& output) {
    output.resize(libdeflate_gzip_compress_bound(compressor.get(), length));
    const size_t written = libdeflate_gzip_compress(compressor.get(), input, length,
                                                   output.data(), output.size());
    if (!written) throw std::runtime_error("read output gzip compression failed");
    output.resize(written);
  }

private:
  std::unique_ptr<libdeflate_compressor, decltype(&libdeflate_free_compressor)>
    compressor{nullptr, libdeflate_free_compressor};
#else
  GzipCompressor() {
    if (deflateInit2(&stream, 1, Z_DEFLATED, 15 + 16, 8, Z_DEFAULT_STRATEGY) != Z_OK) {
      throw std::runtime_error("cannot initialize read output gzip compressor");
    }
  }

  ~GzipCompressor() { deflateEnd(&stream); }

  void compress(const unsigned char* input, size_t length, RawReadBuffer& output) {
    if (deflateReset(&stream) != Z_OK) throw std::runtime_error("read output gzip reset failed");
    output.resize(deflateBound(&stream, length));
    stream.next_in = const_cast<unsigned char*>(input);
    stream.avail_in = static_cast<uInt>(length);
    stream.next_out = reinterpret_cast<unsigned char*>(output.data());
    stream.avail_out = static_cast<uInt>(output.size());
    if (deflate(&stream, Z_FINISH) != Z_STREAM_END) {
      throw std::runtime_error("read output gzip compression failed");
    }
    output.resize(stream.total_out);
  }

private:
  z_stream stream{};
#endif
};

class RawReadWriter {
public:
  RawReadWriter(const std::string& path, bool gzip, int threads) :
    path(path), file(std::filesystem::u8path(path), std::ios::binary), gzip(gzip),
    threads(threads) {
    if (!file) throw std::runtime_error("cannot open read output: " + path);
  }

  void write(const unsigned char* data, size_t length) {
    if (!gzip) {
      file.write(reinterpret_cast<const char*>(data), length);
    } else {
      // Independent gzip members let workers compress in parallel. Write
      // members in input order and keep only one R chunk in flight.
      constexpr size_t block_size = 1024 * 1024;
      const size_t blocks = (length + block_size - 1) / block_size;
      const size_t workers = std::min<size_t>(threads, blocks);
      compressed.resize(blocks);
      while (compressors.size() < workers) compressors.emplace_back(new GzipCompressor);
      auto run_worker = [&](size_t worker) {
        for (size_t block = worker; block < blocks; block += workers) {
          const size_t offset = block * block_size;
          compressors[worker]->compress(data + offset, std::min(block_size, length - offset),
                                        compressed[block]);
        }
      };
      std::vector<std::future<void>> jobs;
      for (size_t worker = 1; worker < workers; ++worker) {
        jobs.push_back(std::async(std::launch::async, run_worker, worker));
      }
      if (workers) run_worker(0);
      for (auto& job : jobs) job.get();
      for (const auto& member : compressed) {
        file.write(reinterpret_cast<const char*>(member.data()), member.size());
      }
    }
    if (!file) throw std::runtime_error("cannot write read output: " + path);
  }

  std::string& text_buffer(size_t capacity) {
    csv.clear();
    csv.reserve(capacity);
    return csv;
  }

  void close() {
    if (!file.is_open()) return;
    file.close();
    if (!file) throw std::runtime_error("cannot close read output: " + path);
  }

private:
  std::string path;
  std::ofstream file;
  bool gzip;
  int threads;
  std::vector<RawReadBuffer> compressed;
  std::string csv;
  std::vector<std::unique_ptr<GzipCompressor>> compressors;
};

}  // namespace

// [[Rcpp::export(rng = false)]]
SEXP raw_open_writer_cpp(Rcpp::CharacterVector path, bool gzip, int threads) {
  return Rcpp::XPtr<RawReadWriter>(
    new RawReadWriter(Rf_translateCharUTF8(STRING_ELT(path, 0)), gzip, threads), true);
}

// [[Rcpp::export(rng = false)]]
void raw_write_fastq_cpp(SEXP writer, SEXP buffer_ptr) {
  const auto& buffer = *Rcpp::XPtr<RawReadBuffer>(buffer_ptr);
  Rcpp::XPtr<RawReadWriter>(writer)->write(buffer.data(), buffer.size());
}

// [[Rcpp::export(rng = false)]]
void raw_close_writer_cpp(SEXP writer) {
  Rcpp::XPtr<RawReadWriter> stream(writer);
  if (!stream.get()) return;
  stream->close();
  stream.release();
}

// [[Rcpp::export(rng = false)]]
void raw_write_truth_cpp(SEXP writer, Rcpp::List columns, std::string id_prefix, int id_width,
                          int first, int n, std::string separator, bool logical01) {
  struct Column {
    SEXP strings;
    const int* values;
    bool logical;
    std::string prefix;
  };
  const char sep = separator[0];
#ifdef _WIN32
  constexpr const char* eol = "\r\n";
#else
  constexpr const char* eol = "\n";
#endif
  auto& output = *Rcpp::XPtr<RawReadWriter>(writer);
  auto& csv = output.text_buffer(static_cast<size_t>(n) * 128);
  auto append_text = [&](std::string& text, const char* value, size_t length) {
    if (length == 0) {
      text += "\"\"";
      return;
    }
    size_t pos = 0;
    while (pos < length && value[pos] != sep && value[pos] != '"' &&
           value[pos] != '\n' && value[pos] != '\r') ++pos;
    if (pos == length) {
      text.append(value, length);
      return;
    }
    text.push_back('"');
    text.append(value, pos);
    for (; pos < length; ++pos) {
      if (value[pos] == '"') text.push_back('"');
      text.push_back(value[pos]);
    }
    text.push_back('"');
  };
  auto append_field = [&](std::string& text, const Column& field, int index) {
    text.push_back(sep);
    if (field.strings != R_NilValue) {
      const SEXP value = STRING_ELT(field.strings, index);
      if (value != NA_STRING) append_text(text, CHAR(value), LENGTH(value));
      return;
    }
    const int value = field.values[index];
    if (value == NA_INTEGER) return;
    if (field.logical && !logical01) {
      text.append(value ? "TRUE" : "FALSE");
    } else {
      char digits[16];
      const auto end = std::to_chars(digits, digits + sizeof(digits), value).ptr;
      text.append(digits, end);
    }
  };

  std::vector<Column> fields;
  std::string constants;
  for (R_xlen_t k = 0; k < columns.size(); ++k) {
    SEXP column = VECTOR_ELT(columns, k);
    Column field{R_NilValue, nullptr, TYPEOF(column) == LGLSXP, {}};
    if (TYPEOF(column) == STRSXP) field.strings = column;
    else field.values = field.logical ? LOGICAL(column) : INTEGER(column);

    // Protocol fields that are constant in the chunk need escaping only once.
    if (XLENGTH(column) == 1) {
      append_field(constants, field, 0);
    } else {
      field.prefix.swap(constants);
      fields.push_back(std::move(field));
    }
  }
  constants.append(eol);
  if (first == 1) {
    csv.append("read_id");
    const SEXP names = Rf_getAttrib(columns, R_NamesSymbol);
    for (R_xlen_t k = 0; k < columns.size(); ++k) {
      csv.push_back(sep);
      const SEXP name = STRING_ELT(names, k);
      append_text(csv, CHAR(name), LENGTH(name));
    }
    csv.append(eol);
  }
  const RawReadId ids(id_prefix, id_width);
  const bool quote_id = id_prefix.find_first_of(std::string{sep, '"', '\n', '\r'}) != std::string::npos;
  std::string id(ids.size(), '\0');
  for (int row = 0; row < n; ++row) {
    ids.write(&id[0], first + row);
    if (quote_id) append_text(csv, id.data(), id.size());
    else csv.append(id);
    for (const auto& field : fields) {
      csv.append(field.prefix);
      append_field(csv, field, row);
    }
    csv.append(constants);
  }
  output.write(reinterpret_cast<const unsigned char*>(csv.data()), csv.size());
}
