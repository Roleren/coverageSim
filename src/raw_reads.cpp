#include <Rcpp.h>
#include <cmath>
#include <cstring>

// FASTQ records for simulate_raw_reads(): cut each read to read_length, draw
// Phred qualities, add substitution errors at the rate the quality implies
// (plus error_rate) and N for the lowest quality, and return the records as one
// raw vector ready to write. Uses R's random numbers, so set.seed() applies.
//
// Quality model: each read has its own level, quality falls along the read and
// about one base in eight comes from a low tail. profile 0 = perfect (Q40, no
// errors unless error_rate), 1 = novaseq (binned to 2/12/23/37), 2 = hiseq.
// insert_start / insert_end are 1-based positions of the insert in the read.
// [[Rcpp::export]]
Rcpp::List raw_read_records_cpp(Rcpp::CharacterVector reads, Rcpp::CharacterVector ids,
                                int read_length, int profile, double error_rate,
                                Rcpp::IntegerVector insert_start,
                                Rcpp::IntegerVector insert_end) {
  const R_xlen_t n = reads.size();
  if (ids.size() != n || insert_start.size() != n || insert_end.size() != n) {
    Rcpp::stop("reads, ids and insert positions must have the same length");
  }
  double error_probability[42];
  for (int q = 0; q < 42; ++q) {
    double p = profile == 0 ? 0.0 : std::pow(10.0, -q / 10.0);
    error_probability[q] = std::min(1.0, p + error_rate);
  }
  std::vector<double> decline(read_length);
  for (int j = 0; j < read_length; ++j) {
    const double position = read_length > 1 ? static_cast<double>(j) / (read_length - 1) : 0.0;
    decline[j] = 5.0 * position * position;
  }
  size_t total = 0;
  for (R_xlen_t i = 0; i < n; ++i) {
    if (static_cast<int>(LENGTH(STRING_ELT(reads, i))) < read_length) {
      Rcpp::stop("a read is shorter than read_length");
    }
    total += static_cast<size_t>(LENGTH(STRING_ELT(ids, i))) + 2 * read_length + 6;
  }
  Rcpp::RawVector records(total);
  Rcpp::IntegerVector n_errors(n), n_errors_insert(n);
  unsigned char* out = RAW(records);
  static const char alphabet[4] = {'A', 'C', 'G', 'T'};
  double sum_errors = 0;
  for (R_xlen_t i = 0; i < n; ++i) {
    SEXP id = STRING_ELT(ids, i);
    const int id_length = LENGTH(id);
    *out++ = '@';
    std::memcpy(out, CHAR(id), id_length);
    out += id_length;
    *out++ = '\n';
    unsigned char* sequence = out;
    std::memcpy(sequence, CHAR(STRING_ELT(reads, i)), read_length);
    out += read_length;
    *out++ = '\n';
    *out++ = '+';
    *out++ = '\n';
    unsigned char* quality = out;
    out += read_length;
    *out++ = '\n';
    double level = 0.0;
    if (profile != 0) {
      level = exp_rand() * 1.2;
      if (unif_rand() < 0.03) level += 10.0;
    }
    int errors = 0, errors_insert = 0;
    const int first = insert_start[i], last = insert_end[i];
    for (int j = 0; j < read_length; ++j) {
      int q = 40;
      if (profile != 0) {
        const double draw = unif_rand();
        double value = 39.5 - level - decline[j];
        if (draw < 0.12) value += std::log(draw / 0.12) * 6.0;
        q = static_cast<int>(std::lround(std::max(2.0, std::min(41.0, value))));
        if (profile == 1) q = q < 7 ? 2 : (q < 18 ? 12 : (q < 30 ? 23 : 37));
      }
      quality[j] = static_cast<unsigned char>(q + 33);
      bool wrong = false;
      if (profile != 0 && q <= 2) {
        sequence[j] = 'N';
        wrong = true;
      } else if (error_probability[q] > 0.0 && unif_rand() < error_probability[q]) {
        int old = -1;
        for (int k = 0; k < 4; ++k) if (sequence[j] == alphabet[k]) old = k;
        if (old >= 0) {
          int step = 1 + static_cast<int>(unif_rand() * 3.0);
          if (step > 3) step = 3;
          sequence[j] = alphabet[(old + step) % 4];
          wrong = true;
        }
      }
      if (wrong) {
        ++errors;
        if (j + 1 >= first && j + 1 <= last) ++errors_insert;
      }
    }
    n_errors[i] = errors;
    n_errors_insert[i] = errors_insert;
    sum_errors += errors;
    if ((i & 0xFFFF) == 0) Rcpp::checkUserInterrupt();
  }
  return Rcpp::List::create(Rcpp::Named("records") = records,
                            Rcpp::Named("n_errors") = n_errors,
                            Rcpp::Named("n_errors_insert") = n_errors_insert,
                            Rcpp::Named("errors") = sum_errors);
}
