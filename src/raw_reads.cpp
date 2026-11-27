#include <Rcpp.h>
#include "raw_read_id.h"
#include "raw_read_buffer.h"
#include <array>
#include <cmath>
#include <cstdint>
#include <cstring>

namespace {

constexpr char alphabet[] = "ACGT";
constexpr double uint32_range = 4294967296.0;

// SplitMix64 output function, used both as a sequential stream and to key
// reproducible shifts for copies of a pre-existing footprint.
uint64_t mix_bits(uint64_t value) {
  value = (value ^ (value >> 30)) * 0xbf58476d1ce4e5b9ULL;
  value = (value ^ (value >> 27)) * 0x94d049bb133111ebULL;
  return value ^ (value >> 31);
}

class RawReadModel {
public:
  RawReadModel(int length, int profile, double error_rate) :
    length(length), profile(profile), perfect(profile == 0 && error_rate == 0),
    decline(length), tail_multiplier(length) {
    for (int q = 0; q < 42; ++q) {
      const double p = profile == 0 ? 0.0 : std::pow(10.0, -q / 10.0);
      error_threshold[q] = std::min(1.0, p + error_rate) * uint32_range;
    }
    for (int j = 0; j < length; ++j) {
      const double position = length > 1 ? static_cast<double>(j) / (length - 1) : 0.0;
      decline[j] = 5.0 * position * position;
      tail_multiplier[j] = std::exp(decline[j] / 6.0);
    }
    if (!perfect) {
      const auto high = static_cast<uint64_t>(unif_rand() * uint32_range);
      state = (high << 32) | static_cast<uint64_t>(unif_rand() * uint32_range);
    }
  }

  void simulate(unsigned char* sequence, unsigned char* quality, int first, int end,
                int& errors, int& errors_insert) {
    errors = 0;
    errors_insert = 0;
    if (perfect) {
      std::memset(quality, 'I', length);
      return;
    }

    // The read-wide Q penalty is exponential with a 3% low-quality tail.
    double level = 0.0;
    if (profile != 0) {
      const uint64_t draw = random_bits();
      level = -1.2 * std::log((static_cast<double>(draw >> 32) + 0.5) / uint32_range);
      if (static_cast<uint32_t>(draw) < 0.03 * uint32_range) level += 10.0;
    }
    const double baseline = 39.5 - level;
    const double tail = profile == 1 ?
      0.12 * uint32_range * std::exp((29.5 - baseline) / 6.0) : 0.0;
    const double last_cutoff = tail * tail_multiplier.back();
    const double error37 = error_threshold[37];
    constexpr double candidate_limit = uint32_range / 8.0;
    if (profile == 1 && baseline - decline.back() >= 29.5 &&
        last_cutoff + (1 - last_cutoff / uint32_range) * error37 <= candidate_limit) {
      // At most one eighth of cycles can differ from a correct Q37 base.
      // AND three independent words to draw Bernoulli(1/8) candidates, then
      // thin these candidates to each cycle's exact quality/error probability.
      std::memset(quality, 'F', length);
      for (int block = 0; block < length;) {
        const int count = std::min(64, length - block);
        uint64_t candidates = random_bits();
        candidates &= random_bits();
        candidates &= random_bits();
        if (count < 64) candidates &= (uint64_t{1} << count) - 1;
        while (candidates) {
          const int j = block + __builtin_ctzll(candidates);
          candidates &= candidates - 1;
          const uint64_t draw = random_bits();
          const double u = (static_cast<double>(draw >> 32) + 0.5) / 8.0;
          const double cutoff = tail * tail_multiplier[j];
          if (u < cutoff) {
            const int q = nova_quality(u, baseline - decline[j], cutoff);
            quality[j] = static_cast<unsigned char>(q + 33);
            if (!apply_error(sequence[j], q, static_cast<uint32_t>(draw))) continue;
          } else if (u >= cutoff + (1 - cutoff / uint32_range) * error37 ||
                     !substitute(sequence[j])) {
            continue;
          }
          ++errors;
          if (j + 1 >= first && j + 1 <= end) ++errors_insert;
        }
        block += count;
      }
      return;
    }

    for (int j = 0; j < length; ++j) {
      const uint64_t draw = random_bits();
      const double quality_draw = static_cast<double>(draw >> 32) + 0.5;
      int q = 40;
      if (profile == 1) {
        q = nova_quality(quality_draw, baseline - decline[j], tail * tail_multiplier[j]);
      } else if (profile == 2) {
        double value = baseline - decline[j];
        if (quality_draw < 0.12 * uint32_range) {
          value += std::log(quality_draw / (0.12 * uint32_range)) * 6.0;
        }
        q = static_cast<int>(std::max(2.0, std::min(41.0, value)) + 0.5);
      }
      quality[j] = static_cast<unsigned char>(q + 33);
      if (!apply_error(sequence[j], q, static_cast<uint32_t>(draw))) continue;
      ++errors;
      if (j + 1 >= first && j + 1 <= end) ++errors_insert;
    }
  }

  const int length;
  const int profile;

private:
  static int nova_quality(double draw, double base, double tail) {
    // Invert the low-quality log tail at the three rounding boundaries.
    // If the unperturbed Q is below a boundary, the tail cannot cross it.
    static const double ratio23 = std::exp((17.5 - 29.5) / 6.0);
    static const double ratio12 = std::exp((6.5 - 29.5) / 6.0);
    if (base >= 29.5 && draw >= tail) return 37;
    if (base >= 17.5 && draw >= tail * ratio23) return 23;
    if (base >= 6.5 && draw >= tail * ratio12) return 12;
    return 2;
  }

  bool substitute(unsigned char& base) {
    const int old = base == 'A' ? 0 : base == 'C' ? 1 : base == 'G' ? 2 : base == 'T' ? 3 : -1;
    if (old < 0) return false;
    const auto next = static_cast<uint32_t>(random_bits());
    const unsigned int step = 1 + ((uint64_t{3} * next) >> 32);
    base = alphabet[(old + step) & 3];
    return true;
  }

  bool apply_error(unsigned char& base, int q, uint32_t draw) {
    if (q == 2) {
      base = 'N';
      return true;
    }
    return draw < error_threshold[q] && substitute(base);
  }

  uint64_t random_bits() {
    // SplitMix64. Seed once from R, then keep the stream independent of gzip.
    return mix_bits(state += 0x9e3779b97f4a7c15ULL);
  }

  const bool perfect;
  std::vector<double> decline, tail_multiplier;
  std::array<double, 42> error_threshold;
  uint64_t state = 0;
};

char complement(unsigned char base) {
  // The full Biostrings DNA alphabet, including ambiguity and mask letters.
  static constexpr char letters[] = "ACGTMRSVWYHKDBN-+.";
  static constexpr char complements[] = "TGCAKYSBWRDMHVN-+.";
  static const auto table = [] {
    std::array<char, 256> result{};
    for (size_t i = 0; i < sizeof(letters) - 1; ++i) result[letters[i]] = complements[i];
    return result;
  }();
  return table[base];
}

struct Bases {
  const unsigned char* data;
  size_t size;
};

Bases bases_from_string(SEXP value) {
  return {reinterpret_cast<const unsigned char*>(CHAR(value)), static_cast<size_t>(LENGTH(value))};
}

void append_bases(unsigned char*& out, size_t& remaining, Bases bases,
                  bool reverse, bool encoded = false) {
  const size_t take = std::min(remaining, bases.size);
  if (take == 0) return;
  // Biostrings DNA encoding packs the four canonical base flags in bits 0–3.
  static constexpr char decoded[] = "-ACMGRSVTWYHKDBN";
  auto decode = [&](size_t index) {
    const unsigned char code = bases.data[index];
    return code < 16 ? decoded[code] : code == 16 ? '-' : code == 32 ? '+' : '.';
  };
  if (encoded && reverse) {
    for (size_t i = 0; i < take; ++i) out[i] = complement(decode(bases.size - 1 - i));
  } else if (encoded) {
    for (size_t i = 0; i < take; ++i) out[i] = decode(i);
  } else if (reverse) {
    for (size_t i = 0; i < take; ++i) out[i] = complement(bases.data[bases.size - 1 - i]);
  } else {
    std::memcpy(out, bases.data, take);
  }
  out += take;
  remaining -= take;
}

}  // namespace

// [[Rcpp::export]]
SEXP raw_read_model_cpp(int read_length, int profile, double error_rate) {
  return Rcpp::XPtr<RawReadModel>(new RawReadModel(read_length, profile, error_rate), true);
}

// [[Rcpp::export(rng = false)]]
SEXP raw_read_buffer_cpp() {
  return Rcpp::XPtr<RawReadBuffer>(new RawReadBuffer, true);
}

// [[Rcpp::export(rng = false)]]
void raw_release_buffer_cpp(SEXP buffer) {
  Rcpp::XPtr<RawReadBuffer>(buffer).release();
}

// [[Rcpp::export(rng = false)]]
Rcpp::List raw_read_records_cpp(SEXP model_ptr, SEXP buffer_ptr, std::string id_prefix,
                                int id_width, int first, int n, Rcpp::List pieces,
                                Rcpp::IntegerVector polya, Rcpp::CharacterVector linker_arg,
                                Rcpp::CharacterVector adapter_arg, Rcpp::List references,
                                Rcpp::List inserts, bool reverse) {
  auto& model = *Rcpp::XPtr<RawReadModel>(model_ptr);
  const RawReadId ids(id_prefix, id_width);
  const Bases linker = bases_from_string(STRING_ELT(linker_arg, 0));
  const Bases adapter = bases_from_string(STRING_ELT(adapter_arg, 0));
  const size_t stride = ids.size() + 2 * static_cast<size_t>(model.length) + 6;
  auto& buffer = *Rcpp::XPtr<RawReadBuffer>(buffer_ptr);
  buffer.resize(stride * n);
  Rcpp::IntegerVector insert_len = Rcpp::no_init(n);
  Rcpp::IntegerVector insert_start = Rcpp::no_init(n), insert_end = Rcpp::no_init(n);
  Rcpp::IntegerVector adapter_start = Rcpp::no_init(n);
  Rcpp::IntegerVector n_errors = Rcpp::no_init(n), n_errors_insert = Rcpp::no_init(n);
  unsigned char* out = buffer.data();
  std::array<SEXP, 6> columns;
  std::array<Bases, 6> constant{};
  for (size_t k = 0; k < columns.size(); ++k) {
    columns[k] = VECTOR_ELT(pieces, k);
    if (XLENGTH(columns[k]) == 1) constant[k] = bases_from_string(STRING_ELT(columns[k], 0));
  }
  const int* tails = INTEGER(polya);
  const bool one_tail = polya.size() == 1;
  std::vector<const unsigned char*> reference_codes;
  for (R_xlen_t k = 0; k < references.size(); ++k) {
    const SEXP sequence = VECTOR_ELT(references, k);
    reference_codes.push_back(sequence == R_NilValue ? nullptr : RAW(sequence));
  }
  const SEXP strings = VECTOR_ELT(inserts, 0);
  const bool one_string = XLENGTH(strings) == 1;
  const SEXP reference_index = VECTOR_ELT(inserts, 1);
  const int* reference = INTEGER(reference_index);
  const bool one_reference = XLENGTH(reference_index) == 1;
  const double* offset = REAL(VECTOR_ELT(inserts, 2));
  const int* left = INTEGER(VECTOR_ELT(inserts, 3));
  const int* right = INTEGER(VECTOR_ELT(inserts, 4));
  double sum_errors = 0;
  for (R_xlen_t i = 0; i < n; ++i) {
    *out++ = '@';
    ids.write(reinterpret_cast<char*>(out), first + i);
    out += ids.size();
    *out++ = '\n';
    unsigned char* const sequence = out;
    size_t remaining = model.length;
    std::array<Bases, 6> part;
    for (size_t k = 0; k < columns.size(); ++k) {
      part[k] = constant[k].data ? constant[k] : bases_from_string(STRING_ELT(columns[k], i));
    }
    const int which_reference = reference[one_reference ? 0 : i];
    const bool encoded = which_reference > 0;
    const Bases insert = encoded ?
      Bases{reference_codes[which_reference - 1] + static_cast<R_xlen_t>(offset[i]),
            static_cast<size_t>(right[i] - left[i] + 1)} :
      bases_from_string(STRING_ELT(strings, one_string ? 0 : i));
    const int tail_bases = tails[one_tail ? 0 : i];
    insert_len[i] = insert.size;
    insert_start[i] = part[0].size + part[1].size + part[2].size + 1;
    insert_end[i] = insert_start[i] + insert_len[i] - 1;
    adapter_start[i] = insert_end[i] + part[3].size + tail_bases +
      part[4].size + part[5].size + 1;

    if (reverse) {
      append_bases(out, remaining, linker, true);
      append_bases(out, remaining, part[5], true);
      append_bases(out, remaining, part[4], true);
    } else {
      for (size_t k = 0; k < 3; ++k) append_bases(out, remaining, part[k], false);
      append_bases(out, remaining, insert, false, encoded);
      append_bases(out, remaining, part[3], false);
    }
    const size_t tail_length = std::min(remaining, static_cast<size_t>(tail_bases));
    std::memset(out, reverse ? 'T' : 'A', tail_length);
    out += tail_length;
    remaining -= tail_length;
    if (reverse) {
      append_bases(out, remaining, part[3], true);
      append_bases(out, remaining, insert, true, encoded);
      for (int k = 2; k >= 0; --k) append_bases(out, remaining, part[k], true);
    } else {
      append_bases(out, remaining, part[4], false);
      append_bases(out, remaining, part[5], false);
      append_bases(out, remaining, linker, false);
    }
    append_bases(out, remaining, adapter, false);
    std::memset(out, model.profile == 1 ? 'G' : 'A', remaining);
    out += remaining;
    *out++ = '\n';
    *out++ = '+';
    *out++ = '\n';
    int errors, errors_insert;
    model.simulate(sequence, out, insert_start[i], insert_end[i], errors, errors_insert);
    out += model.length;
    *out++ = '\n';
    n_errors[i] = errors;
    n_errors_insert[i] = errors_insert;
    sum_errors += errors;
    if ((i & 0xFFFF) == 0) Rcpp::checkUserInterrupt();
  }
  return Rcpp::List::create(Rcpp::Named("insert_len") = insert_len,
                            Rcpp::Named("insert_start_read") = insert_start,
                            Rcpp::Named("insert_end_read") = insert_end,
                            Rcpp::Named("adapter_start_read") = adapter_start,
                            Rcpp::Named("n_errors") = n_errors,
                            Rcpp::Named("n_errors_insert") = n_errors_insert,
                            Rcpp::Named("errors") = sum_errors);
}

// [[Rcpp::export]]
Rcpp::List raw_shift_ends_cpp(Rcpp::RawVector codes, Rcpp::NumericVector offset,
                              Rcpp::IntegerVector width, Rcpp::IntegerVector seq_idx,
                              Rcpp::IntegerVector start, Rcpp::IntegerVector end,
                              Rcpp::NumericVector inner5, Rcpp::NumericVector outer5,
                              Rcpp::NumericVector inner3, Rcpp::NumericVector outer3,
                              int max_shift, int min_length,
                              Rcpp::Nullable<Rcpp::IntegerVector> copy_rank,
                              Rcpp::Nullable<Rcpp::IntegerVector> shift_key) {
  const R_xlen_t n = start.size();
  Rcpp::IntegerVector new_start = Rcpp::no_init(n), new_end = Rcpp::no_init(n);
  const auto* sequence = RAW(codes);
  std::vector<double> cumulative(2 * static_cast<size_t>(max_shift) + 1);
  auto choose = [&](int position, int direction, R_xlen_t lower, int size, double draw,
                    const double* inner, const double* outer) {
    // Only candidates within the footprint's sequence are eligible.
    const int first = std::max(-max_shift, -position);
    const int last = std::min(max_shift, size - 1 - position);
    double sum = 0.0;
    for (int shift = first; shift <= last; ++shift) {
      const int site = position + shift, flank = site + direction;
      const auto outside = flank >= 0 && flank < size ? sequence[lower + flank] : 16;
      double weight = inner[std::min(16, static_cast<int>(sequence[lower + site])) - 1] * outer[std::min(16, outside) - 1];
      if (shift == 0) weight = std::max(weight, 1e-12);
      cumulative[shift - first] = (sum += weight);
    }
    const auto begin = cumulative.begin();
    return first + static_cast<int>(std::lower_bound(begin, begin + last - first + 1,
                                                     draw * sum) - begin);
  };
  Rcpp::IntegerVector rank;
  uint64_t key = 0;
  if (shift_key.isNotNull()) {
    rank = Rcpp::IntegerVector(copy_rank.get());
    const auto words = Rcpp::IntegerVector(shift_key.get());
    key = (static_cast<uint64_t>(words[0]) << 32) | static_cast<uint32_t>(words[1]);
  }
  for (R_xlen_t i = 0; i < n; ++i) {
    double draw5, draw3;
    if (key) {
      // A copy selected repeatedly with replacement retains its shifted ends.
      const uint64_t draw = mix_bits(key + static_cast<uint32_t>(rank[i]));
      draw5 = (static_cast<double>(draw >> 32) + 0.5) / uint32_range;
      draw3 = (static_cast<double>(static_cast<uint32_t>(draw)) + 0.5) / uint32_range;
    } else {
      draw5 = unif_rand();
      draw3 = unif_rand();
    }
    const int idx = seq_idx[i] - 1, size = width[idx];
    const auto lower = static_cast<R_xlen_t>(offset[idx]);
    const int shift5 = choose(start[i] - 1, -1, lower, size, draw5,
                              REAL(inner5), REAL(outer5));
    const int shift3 = choose(end[i] - 1, 1, lower, size, draw3,
                              REAL(inner3), REAL(outer3));
    const bool usable = end[i] + shift3 - (start[i] + shift5) + 1 >= min_length;
    new_start[i] = start[i] + (usable ? shift5 : 0);
    new_end[i] = end[i] + (usable ? shift3 : 0);
  }
  return Rcpp::List::create(Rcpp::Named("start") = new_start, Rcpp::Named("end") = new_end);
}

// [[Rcpp::export]]
Rcpp::List raw_draw_fragments_cpp(int n, Rcpp::IntegerVector usable, Rcpp::IntegerVector width,
                                  Rcpp::IntegerVector lengths, Rcpp::NumericVector source_cdf,
                                  Rcpp::NumericVector length_cdf, int pad) {
  Rcpp::IntegerVector seq_idx = Rcpp::no_init(n), start = Rcpp::no_init(n), end = Rcpp::no_init(n);
  const double* source_prob = REAL(source_cdf);
  const double* length_prob = REAL(length_cdf);
  for (int i = 0; i < n; ++i) {
    const int seq = usable[std::upper_bound(source_prob, source_prob + usable.size(), unif_rand()) -
      source_prob];
    const int seq_width = width[seq - 1];
    const int len = std::min(
      lengths[std::upper_bound(length_prob, length_prob + lengths.size(), unif_rand()) - length_prob],
      seq_width - 2 * pad - 1);
    seq_idx[i] = seq;
    start[i] = pad + static_cast<int>(std::ceil(unif_rand() * (seq_width - len - 2 * pad)));
    end[i] = start[i] + len - 1;
  }
  return Rcpp::List::create(Rcpp::Named("seq_idx") = seq_idx,
                            Rcpp::Named("start") = start, Rcpp::Named("end") = end);
}

// [[Rcpp::export(rng = false)]]
Rcpp::IntegerVector raw_permuted_indices_cpp(int first, int size, int total,
                                              Rcpp::IntegerVector key) {
  // Cycle-walk a keyed Feistel permutation on the smallest even number of
  // bits covering the input. Every supplied footprint is visited once.
  uint32_t bits = 1;
  while ((uint64_t{1} << (2 * bits)) < static_cast<uint32_t>(total)) ++bits;
  const uint32_t mask = (uint32_t{1} << bits) - 1;
  Rcpp::IntegerVector indices = Rcpp::no_init(size);
  for (int i = 0; i < size; ++i) {
    uint32_t rank = static_cast<uint32_t>(first) + i - 1;
    do {
      uint32_t left = rank >> bits;
      uint32_t right = rank & mask;
      for (int round = 0; round < 4; ++round) {
        uint32_t mixed = (right ^ static_cast<uint32_t>(key[round])) * 0x9e3779b9u;
        mixed ^= mixed >> 16;
        const uint32_t next = left ^ (mixed & mask);
        left = right;
        right = next;
      }
      rank = (left << bits) | right;
    } while (rank >= static_cast<uint32_t>(total));
    indices[i] = static_cast<int>(rank + 1);
  }
  return indices;
}
