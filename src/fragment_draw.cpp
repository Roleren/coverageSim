// Joint or per-position fragment draw without building the candidate table
// (carina-claude). Used by make_simulated_rpf_fragments() when no candidate
// weight depends on the fragment sequence: every candidate weight is then
// geometry probability x frame weight x position score. The arithmetic and the
// random draws follow the R path step by step (fragment_candidates(),
// draw_fragment_groups(), stats::rmultinom()), so the drawn fragments, their
// probabilities and the RNG state are identical; only drawn candidates are
// returned.
#include <Rcpp.h>
#include <Rmath.h>
#include <vector>
#include <cmath>
#include <algorithm>
#include <climits>
using namespace Rcpp;

namespace {

struct Pairs {
  const int* length;
  const int* offset;
  const double* prob;
  int n;
};

// Candidate k fits at transcript site `site` of a transcript of length `len`
// (fragment_candidates(): tx_start >= 1 and tx_end <= transcript length).
inline bool fits(int site, double len, const Pairs& pairs, int k) {
  int tx_start = site - pairs.offset[k];
  return tx_start >= 1 && (double)(tx_start + pairs.length[k] - 1) <= len;
}

// probability / sum(probability) of the fitting pairs of one position, as
// data.table's `by = signal_row` with base sum() (long double accumulation).
inline double position_sum(int site, double len, const Pairs& pairs) {
  long double s = 0.0L;
  for (int k = 0; k < pairs.n; ++k) {
    if (fits(site, len, pairs, k)) s += pairs.prob[k];
  }
  return (double)s;
}

}  // namespace

// [[Rcpp::export(rng = true)]]
List fragment_draw_cpp(IntegerVector site_tx, NumericVector tx_length,
                       IntegerVector frame, NumericVector score,
                       IntegerVector group, NumericVector total,
                       IntegerVector pair_length, IntegerVector pair_offset,
                       NumericVector pair_prob, NumericVector frame_weight) {
  const int n_rows = site_tx.size();
  Pairs pairs = {pair_length.begin(), pair_offset.begin(), pair_prob.begin(),
                 (int)pair_length.size()};
  const int n_groups = total.size();

  // Rows of each group, groups in order of first appearance (rows come in
  // signal-row order, as the candidate table does).
  std::vector<int> count(n_groups, 0), order;
  for (int i = 0; i < n_rows; ++i) {
    int g = group[i] - 1;
    if (count[g]++ == 0) order.push_back(g);
  }
  std::vector<int> start(n_groups + 1, 0);
  for (int g = 0; g < n_groups; ++g) start[g + 1] = start[g] + count[g];
  std::vector<int> fill(start.begin(), start.end() - 1), rows(n_rows);
  for (int i = 0; i < n_rows; ++i) rows[fill[group[i] - 1]++] = i;

  // Candidate weight: (geometry probability x bias weight) x score, with bias
  // weight = (1 x 1) x frame weight (no end or codon bias).
  auto weight = [&](int i, int k, double row_sum, double* gp_out, double* fw_out) {
    double gp = pairs.prob[k] / row_sum;
    double fw = frame_weight[frame[i] * pairs.n + k];
    if (gp_out) *gp_out = gp;
    if (fw_out) *fw_out = fw;
    return gp * fw * score[i];
  };

  // As in the R path, every weight is checked before the first draw.
  for (int i = 0; i < n_rows; ++i) {
    double row_sum = position_sum(site_tx[i], tx_length[i], pairs);
    for (int k = 0; k < pairs.n; ++k) {
      if (!fits(site_tx[i], tx_length[i], pairs, k)) continue;
      if (!std::isfinite(weight(i, k, row_sum, nullptr, nullptr))) {
        stop("End-bias strength produces non-finite weights; reduce strength");
      }
    }
  }

  std::vector<int> out_row, out_pair, out_score;
  std::vector<double> out_p, out_gp, out_fw;
  std::vector<double> w, p, fixed;
  std::vector<int> cand_row, cand_pair, drawn;
  for (int g : order) {
    w.clear(); cand_row.clear(); cand_pair.clear();
    for (int r = start[g]; r < start[g + 1]; ++r) {
      int i = rows[r];
      double row_sum = position_sum(site_tx[i], tx_length[i], pairs);
      for (int k = 0; k < pairs.n; ++k) {
        if (!fits(site_tx[i], tx_length[i], pairs, k)) continue;
        w.push_back(weight(i, k, row_sum, nullptr, nullptr));
        cand_row.push_back(i); cand_pair.push_back(k);
      }
    }
    if (w.size() > (size_t)INT_MAX) stop("Too many fragment candidates in one draw group");
    const int K = (int)w.size();
    // R: if (sum(w) <= 0) stop(); w <- w / max(w); p <- w / sum(w)
    long double s = 0.0L;
    for (int j = 0; j < K; ++j) s += w[j];
    if ((double)s <= 0) stop("End-bias weights exclude every feasible fragment");
    double m = w[0];
    for (int j = 1; j < K; ++j) if (w[j] > m) m = w[j];
    p.resize(K);
    for (int j = 0; j < K; ++j) p[j] = w[j] / m;
    s = 0.0L;
    for (int j = 0; j < K; ++j) s += p[j];
    const double ps = (double)s;
    for (int j = 0; j < K; ++j) p[j] = p[j] / ps;
    // stats::rmultinom(1, size, p): FixupProb() rescales by the double sum of
    // the positive entries, then Rf_rmultinom() draws.
    double fsum = 0.0;
    for (int j = 0; j < K; ++j) if (p[j] > 0.0) fsum += p[j];
    if (fsum <= 0.0) stop("no positive probabilities");
    fixed.resize(K);
    for (int j = 0; j < K; ++j) fixed[j] = p[j] / fsum;
    double size = total[g];
    if (ISNAN(size) || size < 0 || size > INT_MAX) stop("invalid second argument 'size'");
    drawn.assign(K, 0);
    rmultinom((int)size, fixed.data(), K, drawn.data());
    for (int j = 0; j < K; ++j) {
      if (drawn[j] <= 0) continue;
      int i = cand_row[j], k = cand_pair[j];
      double gp, fw;
      double row_sum = position_sum(site_tx[i], tx_length[i], pairs);
      weight(i, k, row_sum, &gp, &fw);
      out_row.push_back(i + 1); out_pair.push_back(k + 1);
      out_score.push_back(drawn[j]); out_p.push_back(p[j]);
      out_gp.push_back(gp); out_fw.push_back(fw);
    }
  }
  // Drawn candidates in candidate-table order (signal row, then pair).
  std::vector<int> idx(out_row.size());
  for (size_t j = 0; j < idx.size(); ++j) idx[j] = j;
  std::stable_sort(idx.begin(), idx.end(), [&](int a, int b) {
    return out_row[a] != out_row[b] ? out_row[a] < out_row[b] : out_pair[a] < out_pair[b];
  });
  IntegerVector row_o(idx.size()), pair_o(idx.size()), score_o(idx.size());
  NumericVector p_o(idx.size()), gp_o(idx.size()), fw_o(idx.size());
  for (size_t j = 0; j < idx.size(); ++j) {
    row_o[j] = out_row[idx[j]]; pair_o[j] = out_pair[idx[j]];
    score_o[j] = out_score[idx[j]]; p_o[j] = out_p[idx[j]];
    gp_o[j] = out_gp[idx[j]]; fw_o[j] = out_fw[idx[j]];
  }
  return List::create(_["row"] = row_o, _["pair"] = pair_o, _["score"] = score_o,
                      _["probability"] = p_o, _["geometry_probability"] = gp_o,
                      _["frame_bias_weight"] = fw_o);
}

// Whether any length-and-offset pair fits at each transcript site
// (fragment_signal_fits(); NA sites and lengths never fit).
// [[Rcpp::export(rng = false)]]
LogicalVector fragment_fits_cpp(IntegerVector site_tx, NumericVector tx_length,
                                IntegerVector pair_length, IntegerVector pair_offset) {
  const int n = site_tx.size();
  Pairs pairs = {pair_length.begin(), pair_offset.begin(), nullptr,
                 (int)pair_length.size()};
  LogicalVector out(n);
  for (int i = 0; i < n; ++i) {
    bool any = false;
    if (site_tx[i] != NA_INTEGER && !ISNAN(tx_length[i])) {
      for (int k = 0; k < pairs.n && !any; ++k) any = fits(site_tx[i], tx_length[i], pairs, k);
    }
    out[i] = any;
  }
  return out;
}
