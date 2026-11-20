#include <Rcpp.h>
#include <Rmath.h>
#include <vector>
using namespace Rcpp;

// Uniform multinomial sampling for all regions in a batch. The R API's
// rmultinom uses R's RNG, so the caller's seed and RNG kind remain in force.
// [[Rcpp::export]]
IntegerVector sample_uniform_regions_cpp(IntegerVector widths,
                                          NumericVector budgets,
                                          int padding) {
  if (widths.size() != budgets.size() || padding < 0)
    stop("Invalid region dimensions");
  R_xlen_t total = 0;
  for (R_xlen_t i = 0; i < widths.size(); ++i) {
    if (widths[i] <= 0 || !R_FINITE(budgets[i]) || budgets[i] < 0 ||
        budgets[i] != std::floor(budgets[i]) || budgets[i] > INT_MAX)
      stop("Invalid multinomial region");
    total += static_cast<R_xlen_t>(widths[i]) + 2 * padding;
  }
  IntegerVector out(total);
  R_xlen_t offset = 0;
  for (R_xlen_t i = 0; i < widths.size(); ++i) {
    int width = widths[i];
    std::vector<double> prob(width, 1.0 / width);
    std::vector<int> draws(width);
    ::rmultinom(static_cast<int>(budgets[i]), prob.data(), width, draws.data());
    std::copy(draws.begin(), draws.end(), out.begin() + offset + padding);
    offset += static_cast<R_xlen_t>(width) + 2 * padding;
  }
  return out;
}
