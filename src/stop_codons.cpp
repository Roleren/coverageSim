#include <Rcpp.h>
#include <string>
#include <unordered_set>
#include <vector>

// Positions are one-based, as in R. Return the positions whose three bases
// form a stop codon, in their original order.
// [[Rcpp::export]]
Rcpp::IntegerVector feature_stop_positions(Rcpp::CharacterVector bases,
                                            Rcpp::IntegerVector positions,
                                            Rcpp::CharacterVector stop_codons) {
  std::unordered_set<std::string> stops;
  for (R_xlen_t i = 0; i < stop_codons.size(); ++i) {
    if (stop_codons[i] != NA_STRING) stops.insert(Rcpp::as<std::string>(stop_codons[i]));
  }
  std::vector<int> result;
  for (R_xlen_t i = 0; i < positions.size(); ++i) {
    const int p = positions[i] - 1;
    std::string codon;
    for (int k = 0; k < 3; ++k) {
      SEXP base = bases[p + k];
      codon += base == NA_STRING ? "NA" : CHAR(base);
    }
    if (stops.count(codon)) result.push_back(positions[i]);
  }
  return Rcpp::wrap(result);
}
