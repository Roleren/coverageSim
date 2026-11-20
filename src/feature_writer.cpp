#include <Rcpp.h>
#include <string>
#include <vector>
#include <unordered_set>

// [[Rcpp::export]]
Rcpp::List write_feature_sequences_cpp(Rcpp::CharacterVector tx, Rcpp::IntegerVector leader,
    Rcpp::IntegerVector cds_length, Rcpp::IntegerVector gene, Rcpp::IntegerVector start,
    Rcpp::IntegerVector end, Rcpp::CharacterVector classes, Rcpp::CharacterVector starts, Rcpp::LogicalVector has_start,
    Rcpp::CharacterVector chosen_stops, Rcpp::CharacterVector stop_codons, int max_rounds) {
  const int n = tx.size(), nf = gene.size();
  std::string bases; std::vector<size_t> offsets(n), widths(n);
  for (int i=0; i<n; ++i) { offsets[i]=bases.size(); std::string s=Rcpp::as<std::string>(tx[i]); widths[i]=s.size(); bases+=s; }
  std::vector<unsigned char> fixed(bases.size(), 0), seen(bases.size(), 0);
  std::vector<size_t> cs(n), ce(n), fs(nf), fe(nf), codons;
  for (int i=0; i<n; ++i) { cs[i]=offsets[i]+leader[i]; ce[i]=cs[i]+cds_length[i]-1; for(int k=0;k<3;++k) fixed[cs[i]+k]=fixed[ce[i]-k]=1; }
  for (int i=0;i<nf;++i) { fs[i]=offsets[gene[i]-1]+start[i]-1; fe[i]=offsets[gene[i]-1]+end[i]-1; }
  for (int k=0;k<3;++k) {
    for (int i=0;i<nf;++i) if (has_start[i]) { std::string s=Rcpp::as<std::string>(starts[i]); bases[fs[i]+k]=s[k]; fixed[fs[i]+k]=1; }
    for (int i=0;i<nf;++i) if (Rcpp::as<std::string>(classes[i])!="NTE") { std::string s=Rcpp::as<std::string>(chosen_stops[i]); bases[fe[i]-2+k]=s[k]; fixed[fe[i]-2+k]=1; }
  }
  auto add_codons = [&](size_t a, size_t b) { for(size_t p=a+3; p+5<=b; p+=3) if(!seen[p]) {seen[p]=1; codons.push_back(p);} };
  for(int i=0;i<n;++i) add_codons(cs[i],ce[i]);
  // A feature without its own start codon (frameshift segment, readthrough
  // extension) is stop-free from its first codon on.
  for(int i=0;i<nf;++i) add_codons(has_start[i] ? fs[i] : fs[i]-3, fe[i]);
  std::unordered_set<std::string> stops;
  for (auto s : stop_codons) if(s!=NA_STRING) stops.insert(Rcpp::as<std::string>(s));
  Rcpp::Function sample=Rcpp::Environment::base_env()["sample"];
  Rcpp::CharacterVector alphabet=Rcpp::CharacterVector::create("A","C","G","T");
  int round=0; bool has_bad=false;
  for(round=1;round<=max_rounds;++round) {
    std::vector<size_t> bad, positions;
    for(auto p: codons) if(stops.count(bases.substr(p,3))) bad.push_back(p);
    has_bad=!bad.empty(); if(!has_bad) break;
    std::fill(seen.begin(), seen.end(), 0);
    for(int k=0;k<3;++k) for(auto p:bad) { p+=k; if(!fixed[p] && !seen[p]) {seen[p]=1; positions.push_back(p);} }
    if(positions.empty()) Rcpp::stop("A stop codon inside a feature is made only of fixed start/stop bases; this placement is impossible (please report)");
    Rcpp::CharacterVector draws=sample(alphabet, positions.size(), Rcpp::Named("replace")=true);
    for(size_t i=0;i<positions.size();++i) bases[positions[i]]=CHAR(draws[i])[0];
  }
  if(has_bad) Rcpp::stop("Could not remove internal stop codons in %d rounds", max_rounds);
  Rcpp::CharacterVector result(n);
  for(int i=0;i<n;++i) result[i]=bases.substr(offsets[i],widths[i]);
  return Rcpp::List::create(Rcpp::Named("tx_seqs")=result,Rcpp::Named("rounds")=round);
}
