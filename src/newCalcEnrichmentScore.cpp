#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
double newCalcEnrichmentScoreCPP(IntegerVector Set, NumericVector Eso, double exponent, bool shuffle = false) {
  std::vector<int> sset=Rcpp::as< std::vector<int> >(Set);
  std::vector<double> eso=Rcpp::as< std::vector<double> >(Eso);

//  if(shuffle) {
//    std::random_shuffle(sset.begin(), sset.end());
//  }
  int N=sset.size();
  int nh=std::accumulate(sset.begin(), sset.end(), 0.0);

  if(nh == 0) {
    return 0;
  }

  double n=(-1.0)/((double)(N-nh));
  double nr=0;

  for(int j=0;j<N;j++) {
    if(sset[j]) {
      // nr += (eso[j] > 0 ? eso[j] : -eso[j]);
      nr += std::pow((eso[j] > 0 ? eso[j] : -eso[j]), exponent);
    }
  }

  double smax=0; double smin=0; double cs=0;
  for(int j=0;j<N;j++) {
    if(sset[j]) {
      // cs += (eso[j] > 0 ? eso[j] : -eso[j]) / nr;
      cs += std::pow((eso[j] > 0 ? eso[j] : -eso[j]), exponent) / nr;
    } else {
      cs += n;
    }
    if(cs>smax) {
      smax = cs;
    } else if(cs<smin) {
      smin = cs;
    }
  }
  return std::abs(smax) > std::abs(smin) ? smax : smin;
}

// [[Rcpp::export]]
double single_newCalcEnrichmentScoreCPP(IntegerVector Set, NumericVector Eso, double exponent, int times = 1000) {
  double es;
//  es = newCalcEnrichmentScoreCPP(Set, Eso, exponent);
  std::vector<int> sset=Rcpp::as< std::vector<int> >(Set);
  std::vector<double> eso=Rcpp::as< std::vector<double> >(Eso);

  int N=sset.size();
  int nh=std::accumulate(sset.begin(), sset.end(), 0.0);

  if(nh == 0) {
    return 0;
  }

  double n=(-1.0)/((double)(N-nh));
  double nr=0;

  for(int j=0;j<N;j++) {
    if(sset[j]) {
      eso[j] = std::pow((eso[j] > 0 ? eso[j] : -eso[j]), exponent);
      // nr += (eso[j] > 0 ? eso[j] : -eso[j]);
      nr += eso[j];
    } else {
      eso[j] = n;
    }
  }

  double smax=0; double smin=0; double cs=0;
  for(int j=0;j<N;j++) {
    if(sset[j]) {
      // cs += (eso[j] > 0 ? eso[j] : -eso[j]) / nr;
      eso[j] /= nr;
    }
    cs += eso[j];

    if(cs>smax) {
      smax = cs;
    } else if(cs<smin) {
      smin = cs;
    }
  }
  es = std::abs(smax) > std::abs(smin) ? smax : smin;

  int num = 0; double sum_es = 0;
  for(int i = 0; i < times; i++){
    double new_es;
    std::random_shuffle(eso.begin(), eso.end());
    smax = 0; smin = 0; cs = 0;
    for(int i = 0; i < eso.size(); i++){
      cs += eso[i];
      if(cs>smax) {
        smax = cs;
      } else if(cs<smin) {
        smin = cs;
      }
    }

    new_es = std::abs(smax) > std::abs(smin) ? smax : smin;
    if(new_es * es > 0){
      sum_es += new_es;
      num ++;
    }
  }
  if(num == 0) {
    return (sum_es/num);
  } else {
    return es * num / std::abs(sum_es);
  }
}
