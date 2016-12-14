// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;
#include "random.h"

// This function is fast for large L, m and small p.
//' @export
// [[Rcpp::export(".sim_crossovers")]]
arma::vec sim_crossovers(const double L,
                         const int m,
                         const double p,
                         const bool obligate_chiasma,
                         const double Lstar)
{
  if(m == 0) { // no-interference model is a lot easier
    int n_xo;

    if (obligate_chiasma) {
      // rejection sampling to get at least one chiasma
      while ((n_xo = R::rpois(Lstar / 50)) == 0);

      n_xo = R::rbinom((double) n_xo, 0.5);
    }
    else {
      n_xo = R::rpois(L / 100);
    }
    // SORT HERE!
    // arma::vec tmp = arma::randu<arma::vec>(n_xo) * L;
    // arma::vec tmp = Rcpp::as<arma::vec>(runif(n_xo, 0.0, L));
    // std::sort(tmp.begin(), tmp.end());
    // return tmp;

    // return arma::sort(arma::randu<arma::vec>(n_xo)) * L; // not bad

    // This weiredly works, probably automatic conversion to arma::vec
    Rcpp::NumericVector tmp = Rcpp::runif(n_xo, 0.0, L);
    std::sort(tmp.begin(), tmp.end());
    return tmp;
  }

  int n_points, first, n_nichi, n_ichi;

  double lambda1 = Lstar / 50 * (double) (m + 1) * (1 - p);
  double lambda2 = Lstar / 50 * p;

  while(true) {
    // chiasma and intermediate points
    n_points = R::rpois(lambda1);

    // which point is the first chiasma?
    first = random_int(0, m);
    if (first > n_points) n_ichi = 0;
    else n_ichi = n_points / (m + 1) + (int)(first < (n_points % (m + 1)));

    // no. chiasma from no interference process
    if (p > 0) n_nichi = R::rpois(lambda2);
    else n_nichi = 0;

    if (!obligate_chiasma || n_ichi + n_nichi > 0) break;
  }

  // locations of chiasmata and intermediate points for process w/ interference
  // arma::vec point_locations = arma::sort(arma::randu<arma::vec>(n_points)) * L;
  arma::vec point_locations = arma::randu<arma::vec>(n_points) * L;
  std::sort(point_locations.begin(), point_locations.end());

  // move every (m+1)st point back to front
  arma::uword n_chi = 0;
  for (arma::uword j = first; j < n_points; j += (m + 1), ++n_chi)
    point_locations(n_chi) = point_locations(j);

  // chiasma locations from non-interference process
  const arma::vec nichi_locations = arma::randu<arma::vec>(n_nichi) * L;

  // combine interference and no interference chiasma locations
  arma::vec chi_locations(n_chi + n_nichi);
  // chi_locations(span(0, n_chi - 1)) = point_locations(span(0, n_chi - 1));
  // chi_locations(span(n_chi, n_chi + n_nichi - 1)) = nichi_locations;
  std::copy(point_locations.begin(), point_locations.begin() + n_chi, chi_locations.begin());
  std::copy(nichi_locations.begin(), nichi_locations.end(), chi_locations.begin() + n_chi);
  // chi_locations = arma::sort(chi_locations);
  std::sort(chi_locations.begin(), chi_locations.end());

  // thining by 0.5
  arma::uword n_xo = 0;
  const arma::vec coins = arma::randu<arma::vec>(n_chi + n_nichi);
  for(arma::uword i = 0; i < n_chi + n_nichi; ++i) {
    if (coins[i] < 0.5) {
      chi_locations(n_xo) = chi_locations(i);
      ++n_xo;
    }
  }

  const arma::vec xo_locations(chi_locations.memptr(), n_xo, false, false);
  return xo_locations;
}

// // This function is a fall-back. It is fastest for m = 0.
// //' @export
// // [[Rcpp::export(".sim_crossovers_rcpp")]]
// Rcpp::NumericVector sim_crossovers_rcpp(const double L,
//                                         const int m,
//                                         const double p,
//                                         const bool obligate_chiasma,
//                                         const double Lstar)
// {
//   if(m==0) { // no-interference model is a lot easier
//     int n_xo;
//
//     if(obligate_chiasma) {
//       // rejection sampling to get at least one chiasma
//       while((n_xo = R::rpois(Lstar/50.0)) == 0);
//
//       n_xo = R::rbinom((double)n_xo, 0.5);
//     }
//     else {
//       n_xo = R::rpois(L/100.0);
//     }
//     // SORT HERE!
//     NumericVector tmp = runif(n_xo, 0.0, L);
//     std::sort(tmp.begin(), tmp.end()); // faster
//     return tmp;
//   }
//
//   int n_points, first, n_nichi, n_ichi;
//
//   double lambda1 = Lstar/50.0 * (double)(m+1) * (1.0 - p);
//   double lambda2 = Lstar/50.0 * p;
//
//   while(1) {
//     // chiasma and intermediate points
//     n_points = R::rpois(lambda1);
//
//     // which point is the first chiasma?
//     first = random_int(0, m);
//     if(first > n_points) n_ichi = 0;
//     else n_ichi = n_points/(m+1) + (int)(first < (n_points % (m+1)));
//
//     // no. chiasma from no interference process
//     if(p > 0) n_nichi = R::rpois(lambda2);
//     else n_nichi = 0;
//
//     if(!obligate_chiasma || n_ichi + n_nichi > 0) break;
//   }
//
//   // locations of chiasmata and intermediate points for process w/ interference
//   NumericVector point_locations = runif(n_points, 0.0, L);
//   point_locations.sort();
//
//   // move every (m+1)st point back to front
//   int n_chi=0;
//   for(int j=first; j < n_points; j += (m+1), n_chi++)
//     point_locations[n_chi] = point_locations[j];
//
//   // chiasma locations from non-interference process
//   NumericVector nichi_locations = runif(n_nichi, 0.0, L);
//
//   // combine interference and no interference chiasma locations
//   NumericVector chi_locations(n_chi + n_nichi);
//   std::copy(point_locations.begin(), point_locations.begin()+n_chi, chi_locations.begin());
//   std::copy(nichi_locations.begin(), nichi_locations.end(), chi_locations.begin()+n_chi);
//   chi_locations.sort();
//
//   // thin by 1/2
//   int n_xo=0;
//   for(int i=0; i<n_chi+n_nichi; i++) {
//     if(R::unif_rand() < 0.5) { // flip coin -> chiasma
//       chi_locations[n_xo] = chi_locations[i];
//       n_xo++;
//     }
//   }
//
//   NumericVector xo_locations(n_xo);
//   std::copy(chi_locations.begin(), chi_locations.begin()+n_xo, xo_locations.begin());
//   return xo_locations;
// }
