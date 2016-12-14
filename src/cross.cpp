// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;
#include "sim_meiosis.h"

//' @export
//[[Rcpp::export(".gamete")]]
Rcpp::List gamete(const List& parent,
                  const List& params)
{
  const unsigned int n_chrom = parent.size();
  Rcpp::List ret(n_chrom);
  const bool homozygous = parent.attr("homozygous");
  const arma::vec L = params["L"];
  const arma::vec Lstar = params["Lstar"];
  const unsigned int m = params["m"];
  const double p = params["p"];
  const bool obligate_chiasma = params["obligate_chiasma"];
  for (unsigned int i = 0; i < n_chrom; ++i) {
    Rcpp::List gam;
    if (homozygous) {
      ret[i] = parent[i];
    } else {
      ret[i] = sim_meiosis(parent[i], L[i], m, p, obligate_chiasma, Lstar[i]);
    }
  }
  ret.attr("names") = parent.attr("names");
  ret.attr("individual") = false;
  return ret;
}

//' @export
//[[Rcpp::export(".cross")]]
Rcpp::List cross(const List& mother,
                 const List& father,
                 const List& params)
{
  const unsigned int n_chrom = mother.size();
  Rcpp::List ret(n_chrom);
  const bool mother_homo = mother.attr("homozygous");
  const bool father_homo = father.attr("homozygous");
  const arma::vec L = params["L"];
  const arma::vec Lstar = params["Lstar"];
  const unsigned int m = params["m"];
  const double p = params["p"];
  bool obligate_chiasma = params["obligate_chiasma"];
  for (unsigned int i = 0; i < n_chrom; ++i) {
    Rcpp::List mat, pat;
    if (mother_homo) {
      mat = mother[i];
    } else {
      mat = sim_meiosis(mother[i], L[i], m, p, obligate_chiasma, Lstar[i]);
    }
    if (father_homo) {
      pat = father[i];
    } else {
      pat = sim_meiosis(father[i], L[i], m, p, obligate_chiasma, Lstar[i]);
    }

    ret[i] = Rcpp::List::create(Rcpp::Named("mat") = mat,
                                Rcpp::Named("pat") = pat);
  }
  ret.attr("names") = mother.attr("names");
  ret.attr("homozygous") = false;
  ret.attr("individual") = true;
  return ret;
}


// //' @export
// //[[Rcpp::export(".cross_exp")]]
// Rcpp::List cross_exp(const List& mother,
//                  const List& father,
//                  const NumericVector& L,
//                  const int m,
//                  const double p,
//                  const bool obligate_chiasma,
//                  const NumericVector& Lstar)
// {
//   int n_chrom = mother.size();
//   Rcpp::List ret(n_chrom);
//   bool mother_homo = mother.attr("homozygous");
//   bool father_homo = father.attr("homozygous");
//   for (int i = 0; i < n_chrom; ++i) {
//     List mat, pat;
//     if (mother_homo) {
//       mat = as<Rcpp::List>(mother[i])[0];
//     } else {
//       mat = sim_meiosis(mother[i], L[i], m, p, obligate_chiasma, Lstar[i]);
//     }
//     if (father_homo) {
//       pat = as<Rcpp::List>(father[i])[0];
//     } else {
//       pat = sim_meiosis(father[i], L[i], m, p, obligate_chiasma, Lstar[i]);
//     }
//
//     List tmp = List::create(Rcpp::Named("mat") = mat,
//                             Rcpp::Named("pat") = pat);
//     // tmp.attr("homozygous") = false;
//     ret[i] = tmp;
//   }
//   ret.attr("names") = mother.attr("names");
//   ret.attr("homozygous") = false;
//   return ret;
// }

// List dh(const List parent,
//         const NumericVector& L,
//         const int m,
//         const double p,
//         const bool obligate_chiasma,
//         const NumericVector Lstar) {
//
//   int n_chrom = parent.size();
//   bool homozygous = parent.attr("homozygous");
//   List ret(n_chrom);
//   for (int i = 0; i < n_chrom; i++) {
//     List gam;
//     if (homozygous) {
//       gam = parent[i];
//     } else {
//       List gam = sim_meiosis(parent[i], L[i], m, p, obligate_chiasma, Lstar[i]);
//     }
//     List tmp = List::create(Rcpp::Named("par") = gam);
//     tmp.attr("homozygous") = true;
//     ret[i] = tmp;
//   }
//   ret.attr("names") = parent.attr("names");
//   ret.attr("homozygous") = true;
//   return ret;
// }

// //[[Rcpp::export]]
// List self(const List parent,
//           const NumericVector L, const int m,
//           const double p, const bool obligate_chiasma,
//           const NumericVector Lstar) {
//   return cross(parent, parent, L, m, p, obligate_chiasma, Lstar);
// }

// //[[Rcpp::export]]
// void test(const List params)
// {
//   int n_chrom = params.size();
//   double L, p, Lstar;
//   int m;
//   bool obligate_chiasma;
//   for (int i = 0; i < n_chrom; i++) {
//     L = Rcpp::as<double>( Rcpp::as<Rcpp::List>(params[i])["L"] );
//     m = Rcpp::as<double>( Rcpp::as<Rcpp::List>(params[i])["m"] );
//     p = Rcpp::as<double>( Rcpp::as<Rcpp::List>(params[i])["p"] );
//     obligate_chiasma = Rcpp::as<double>( Rcpp::as<Rcpp::List>(params[i])["obligate_chiasma"] );
//     Lstar = Rcpp::as<double>( Rcpp::as<Rcpp::List>(params[i])["Lstar"] );
//   }
// }
