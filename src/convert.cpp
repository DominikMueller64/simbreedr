// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;
using namespace std;
// #include "convert.h"

// arma::rowvec xo2geno_chromatid_reliable(const List& xodat,
//                                         const arma::vec& map,
//                                         const arma::mat& founder)
// {
//   ivec alleles = xodat["alleles"];
//   vec locations = xodat["locations"];
//   // Minimal speedup possible by using indicies, but more dangerous.
//   // ivec alleles = xodat[0];
//   // vec locations = xodat[1];
//
//   int n_mar = map.n_elem;
//   rowvec out(n_mar);
//
//   int j = 0;
//   for(int i = 0; i < n_mar; i++) {
//     // advance/
//     double pos = map[i];
//     while(locations[j] < pos) {
//       j++;
//     }
//     out[i] = founder(alleles[j], i);
//   }
//
//   return out;
// }

// arma::rowvec xo2geno_chromatid_dev(const List& xodat,
//                                    const arma::vec& map,
//                                    const arma::mat& founder)
// {
//   const arma::ivec alleles = xodat["alleles"];
//   const arma::vec locations = xodat["locations"];
//
//   unsigned int n_mar = map.n_elem;
//   arma::rowvec out(n_mar);
//   arma::uword j = 0;
//   double loc = locations[j];
//   for(arma::uword i = 0; i < n_mar; ++i) {
//     while(loc < map[i]) {
//       loc = locations[++j];
//     }
//     out[i] = founder.at(i, alleles[j]);
//   }
//   return out;
// }


//' @export
//[[Rcpp::export(".xo2geno_chromatid")]]
arma::rowvec xo2geno_chromatid(const List& xodat,
                               const arma::vec& map,
                               const arma::mat& founder)
{
  // const arma::ivec alleles = xodat["alleles"];
  // const arma::vec locations = xodat["locations"];
  const arma::ivec alleles = xodat[0];
  const arma::vec locations = xodat[1];

  arma::rowvec out(map.n_elem);
  arma::uword ix1 = 0;
  arma::uword ix2;
  for (arma::uword j = 0; j < locations.n_elem; ++j) {
    auto it = std::upper_bound(map.begin() + ix1, map.end(), locations[j]);
    ix2 = it - map.begin();
    for (arma::uword i = ix1; i < ix2; ++i) {
      out[i] = founder.at(i, alleles[j]);
    }
    ix1 = ix2;
  }
  return out;
}

// //' @export
// //[[Rcpp::export(".xo2geno_chromatid_vert")]]
// arma::vec xo2geno_chromatid_vert(const List& xodat,
//                                  const arma::vec& map,
//                                  const arma::mat& founder)
// {
//   // const arma::ivec alleles = xodat["alleles"];
//   // const arma::vec locations = xodat["locations"];
//   const arma::ivec alleles = xodat[0];
//   const arma::vec locations = xodat[1];
//
//   arma::vec out(map.n_elem);
//   arma::uword ix1 = 0;
//   arma::uword ix2;
//   for (arma::uword j = 0; j < locations.n_elem; ++j) {
//     auto it = std::upper_bound(map.begin() + ix1, map.end(), locations[j]);
//     ix2 = it - map.begin();
//     for (arma::uword i = ix1; i < ix2; ++i) {
//       out[i] = founder.at(i, alleles[j]);
//     }
//     ix1 = ix2;
//   }
//   return out;
// }



//' @export
// [[Rcpp::export('.xo2geno_chromosome')]]
arma::mat xo2geno_chromosome(const Rcpp::List& xodat,
                              const arma::vec& map,
                              const arma::mat& founder,
                              const bool homozygous)
{
  const unsigned int n_loci = map.size();
  if (homozygous) {
    arma::rowvec tmp = xo2geno_chromatid(xodat, map, founder);
    return arma::join_vert(tmp, tmp);
  } else {
    return arma::join_vert(xo2geno_chromatid(xodat[0], map, founder),
                           xo2geno_chromatid(xodat[1], map, founder));
  }
}

// //' @export
// // [[Rcpp::export('.xo2geno_chromosome_vert')]]
// arma::mat xo2geno_chromosome_vert(const Rcpp::List& xodat,
//                                   const arma::vec& map,
//                                   const arma::mat& founder,
//                                   const bool homozygous)
// {
//   const unsigned int n_loci = map.size();
//   if (homozygous) {
//     arma::rowvec tmp = xo2geno_chromatid_vert(xodat, map, founder);
//     return arma::join_horiz(tmp, tmp);
//   } else {
//     return arma::join_horiz(xo2geno_chromatid_vert(xodat[0], map, founder),
//                             xo2geno_chromatid_vert(xodat[1], map, founder));
//   }
// }



//' @export
// [[Rcpp::export('.xo2geno_individual')]]
arma::mat xo2geno_individual(const Rcpp::List& xodat,
                             const Rcpp::List& map,
                             const Rcpp::List& founder)
{
  const bool homozygous = xodat.attr("homozygous");
  const unsigned int n_chr = map.size();

  arma::uvec map_len(n_chr);
  for (arma::uword i = 0; i < n_chr; ++i)
    map_len[i] = (Rcpp::as<Rcpp::NumericVector>(map[i])).size();

  const unsigned int n_total = arma::sum(map_len);
  arma::mat matr(2, n_total);
  arma::uword fi = 0;
  arma::uword la = 0;
  for (arma::uword i = 0; i < n_chr; ++i) {
    la = fi + map_len[i] - 1;
    matr.cols(fi, la) = xo2geno_chromosome(xodat[i], map[i], founder[i], homozygous);
    fi = la + 1;
  }
  return matr;
}

// //' @export
// // [[Rcpp::export('.xo2geno_individual_vert')]]
// arma::mat xo2geno_individual_vert(const Rcpp::List& xodat,
//                                   const Rcpp::List& map,
//                                   const Rcpp::List& founder)
// {
//   const bool homozygous = xodat.attr("homozygous");
//   const unsigned int n_chr = map.size();
//
//   arma::uvec map_len(n_chr);
//   for (arma::uword i = 0; i < n_chr; ++i)
//     map_len[i] = (Rcpp::as<Rcpp::NumericVector>(map[i])).size();
//
//   const unsigned int n_total = arma::sum(map_len);
//   arma::mat matr(n_total, 2);
//   arma::uword fi = 0;
//   arma::uword la = 0;
//   for (arma::uword i = 0; i < n_chr; ++i) {
//     la = fi + map_len[i] - 1;
//     matr.rows(fi, la) = xo2geno_chromosome_vert(xodat[i], map[i], founder[i], homozygous);
//     fi = la + 1;
//   }
//   return matr;
// }


//' @export
// [[Rcpp::export('.xo2geno_gamete')]]
arma::mat xo2geno_gamete(const Rcpp::List& xodat,
                         const Rcpp::List& map,
                         const Rcpp::List& founder)
{
  const unsigned int n_chr = map.size();

  arma::uvec map_len(n_chr);
  for (arma::uword i = 0; i < n_chr; ++i)
    map_len[i] = (Rcpp::as<Rcpp::NumericVector>(map[i])).size();

  const unsigned int n_total = arma::sum(map_len);
  arma::rowvec matr(n_total);

  arma::uword fi = 0;
  arma::uword la = 0;
  for (arma::uword i = 0; i < n_chr; ++i) {
    la = fi + map_len[i] - 1;
    matr.subvec(fi, la) = xo2geno_chromatid(xodat[i], map[i], founder[i]);
    fi = la + 1;
  }
  return matr;
}

//' @export
// [[Rcpp::export(".xo2geno_population")]]
arma::mat xo2geno_population(const Rcpp::List& xodat,
                             const Rcpp::List& map,
                             const Rcpp::List& founder)
{
  const unsigned int n_ind = xodat.size();
  const unsigned int n_chr = map.size();
  arma::uvec map_len(n_chr);
  for (arma::uword i = 0; i < n_chr; ++i)
    map_len[i] = (Rcpp::as<Rcpp::NumericVector>(map[i])).size();

  const unsigned int n_total = arma::sum(map_len);
  arma::mat matr(2 * n_ind, n_total);

  arma::uword row = 0;
  for (arma::uword i = 0; i < n_ind; ++i, row += 2) {
    matr.rows(row, row + 1) = xo2geno_individual(xodat[i], map, founder);
  }
  return matr;
}


//// [[Rcpp::export('.xo2geno_mat')]]
// mat xo2geno_Matrix(const List xodat, const List map, const List founder) {
//   int n_ind = xodat.size();
//   int n_chrom = map.size();
//   sapply(xodat, Rcpp::length);
//
// //   int n_ind = ind.size();
// //   NumericMatrix out(2 * n_ind, map.size());
// //   int j = 0;
// //   for(int i = 0; i < n_ind; i++) {
// //     List ind_ = ind[i];
// //     out(j++,_) = xo2geno_chromatid(ind_[0], map, founder);
// //     out(j++,_) = xo2geno_chromatid(ind_[1], map, founder);
// //   }
// //   return out;
//    return mat(2,2);
// }

// //' @export
// //[[Rcpp::export]]
// List xo2geno(const List xodat, const List map,
//              const List founder) {
//   int n = xodat.size();
//   List out(n);
//   for(int i = 0; i < n; i++) {
//     out[i] = xo2geno_individual_List(xodat[i], map, founder);
//   }
//   return out;
// }
//
//
// //[[Rcpp::export]]
// List xo2geno_homozygous(const List xodat, const List map,
//                         const List founder, const bool aggregate) {
//   int n = xodat.size();
//   List out(n);
//   for(int i = 0; i < n; i++) {
//     List ind_ = xodat[i];
//     rowvec gam = xo2geno_chromatid(ind_["gam"], map[i], founder[i]);
//     if (aggregate) {
//       out[i] = List::create(Rcpp::Named("geno") = 2.0 * gam);
//     } else {
//       out[i] = List::create(Rcpp::Named("mat") = gam,
//                             Rcpp::Named("pat") = gam);
//     }
//   }
//   out.attr("homozygous") = true;
//   return out;
// }
//
//
//' @export
//[[Rcpp::export]]
Rcpp::List xo2geno(const Rcpp::List& xodat,
                   const Rcpp::List& map,
                   const Rcpp::List& founder) {

  const bool homozygous = xodat.attr("homozygous");
  const unsigned int n = xodat.size();
  Rcpp::List ret(n);
  if (homozygous) {
    for(int i = 0; i < n; i++) {
      ret[i] = xo2geno_chromatid(xodat[i], map[i], founder[i]);
    }
  }
  else {
    for(int i = 0; i < n; i++) {
      ret[i] = Rcpp::List::create(Rcpp::Named("mat") = xo2geno_chromatid(Rcpp::as<Rcpp::List>(xodat[i])["mat"], map[i], founder[i]),
                                  Rcpp::Named("pat") = xo2geno_chromatid(Rcpp::as<Rcpp::List>(xodat[i])["pat"], map[i], founder[i]));
    }
  }
  return ret;
}

// // Second lowest level.
// //' @export
// //[[Rcpp::export]]
// List xo2geno_homolog_List(const List& xodat,
//                            const arma::vec& map,
//                            const arma::mat& founder) {
//
//   bool homozygous = xodat.attr("homozygous");
//   rowvec mat, pat;
//   if (homozygous) {
//     mat = xo2geno_chromatid(xodat["par"], map, founder);
//     pat = mat;
//   } else {
//     mat = xo2geno_chromatid(xodat["mat"], map, founder);
//     pat = xo2geno_chromatid(xodat["pat"], map, founder);
//   }
//   return List::create(Rcpp::Named("mat") = mat,
//                       Rcpp::Named("pat") = pat);
// }

// //' @export
// //[[Rcpp::export]]
// arma::mat xo2geno_homolog_Matrix(const List xodat,
//                            const arma::vec& map,
//                            const arma::mat& founder) {
//
//   bool homozygous = xodat.attr("homozygous");
//   mat ret(2, map.n_elem);
//
//   rowvec mat, pat;
//   if (homozygous) {
//     mat = xo2geno_chromatid(xodat["par"], map, founder);
//     pat = mat;
//   } else {
//     mat = xo2geno_chromatid(xodat["mat"], map, founder);
//     pat = xo2geno_chromatid(xodat["pat"], map, founder);
//   }
//   ret.row(0) = mat;
//   ret.row(1) = pat;
//   return ret;
// }

