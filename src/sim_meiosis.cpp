// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

// #include "random.h"
#include "sim_crossovers.h"

//' @export
// [[Rcpp::export]]
List sim_meiosis_xo(const List parent,
                 arma::vec& Xlocations) // weiredly, arma:: is necessary for the package to compile.
{
  const double tol = 1e-12; // for comparison of chr lengths in parents

  List mat = parent[0];
  List pat = parent[1];

  ivec matalle = mat[0];
  vec matloc  = mat[1];

  ivec patalle = pat[0];
  vec patloc  = pat[1];

  int matloc_size = matloc.n_elem;
  int patloc_size = patloc.n_elem;
  double L = matloc[matloc_size - 1];

  // simulate crossover locations; add -1 to the beginning
  vec tmp = Xlocations;
  int product_size = tmp.n_elem + 1;
  vec product(product_size);
  product[0] = -1.0;
  std::copy(tmp.begin(), tmp.end(), product.begin() + 1);

  // int cur_allele = random_int(0, 1); // first allele (0 or 1)
  int cur_allele = (int) R::runif(0, 2);

  int biggest_length = product_size + matloc_size + patloc_size;
  vec loc(biggest_length);
  ivec alle(biggest_length);

  int curpos = 0;
    if(product_size == 1) {
        if(cur_allele == 0) return mat;
        else return pat;
    }
    else {
        int i, j;
        int jmatstart = 0;
        int jpatstart = 0;
        for(i = 1; i < product_size; i++) {

          if(cur_allele == 0) { // mat chr
            for(j = jmatstart; j < matloc.size(); j++) {
              if(matloc[j] >= product[i-1] && matloc[j] < product[i]) {
                loc[curpos] = matloc[j];
                alle[curpos] = matalle[j];
                curpos++;
                jmatstart++;
              }
              else if(matloc[j] > product[i]) break;
            }
            loc[curpos] = product[i];
            alle[curpos] = matalle[j];
            curpos++;
          }
          else { // pat chr
            for(j = jpatstart; j < patloc.size(); j++) {
              if(patloc[j] >= product[i-1] && patloc[j] < product[i]) {
                loc[curpos] = patloc[j];
                alle[curpos] = patalle[j];
                curpos++;
                jpatstart++;
              }
              else if(patloc[j] > product[i]) break;
            }
            loc[curpos] = product[i];
            alle[curpos] = patalle[j];
            curpos++;
          }

          cur_allele = 1 - cur_allele;
        }

        double lastxo = product[product_size - 1];

        if(cur_allele == 0) { // mat chr
            for(j = 0; j < matloc_size; j++) {
                if(matloc[j] > lastxo) {
                    loc[curpos] = matloc[j];
                    alle[curpos] = matalle[j];
                    curpos++;
                }
            }
        }
        else { // pat chr
            for(j = 0; j < patloc_size; j++) {
                if(patloc[j] > lastxo) {
                    loc[curpos] = patloc[j];
                    alle[curpos] = patalle[j];
                    curpos++;
                }
            }
        }
    }

    if(curpos > 1) { // clean up repeated alleles

        vec loc_clean(curpos);
        ivec alle_clean(curpos);

        loc_clean[0] = loc[0];
        alle_clean[0] = alle[0];
        int lastpos = 0;

        for(int i = 1; i < curpos; i++) {
            if(alle_clean[lastpos] == alle[i]) {
                loc_clean[lastpos] = loc[i];
            }
            else {
                lastpos++;
                loc_clean[lastpos] = loc[i];
                alle_clean[lastpos] = alle[i];
            }
        }
        curpos = lastpos + 1;
        loc = loc_clean;
        alle = alle_clean;
    }

    // copy over to short vectors
    NumericVector loc_result(loc.begin(), loc.begin() + curpos);
    IntegerVector alle_result(alle.begin(), alle.begin() + curpos);
//     NumericVector loc_result(curpos);
//     IntegerVector alle_result(curpos);
//     std::copy(loc.begin(), loc.begin() + curpos, loc_result.begin());
//     std::copy(alle.begin(), alle.begin() + curpos, alle_result.begin());

    return List::create(Named("alleles") = alle_result,
                        Named("locations") = loc_result);

}


//' @export
// [[Rcpp::export]]
List sim_meiosis_xo_test(const List parent,
                 arma::vec& Xlocations) // weiredly, arma:: is necessary for the package to compile.
{
  const double tol = 1e-12; // for comparison of chr lengths in parents

  List mat = parent[0];
  List pat = parent[1];

  ivec matalle = mat[0];
  vec matloc  = mat[1];

  ivec patalle = pat[0];
  vec patloc  = pat[1];

  int matloc_size = matloc.n_elem;
  int patloc_size = patloc.n_elem;
  double L = matloc[matloc_size - 1];

  // simulate crossover locations; add -1 to the beginning
  vec tmp = Xlocations;
  int product_size = tmp.n_elem + 1;
  vec product(product_size);
  product[0] = -1.0;
  std::copy(tmp.begin(), tmp.end(), product.begin() + 1);

  // int cur_allele = random_int(0, 1); // first allele (0 or 1)
  int cur_allele = (int) R::runif(0, 2);

  int biggest_length = product_size + matloc_size + patloc_size;
  vec loc(biggest_length);
  ivec alle(biggest_length);

  int curpos = 0;
    if(product_size == 1) {
        if(cur_allele == 0) return mat;
        else return pat;
    }
    else {
        int jpatstart = 0;
        int jmat = 0;
        int jpat = 0;
        double prod0, prod1;
        for(int i = 1; i < product_size; i++) {
          prod0 = product[i - 1];
          prod1 = product[i];
          if(cur_allele == 0) { // mat chr
            while(matloc[jmat] < prod0) {
              ++jmat;
            }
            while(matloc[jmat] < prod1 && jmat < matloc.size()) {
              loc[curpos] = matloc[jmat];
              alle[curpos] = matalle[jmat];
              jmat++;
              curpos++;
            }
            loc[curpos] = prod1;
            alle[curpos] = matalle[jmat];
            curpos++;
          }
          else { // pat chr
            while(patloc[jpat] < prod0) {
              ++jpat;
            }
            while(patloc[jpat] < prod1 && jpat < patloc.size()) {
              loc[curpos] = patloc[jpat];
              alle[curpos] = patalle[jpat];
              jpat++;
              curpos++;
            }
            loc[curpos] = prod1;
            alle[curpos] = patalle[jpat];
            curpos++;
          }
          cur_allele = 1 - cur_allele;
        }

        double lastxo = product[product_size - 1];

        if(cur_allele == 0) { // mat chr
            for(int j = 0; j < matloc_size; j++) {
                if(matloc[j] > lastxo) {
                    loc[curpos] = matloc[j];
                    alle[curpos] = matalle[j];
                    curpos++;
                }
            }
        }
        else { // pat chr
            for(int j = 0; j < patloc_size; j++) {
                if(patloc[j] > lastxo) {
                    loc[curpos] = patloc[j];
                    alle[curpos] = patalle[j];
                    curpos++;
                }
            }
        }
    }

    if(curpos > 1) { // clean up repeated alleles

      // vec loc_clean(curpos);
      // ivec alle_clean(curpos);

      // loc_clean[0] = loc[0];
      // alle_clean[0] = alle[0];
      int lastpos = 0;

      for(int i = 1; i < curpos; i++) {
        if(alle[lastpos] == alle[i]) {
          loc[lastpos] = loc[i];
        }
        else {
          lastpos++;
          loc[lastpos] = loc[i];
          alle[lastpos] = alle[i];
        }
      }
      curpos = lastpos + 1;
      // loc = loc_clean;
      // alle = alle_clean;
    }

    // copy over to short vectors
    NumericVector loc_result(loc.begin(), loc.begin() + curpos);
    IntegerVector alle_result(alle.begin(), alle.begin() + curpos);
    //     NumericVector loc_result(curpos);
    //     IntegerVector alle_result(curpos);
    //     std::copy(loc.begin(), loc.begin() + curpos, loc_result.begin());
    //     std::copy(alle.begin(), alle.begin() + curpos, alle_result.begin());

    return List::create(Named("alleles") = alle_result,
                        Named("locations") = loc_result);

}




// This function is nearly optimized and should work reliably.
//' @export
// [[Rcpp::export]]
Rcpp::List sim_meiosis_reliable(const Rcpp::List& parent,
                       const double L,
                       const int m,
                       const double p,
                       const bool obligate_chiasma,
                       const double Lstar)
{
  const Rcpp::List mat = parent["mat"];
  const Rcpp::List pat = parent["pat"];

  const arma::ivec matalle = mat["alleles"];
  const arma::vec matloc  = mat["locations"];

  const arma::ivec patalle = pat["alleles"];
  const arma::vec patloc  = pat["locations"];

  const int matloc_size = matloc.n_elem;
  const int patloc_size = patloc.n_elem;

  // simulate crossover locations; add -1 to the beginning
  arma::vec tmp = sim_crossovers(L, m, p, obligate_chiasma, Lstar);
  const int product_size = tmp.n_elem + 1;
  vec product(product_size);
  product[0] = -1.0;
  std::copy(tmp.begin(), tmp.end(), product.begin() + 1);

  int cur_allele = (int) R::runif(0, 2);
  int biggest_length = product_size + matloc_size + patloc_size;
  arma::vec loc(biggest_length);
  arma::ivec alle(biggest_length);

  int curpos = 0;
  if (product_size == 1) {
    if (cur_allele == 0) return mat;
    else return pat;
  }
  else {
     int i, j;
    int jmatstart = 0;
    int jpatstart = 0;
    for(i = 1; i < product_size; i++) {

      if(cur_allele == 0) { // mat chr
        for(j = jmatstart; j < matloc.size(); j++) {
          if(matloc[j] >= product[i-1] && matloc[j] < product[i]) {
            loc[curpos] = matloc[j];
            alle[curpos] = matalle[j];
            curpos++;
            jmatstart++;
          }
          else if(matloc[j] > product[i]) break;
        }
        loc[curpos] = product[i];
        alle[curpos] = matalle[j];
        curpos++;
      }
      else { // pat chr
        for(j = jpatstart; j < patloc.size(); j++) {
          if(patloc[j] >= product[i-1] && patloc[j] < product[i]) {
            loc[curpos] = patloc[j];
            alle[curpos] = patalle[j];
            curpos++;
            jpatstart++;
          }
          else if(patloc[j] > product[i]) break;
        }
        loc[curpos] = product[i];
        alle[curpos] = patalle[j];
        curpos++;
      }

      cur_allele = 1 - cur_allele;
    }

    const double lastxo = product[product_size - 1];

    if(cur_allele == 0) { // mat chr
      for(j = 0; j < matloc_size; j++) {
        if(matloc[j] > lastxo) {
          loc[curpos] = matloc[j];
          alle[curpos] = matalle[j];
          curpos++;
        }
      }
    }
    else{ // pat chr
      for(j = 0; j < patloc_size; j++) {
        if(patloc[j] > lastxo) {
          loc[curpos] = patloc[j];
          alle[curpos] = patalle[j];
          curpos++;
        }
      }
    }
  }

  if(curpos > 1) { // clean up repeated alleles

    arma::vec loc_clean(curpos);
    arma::ivec alle_clean(curpos);

    loc_clean[0] = loc[0];
    alle_clean[0] = alle[0];
    int lastpos = 0;

    for(int i = 1; i < curpos; i++) {
      if(alle_clean[lastpos] == alle[i]) {
        loc_clean[lastpos] = loc[i];
      }
      else {
        lastpos++;
        loc_clean[lastpos] = loc[i];
        alle_clean[lastpos] = alle[i];
      }
    }
    curpos = lastpos + 1;
    loc = loc_clean;
    alle = alle_clean;
  }

  // copy over to short vectors
  Rcpp::NumericVector loc_result(loc.begin(), loc.begin() + curpos);
  Rcpp::IntegerVector alle_result(alle.begin(), alle.begin() + curpos);
  return Rcpp::List::create(Rcpp::Named("alleles") = alle_result,
                            Rcpp::Named("locations") = loc_result);
}

// // This is the lastest function without the while-solution
// //' @export
// // [[Rcpp::export]]
// Rcpp::List sim_meiosis(const Rcpp::List& parent,
//                        const double L,
//                        const int m,
//                        const double p,
//                        const bool obligate_chiasma,
//                        const double Lstar)
// {
//   // Simulate crossover locations.
//   const arma::vec tmp = sim_crossovers(L, m, p, obligate_chiasma, Lstar);
//   // Draw the parent where to start from.
//   int cur_allele = (int) R::runif(0, 2);
//   // Performance gain by using index-based subsetting is substantial.
//   const Rcpp::List mat = parent[0];
//   const Rcpp::List pat = parent[1];
//   // If no crossover occurs, early exit.
//   if (tmp.size() == 0) {
//     if (cur_allele == 0) return mat;
//     else return pat;
//   }
//
//   const arma::ivec matalle = mat[0];
//   const arma::vec matloc  = mat[1];
//   const arma::ivec patalle = pat[0];
//   const arma::vec patloc  = pat[1];
//
//   const int matloc_size = matloc.size();
//   const int patloc_size = patloc.size();
//
//   const unsigned int product_size = tmp.size() + 1;
//   arma::vec product(product_size);
//   product[0] = -1.0;
//   std::copy(tmp.begin(), tmp.end(), product.begin() + 1);
//
//   const unsigned int biggest_length = product_size + matloc_size + patloc_size;
//   arma::vec loc(biggest_length);
//   arma::ivec alle(biggest_length);
//
//   int curpos = 0;
//   arma::uword i, j;
//   arma::uword jmatstart = 0;
//   arma::uword jpatstart = 0;
//   for(i = 1; i < product_size; ++i) {
//
//     if(cur_allele == 0) { // mat
//       for(j = jmatstart; j < matloc_size; ++j) {
//         if(matloc[j] >= product[i - 1] && matloc[j] < product[i]) {
//           loc[curpos] = matloc[j];
//           alle[curpos] = matalle[j];
//           ++curpos;
//           ++jmatstart;
//         }
//         else if(matloc[j] > product[i]) break;
//       }
//       loc[curpos] = product[i];
//       alle[curpos] = matalle[j];
//       ++curpos;
//     }
//     else { // pat
//       for(j = jpatstart; j < patloc_size; ++j) {
//         if(patloc[j] >= product[i - 1] && patloc[j] < product[i]) {
//           loc[curpos] = patloc[j];
//           alle[curpos] = patalle[j];
//           ++curpos;
//           ++jpatstart;
//         }
//         else if(patloc[j] > product[i]) break;
//       }
//       loc[curpos] = product[i];
//       alle[curpos] = patalle[j];
//       ++curpos;
//     }
//
//     cur_allele = 1 - cur_allele;
//   }
//
//   const double lastxo = product[product_size - 1];
//
//   if(cur_allele == 0) { // mat chr
//     for(j = 0; j < matloc_size; ++j) {
//       if(matloc[j] > lastxo) {
//         loc[curpos] = matloc[j];
//         alle[curpos] = matalle[j];
//         ++curpos;
//       }
//     }
//   }
//   else{ // pat chr
//     for(j = 0; j < patloc_size; ++j) {
//       if(patloc[j] > lastxo) {
//         loc[curpos] = patloc[j];
//         alle[curpos] = patalle[j];
//         ++curpos;
//       }
//     }
//   }
//
//   if(curpos > 1) { // clean up repeated alleles
//
//     arma::vec loc_clean(curpos);
//     arma::ivec alle_clean(curpos);
//
//     loc_clean[0] = loc[0];
//     alle_clean[0] = alle[0];
//     arma::uword lastpos = 0;
//
//     for(i = 1; i < curpos; ++i) {
//       if(alle_clean[lastpos] == alle[i]) {
//         loc_clean[lastpos] = loc[i];
//       }
//       else {
//         ++lastpos;
//         loc_clean[lastpos] = loc[i];
//         alle_clean[lastpos] = alle[i];
//       }
//     }
//     curpos = lastpos + 1;
//     loc = loc_clean;
//     alle = alle_clean;
//   }
//
//   // copy over to short vectors
//   Rcpp::NumericVector loc_result(loc.begin(), loc.begin() + curpos);
//   Rcpp::IntegerVector alle_result(alle.begin(), alle.begin() + curpos);
//   // arma::vec loc_result(loc.memptr(), curpos, false, false);
//   // arma::ivec alle_result(alle.memptr(), curpos, false, false);
//   return Rcpp::List::create(Rcpp::Named("alleles") = alle_result,
//                             Rcpp::Named("locations") = loc_result);
// }

// This seems to be the best algorithm, combining advantages. Should
// get standard.
//' @export
// [[Rcpp::export(".sim_meiosis")]]
Rcpp::List sim_meiosis(const Rcpp::List& parent,
                       const double L,
                       const int m,
                       const double p,
                       const bool obligate_chiasma,
                       const double Lstar)
{
  // Simulate crossover locations.
  const arma::vec tmp = sim_crossovers(L, m, p, obligate_chiasma, Lstar);
  // Draw the parent where to start from.
  int cur_allele = (int) R::runif(0, 2);
  // Performance gain by using index-based subsetting is substantial.
  const Rcpp::List mat = parent[0];
  const Rcpp::List pat = parent[1];
  // If no crossover occurs, early exit.
  if (tmp.size() == 0) {
    if (cur_allele == 0) return mat;
    else return pat;
  }

  const arma::ivec matalle = mat[0];
  const arma::vec matloc  = mat[1];
  const arma::ivec patalle = pat[0];
  const arma::vec patloc  = pat[1];

  const int matloc_size = matloc.size();
  const int patloc_size = patloc.size();

  const unsigned int product_size = tmp.size() + 1;
  arma::vec product(product_size);
  product[0] = -1.0;
  std::copy(tmp.begin(), tmp.end(), product.begin() + 1);

  const unsigned int biggest_length = product_size + matloc_size + patloc_size;
  arma::vec loc(biggest_length);
  arma::ivec alle(biggest_length);

  int curpos = 0;
  arma::uword i, j;
  arma::uword jmat = 0;
  arma::uword jpat = 0;
  double prod0, prod1;
  for(i = 1; i < product_size; ++i) {
    prod0 = product[i - 1];
    prod1 = product[i];
    if(cur_allele == 0) { // mat chr
      while(matloc[jmat] < prod0) {
        ++jmat;
      }
      while(matloc[jmat] < prod1 && jmat < matloc.size()) {
        loc[curpos] = matloc[jmat];
        alle[curpos] = matalle[jmat];
        ++jmat;
        ++curpos;
      }
      loc[curpos] = prod1;
      alle[curpos] = matalle[jmat];
      ++curpos;
    }
    else { // pat chr
      while(patloc[jpat] < prod0) {
        ++jpat;
      }
      while(patloc[jpat] < prod1 && jpat < patloc.size()) {
        loc[curpos] = patloc[jpat];
        alle[curpos] = patalle[jpat];
        ++jpat;
        ++curpos;
      }
      loc[curpos] = prod1;
      alle[curpos] = patalle[jpat];
      ++curpos;
    }
    cur_allele = 1 - cur_allele;
  }

  const double lastxo = product[product_size - 1];

  if(cur_allele == 0) { // mat chr
    for(j = 0; j < matloc_size; ++j) {
      if(matloc[j] > lastxo) {
        loc[curpos] = matloc[j];
        alle[curpos] = matalle[j];
        ++curpos;
      }
    }
  }
  else{ // pat chr
    for(j = 0; j < patloc_size; ++j) {
      if(patloc[j] > lastxo) {
        loc[curpos] = patloc[j];
        alle[curpos] = patalle[j];
        ++curpos;
      }
    }
  }

  if(curpos > 1) { // clean up repeated alleles
    arma::uword lastpos = 0;
    for(i = 1; i < curpos; ++i) {
      if(alle[lastpos] == alle[i]) {
        loc[lastpos] = loc[i];
      }
      else {
        ++lastpos;
        loc[lastpos] = loc[i];
        alle[lastpos] = alle[i];
      }
    }
    curpos = lastpos + 1;
  }
  // copy over to short vectors
  Rcpp::NumericVector loc_result(loc.begin(), loc.begin() + curpos);
  Rcpp::IntegerVector alle_result(alle.begin(), alle.begin() + curpos);
  return Rcpp::List::create(Rcpp::Named("alleles") = alle_result,
                            Rcpp::Named("locations") = loc_result);
}




// // This was a nice idea of optimization, but its acutally slower...
// //' @export
// // [[Rcpp::export]]
// Rcpp::List sim_meiosis_WTF(const Rcpp::List& parent,
//                                    arma::vec& Xlocations)
// {
//   // Simulate crossover locations.
//   const arma::vec tmp = Xlocations;
//   // const arma::vec tmp = sim_crossovers(L, m, p, obligate_chiasma, Lstar);
//   // Draw the parent where to start from.
//   int cur_allele = (int) R::runif(0, 2);
//   // Performance gain by using index-based subsetting is substantial.
//   const Rcpp::List mat = parent[0];
//   const Rcpp::List pat = parent[1];
//   // If no crossover occurs, early exit.
//   if (tmp.size() == 0) {
//     if (cur_allele == 0) return mat;
//     else return pat;
//   }
//
//   const arma::ivec matalle = mat[0];
//   const arma::vec matloc  = mat[1];
//   const arma::ivec patalle = pat[0];
//   const arma::vec patloc  = pat[1];
//
//   const int matloc_size = matloc.size();
//   const int patloc_size = patloc.size();
//
//   const unsigned int product_size = tmp.size() + 1;
//   arma::vec product(product_size);
//   product[0] = -1.0;
//   std::copy(tmp.begin(), tmp.end(), product.begin() + 1);
//
//   const unsigned int biggest_length = product_size + matloc_size + patloc_size;
//   arma::vec loc(biggest_length);
//   arma::ivec alle(biggest_length);
//
//   int curpos = 0;
//   arma::uword i, j;
//   arma::uword jmatstart = 0;
//   arma::uword jpatstart = 0;
//   int last_allele = -100;
//   for(i = 1; i < product_size; ++i) {
//     if(cur_allele == 0) { // mat
//       for(j = jmatstart; j < matloc_size; ++j) {
//         if(matloc[j] >= product[i - 1] && matloc[j] < product[i]) {
//           if (last_allele == matalle[j]) {
//             loc[curpos - 1] = matloc[j];
//             // alle[curpos - 1] = matalle[j];
//           } else {
//             loc[curpos] = matloc[j];
//             alle[curpos] = matalle[j];
//             last_allele = matalle[j];
//             ++curpos;
//           }
//           ++jmatstart;
//         }
//         else if(matloc[j] > product[i]) {
//           break;
//         }
//       }
//       if (last_allele == matalle[j]) {
//         loc[curpos - 1] = product[i];
//         // alle[curpos - 1] = matalle[j];
//       } else {
//         loc[curpos] = product[i];
//         alle[curpos] = matalle[j];
//         last_allele = matalle[j];
//         ++curpos;
//       }
//     }
//     else { // pat
//       for(j = jpatstart; j < patloc_size; ++j) {
//         if(patloc[j] >= product[i - 1] && patloc[j] < product[i]) {
//           if (last_allele == patalle[j]) {
//             loc[curpos - 1] = patloc[j];
//             // alle[curpos - 1] = patalle[j];
//           } else {
//             loc[curpos] = patloc[j];
//             alle[curpos] = patalle[j];
//             last_allele = patalle[j];
//             ++curpos;
//           }
//           ++jpatstart;
//         }
//         else if(patloc[j] > product[i]) {
//           break;
//         }
//       }
//       if (last_allele == patalle[j]) {
//         loc[curpos - 1] = product[i];
//         // alle[curpos - 1] = patalle[j];
//       } else {
//         loc[curpos] = product[i];
//         alle[curpos] = patalle[j];
//         last_allele = patalle[j];
//         ++curpos;
//       }
//     }
//     cur_allele = 1 - cur_allele;
//   }
//
//   const double lastxo = product[product_size - 1];
//
//   if(cur_allele == 0) { // mat chr
//     for(j = 0; j < matloc_size; ++j) {
//       if(matloc[j] > lastxo) {
//         if (last_allele == matalle[j]) {
//           loc[curpos - 1] = matloc[j];
//           // alle[curpos - 1] = matalle[j];
//         } else {
//           loc[curpos] = matloc[j];
//           alle[curpos] = matalle[j];
//           last_allele = matalle[j];
//           ++curpos;
//         }
//       }
//     }
//   }
//   else{ // pat chr
//     for(j = 0; j < patloc_size; ++j) {
//       if(patloc[j] > lastxo) {
//         if (last_allele == patalle[j]) {
//           loc[curpos - 1] = patloc[j];
//           // alle[curpos - 1] = patalle[j];
//         } else {
//           loc[curpos] = patloc[j];
//           alle[curpos] = patalle[j];
//           last_allele = patalle[j];
//           ++curpos;
//         }
//       }
//     }
//   }
//
//   // copy over to short vectors
//   Rcpp::NumericVector loc_result(loc.begin(), loc.begin() + curpos);
//   Rcpp::IntegerVector alle_result(alle.begin(), alle.begin() + curpos);
//   // arma::vec loc_result(loc.memptr(), curpos, false, false);
//   // arma::ivec alle_result(alle.memptr(), curpos, false, false);
//   return Rcpp::List::create(Rcpp::Named("alleles") = alle_result,
//                             Rcpp::Named("locations") = loc_result);
// }


// // [[Rcpp::export("sim_meiosis2")]]
// List sim_meiosis2(const List parent, const NumericVector Xlocations)
// {
//   const double tol=1e-12; // for comparison of chr lengths in parents
//
//   List mat, pat;
//   mat = parent[0];
//   pat = parent[1];
//
//   IntegerVector matalle = mat[0];
//   NumericVector matloc  = mat[1];
//
//   IntegerVector patalle = pat[0];
//   NumericVector patloc  = pat[1];
//
//   double L = max(matloc);
//   if(fabs(L - max(patloc)) > tol)
//     throw std::range_error("parent's two chromosomes are not the same length");
//
//   // simulate crossover locations; add -1 to the beginning
//   // NumericVector tmp = sim_crossovers(L, m, p, obligate_chiasma, Lstar);
//   // NumericVector Xlocations = sim_crossovers(L, m, p, obligate_chiasma, Lstar);
//
//   // potential for speed up: do not copy this vector at all
//   // NumericVector product(tmp.size() + 1);
//   // product[0] = -1.0;
//   // std::copy(tmp.begin(), tmp.end(), product.begin()+1);
//
//   // int cur_allele = random_int(0, 1); // first allele (0 or 1)
//   int cur_allele = (int) R::runif(0, 2);
//
//   // int biggest_length = product.size() + matloc.size() + patloc.size();
//   int biggest_length = Xlocations.size() + matloc.size() + patloc.size();
//   NumericVector loc(biggest_length);
//   IntegerVector alle(biggest_length);
//   int curpos = 0;
//
//   if(Xlocations.size() == 0) {
//     if(cur_allele == 0) return mat;
//     else return pat;
//   }
//   else {
//     int ct = 0;
//     double leftLoc = -1;
//
//     // NumericVector parloc;
//     // IntegerVector paralle;
//     // IntegerVector matstart = IntegerVector::create(0);
//     // IntegerVector patstart = IntegerVector::create(0);
//     // IntegerVector parstart;
//     // for(int i = 0; i < Xlocations.size(); i++) {
//     //   double rightLoc = Xlocations[i];
//     //   if(cur_allele == 0) { // mat chr
//     //     parloc = matloc; paralle = matalle; parstart = matstart;
//     //   } else {
//     //     parloc = patloc; paralle = patalle; parstart = patstart;
//     //   }
//     //
//     //   int pa;
//     //   for(int j = parstart[0]; j < parloc.size(); j++) {
//     //     double pl = parloc[j];
//     //     pa = paralle[j];
//     //     if (leftLoc <= pl && pl < rightLoc) {
//     //       loc[ct] = pl;
//     //       alle[ct] = pa;
//     //       ct++;
//     //       parstart[0]++;
//     //     } else if (rightLoc < pl) {
//     //       break;
//     //     }
//     //   }
//     //   loc[ct] = rightLoc;
//     //   alle[ct] = pa;
//     //   cur_allele = (int)(1 - cur_allele);
//     //   ct++;
//     //   leftLoc = rightLoc;
//     // }
//
//     // original
//     NumericVector parloc;
//     IntegerVector paralle;
//     int i, j;
//     for(i = 0; i < Xlocations.size(); i++) {
//       double rightLoc = Xlocations[i];
//       if(cur_allele == 0) { // mat chr
//         parloc = matloc; paralle = matalle;
//       } else {
//         parloc = patloc; paralle = patalle;
//       }
//
//       for(j = 0; j < parloc.size(); j++) {
//         double pl = parloc[j];
//         // pa = paralle[j];
//         if (leftLoc <= pl && pl < rightLoc) {
//           loc[ct] = pl;
//           alle[ct] = paralle[j];
//           ct++;
//         } else if (rightLoc < pl) {
//           break;
//         }
//       }
//       loc[ct] = rightLoc;
//       alle[ct] = paralle[j];
//       cur_allele = 1 - cur_allele;
//       ct++;
//       leftLoc = rightLoc;
//     }
//
//     // optimize with: double lastxo = Xlocations[Xlocations.size() - 1];
//     double lastxo = max(Xlocations);
//
//     if(cur_allele == 0) { // mat chr
//       parloc = matloc; paralle = matalle;
//     } else {
//       parloc = patloc; paralle = patalle;
//     }
//
//     for(int j = 0; j < parloc.size(); j++) {
//       if(parloc[j] > lastxo) {
//         loc[ct] = parloc[j];
//         alle[ct] = paralle[j];
//         ct++;
//       }
//     }
//
//     curpos = ct;
//   }
//
//   if(curpos > 1) { // clean up repeated alleles
//
//     NumericVector loc_clean(curpos);
//     IntegerVector alle_clean(curpos);
//     loc_clean[0] = loc[0];
//     alle_clean[0] = alle[0];
//     int lastpos=0;
//
//     for(int i = 1; i < curpos; i++) {
//       if(alle_clean[lastpos] == alle[i]) {
//         loc_clean[lastpos] = loc[i];
//       }
//       else {
//         lastpos++;
//         loc_clean[lastpos] = loc[i];
//         alle_clean[lastpos] = alle[i];
//       }
//     }
//     curpos = lastpos + 1;
//     loc = loc_clean;
//     alle = alle_clean;
//   }
//
//   // copy over to short vectors
//   NumericVector loc_result(curpos);
//   IntegerVector alle_result(curpos);
//   std::copy(loc.begin(), loc.begin()+curpos, loc_result.begin());
//   std::copy(alle.begin(), alle.begin()+curpos, alle_result.begin());
//
//   return List::create(Named("alleles") = alle_result, Named("locations") = loc_result);
// }

// //' @export
// // [[Rcpp::export]]
// List sim_meiosis_bro(const List parent, const int m, const double p,
//                  const bool obligate_chiasma, const double Lstar)
// {
//     const double tol=1e-12; // for comparison of chr lengths in parents
//
//     List mat, pat;
//     mat = parent[0];
//     pat = parent[1];
//
//     IntegerVector matalle = mat[0];
//     NumericVector matloc  = mat[1];
//
//     IntegerVector patalle = pat[0];
//     NumericVector patloc  = pat[1];
//
//     double L = max(matloc);
//     if(fabs(L - max(patloc)) > tol)
//         throw std::range_error("parent's two chromosomes are not the same length");
//
//     // simulate crossover locations; add -1 to the beginning
//     NumericVector tmp = sim_crossovers_rcpp(L, m, p, obligate_chiasma, Lstar);
//
//     NumericVector product(tmp.size() + 1);
//     product[0] = -1.0;
//     std::copy(tmp.begin(), tmp.end(), product.begin()+1);
//
//     int cur_allele = random_int(0, 1); // first allele (0 or 1)
//
//     int biggest_length = product.size() + matloc.size() + patloc.size();
//     NumericVector loc(biggest_length);
//     IntegerVector alle(biggest_length);
//
//     int curpos = 0;
//     if(product.size()==1) {
//         if(cur_allele==0) return mat;
//         else return pat;
//     }
//     else {
//         int i, j;
//         for(i=1; i<product.size(); i++) {
//
//             if(cur_allele==0) { // mat chr
//                 for(j=0; j<matloc.size(); j++) {
//                     if(matloc[j] >= product[i-1] && matloc[j] < product[i]) {
//                         loc[curpos] = matloc[j];
//                         alle[curpos] = matalle[j];
//                         curpos++;
//                     }
//                     else if(matloc[j] > product[i]) break;
//                 }
//                 loc[curpos] = product[i];
//                 alle[curpos] = matalle[j];
//                 curpos++;
//             }
//             else { // pat chr
//                 for(j=0; j<patloc.size(); j++) {
//                     if(patloc[j] >= product[i-1] && patloc[j] < product[i]) {
//                         loc[curpos] = patloc[j];
//                         alle[curpos] = patalle[j];
//                         curpos++;
//                     }
//                     else if(patloc[j] > product[i]) break;
//                 }
//                 loc[curpos] = product[i];
//                 alle[curpos] = patalle[j];
//                 curpos++;
//             }
//
//             cur_allele = 1 - cur_allele;
//
//         }
//
//         double lastxo = max(product);
//
//         if(cur_allele==0) { // mat chr
//             for(j=0; j<matloc.size(); j++) {
//                 if(matloc[j] > lastxo) {
//                     loc[curpos] = matloc[j];
//                     alle[curpos] = matalle[j];
//                     curpos++;
//                 }
//             }
//         }
//         else { // pat chr
//             for(j=0; j<patloc.size(); j++) {
//                 if(patloc[j] > lastxo) {
//                     loc[curpos] = patloc[j];
//                     alle[curpos] = patalle[j];
//                     curpos++;
//                 }
//             }
//         }
//     }
//
//     if(curpos > 1) { // clean up repeated alleles
//
//         NumericVector loc_clean(curpos);
//         IntegerVector alle_clean(curpos);
//
//         loc_clean[0] = loc[0];
//         alle_clean[0] = alle[0];
//         int lastpos=0;
//
//         for(int i=1; i<curpos; i++) {
//             if(alle_clean[lastpos] == alle[i]) {
//                 loc_clean[lastpos] = loc[i];
//             }
//             else {
//                 lastpos++;
//                 loc_clean[lastpos] = loc[i];
//                 alle_clean[lastpos] = alle[i];
//             }
//         }
//         curpos = lastpos+1;
//         loc = loc_clean;
//         alle = alle_clean;
//     }
//
//     // copy over to short vectors
//     NumericVector loc_result(curpos);
//     IntegerVector alle_result(curpos);
//     std::copy(loc.begin(), loc.begin()+curpos, loc_result.begin());
//     std::copy(alle.begin(), alle.begin()+curpos, alle_result.begin());
//
//     return List::create(Named("alleles")= alle_result, Named("locations")=loc_result);
// }

// // [[Rcpp::export]]
// List sim_meiosis_bro_loc(const List parent, NumericVector Xlocations)
// {
//     const double tol=1e-12; // for comparison of chr lengths in parents
//
//     List mat, pat;
//     mat = parent[0];
//     pat = parent[1];
//
//     IntegerVector matalle = mat[0];
//     NumericVector matloc  = mat[1];
//
//     IntegerVector patalle = pat[0];
//     NumericVector patloc  = pat[1];
//
//     double L = max(matloc);
//     if(fabs(L - max(patloc)) > tol)
//         throw std::range_error("parent's two chromosomes are not the same length");
//
//     // simulate crossover locations; add -1 to the beginning
//     NumericVector tmp = Xlocations;
//     // NumericVector tmp = sim_crossovers(L, m, p, obligate_chiasma, Lstar);
//
//     NumericVector product(tmp.size() + 1);
//     product[0] = -1.0;
//     std::copy(tmp.begin(), tmp.end(), product.begin()+1);
//
//     int cur_allele = random_int(0, 1); // first allele (0 or 1)
//
//     int biggest_length = product.size() + matloc.size() + patloc.size();
//     NumericVector loc(biggest_length);
//     IntegerVector alle(biggest_length);
//
//     int curpos = 0;
//     if(product.size()==1) {
//         if(cur_allele==0) return mat;
//         else return pat;
//     }
//     else {
//         int i, j;
//         for(i=1; i<product.size(); i++) {
//
//             if(cur_allele==0) { // mat chr
//                 for(j=0; j<matloc.size(); j++) {
//                     if(matloc[j] >= product[i-1] && matloc[j] < product[i]) {
//                         loc[curpos] = matloc[j];
//                         alle[curpos] = matalle[j];
//                         curpos++;
//                     }
//                     else if(matloc[j] > product[i]) break;
//                 }
//                 loc[curpos] = product[i];
//                 alle[curpos] = matalle[j];
//                 curpos++;
//             }
//             else { // pat chr
//                 for(j=0; j<patloc.size(); j++) {
//                     if(patloc[j] >= product[i-1] && patloc[j] < product[i]) {
//                         loc[curpos] = patloc[j];
//                         alle[curpos] = patalle[j];
//                         curpos++;
//                     }
//                     else if(patloc[j] > product[i]) break;
//                 }
//                 loc[curpos] = product[i];
//                 alle[curpos] = patalle[j];
//                 curpos++;
//             }
//
//             cur_allele = 1 - cur_allele;
//
//         }
//
//         double lastxo = max(product);
//
//         if(cur_allele==0) { // mat chr
//             for(j=0; j<matloc.size(); j++) {
//                 if(matloc[j] > lastxo) {
//                     loc[curpos] = matloc[j];
//                     alle[curpos] = matalle[j];
//                     curpos++;
//                 }
//             }
//         }
//         else { // pat chr
//             for(j=0; j<patloc.size(); j++) {
//                 if(patloc[j] > lastxo) {
//                     loc[curpos] = patloc[j];
//                     alle[curpos] = patalle[j];
//                     curpos++;
//                 }
//             }
//         }
//     }
//
//     if(curpos > 1) { // clean up repeated alleles
//
//         NumericVector loc_clean(curpos);
//         IntegerVector alle_clean(curpos);
//
//         loc_clean[0] = loc[0];
//         alle_clean[0] = alle[0];
//         int lastpos=0;
//
//         for(int i=1; i<curpos; i++) {
//             if(alle_clean[lastpos] == alle[i]) {
//                 loc_clean[lastpos] = loc[i];
//             }
//             else {
//                 lastpos++;
//                 loc_clean[lastpos] = loc[i];
//                 alle_clean[lastpos] = alle[i];
//             }
//         }
//         curpos = lastpos+1;
//         loc = loc_clean;
//         alle = alle_clean;
//     }
//
//     // copy over to short vectors
//     NumericVector loc_result(curpos);
//     IntegerVector alle_result(curpos);
//     std::copy(loc.begin(), loc.begin()+curpos, loc_result.begin());
//     std::copy(alle.begin(), alle.begin()+curpos, alle_result.begin());
//
//     return List::create(Named("alleles")= alle_result, Named("locations")=loc_result);
// }
//
//
//
//
//
// // [[Rcpp::export('.sim_meiosis_rcpp')]]
// List sim_meiosis_rcpp(const List parent,
//                          NumericVector Xlocations,
//                          bool check = true)
// {
//     const double tol = 1e-12; // for comparison of chr lengths in parents
//
//     List mat, pat;
//     mat = parent[0];
//     pat = parent[1];
//
//     IntegerVector matalle = mat[0];
//     NumericVector matloc  = mat[1];
//
//     IntegerVector patalle = pat[0];
//     NumericVector patloc  = pat[1];
//
//     // double L = max(matloc);
//     // if(fabs(L - max(patloc)) > tol)
//         // throw std::range_error("parent's two chromosomes are not the same length");
//
//     // We natuerally assume the location vectors are sorted.
//     int matloc_size = matloc.size();
//     int patloc_size = patloc.size();
//     double L = matloc[matloc_size - 1];
//     if (check) {
//       if(fabs(L - patloc[patloc_size - 1]) > tol)
//         throw std::range_error("parent's two chromosomes are not of the same length");
//     }
//
//     // simulate crossover locations; add -1 to the beginning
//     NumericVector tmp = Xlocations;
//     // NumericVector tmp = sim_crossovers(L, m, p, obligate_chiasma, Lstar);
//
//     int product_size = tmp.size() + 1;
//     NumericVector product(product_size);
//     product[0] = -1.0;
//     std::copy(tmp.begin(), tmp.end(), product.begin() + 1);
//
//     // int cur_allele = random_int(0, 1); // first allele (0 or 1)
//     int cur_allele = (int) R::runif(0, 2);
//
//     int biggest_length = product_size + matloc_size + patloc_size;
//     NumericVector loc(biggest_length);
//     IntegerVector alle(biggest_length);
//
//     int curpos = 0;
//     if(product_size == 1) {
//         if(cur_allele == 0) return mat;
//         else return pat;
//     }
//     else {
//         int i, j;
//         int jmatstart = 0;
//         int jpatstart = 0;
//         for(i = 1; i < product_size; i++) {
//
//             if(cur_allele == 0) { // mat chr
//                 for(j = jmatstart; j < matloc.size(); j++) {
//                     if(matloc[j] >= product[i-1] && matloc[j] < product[i]) {
//                         loc[curpos] = matloc[j];
//                         alle[curpos] = matalle[j];
//                         curpos++;
//                         jmatstart++;
//                     }
//                     else if(matloc[j] > product[i]) break;
//                 }
//                 loc[curpos] = product[i];
//                 alle[curpos] = matalle[j];
//                 curpos++;
//             }
//             else { // pat chr
//                 for(j = jpatstart; j < patloc.size(); j++) {
//                     if(patloc[j] >= product[i-1] && patloc[j] < product[i]) {
//                         loc[curpos] = patloc[j];
//                         alle[curpos] = patalle[j];
//                         curpos++;
//                         jpatstart++;
//                     }
//                     else if(patloc[j] > product[i]) break;
//                 }
//                 loc[curpos] = product[i];
//                 alle[curpos] = patalle[j];
//                 curpos++;
//             }
//
//             cur_allele = 1 - cur_allele;
//
//         }
//
//         // double lastxo = max(product);
//         double lastxo = product[product_size - 1];
//
//         if(cur_allele == 0) { // mat chr
//             for(j = 0; j < matloc_size; j++) {
//                 if(matloc[j] > lastxo) {
//                     loc[curpos] = matloc[j];
//                     alle[curpos] = matalle[j];
//                     curpos++;
//                 }
//             }
//         }
//         else { // pat chr
//             for(j = 0; j < patloc_size; j++) {
//                 if(patloc[j] > lastxo) {
//                     loc[curpos] = patloc[j];
//                     alle[curpos] = patalle[j];
//                     curpos++;
//                 }
//             }
//         }
//     }
//
//     if(curpos > 1) { // clean up repeated alleles
//
//         NumericVector loc_clean(curpos);
//         IntegerVector alle_clean(curpos);
//
//         loc_clean[0] = loc[0];
//         alle_clean[0] = alle[0];
//         int lastpos = 0;
//
//         for(int i = 1; i < curpos; i++) {
//             if(alle_clean[lastpos] == alle[i]) {
//                 loc_clean[lastpos] = loc[i];
//             }
//             else {
//                 lastpos++;
//                 loc_clean[lastpos] = loc[i];
//                 alle_clean[lastpos] = alle[i];
//             }
//         }
//         curpos = lastpos + 1;
//         loc = loc_clean;
//         alle = alle_clean;
//     }
//
//     // copy over to short vectors
//     NumericVector loc_result(curpos);
//     IntegerVector alle_result(curpos);
//     std::copy(loc.begin(), loc.begin() + curpos, loc_result.begin());
//     std::copy(alle.begin(), alle.begin() + curpos, alle_result.begin());
//
//     return List::create(Named("alleles") = alle_result,
//                         Named("locations") = loc_result);
// }
