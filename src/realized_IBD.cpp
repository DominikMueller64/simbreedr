// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;

//[[Rcpp::export(".realized_IBD_chromatids")]]
double realized_IBD_chromatids(List chromatid1, List chromatid2)
{
  double tot = 0;
  List mat = chromatid1;
  List pat = chromatid2;
  vec matloc = mat["locations"];
  vec patloc = pat["locations"];
  ivec matalle = mat["alleles"];
  ivec patalle = pat["alleles"];
  int ipatmax = patloc.n_elem;
  int imatmax = matloc.n_elem;
  int imat = 0;
  int ipat = 0;
  int cur_alle = (matloc[ipat] <= patloc[imat]);
  // Rcout << (matloc[ipat] <= patloc[imat]) << std::endl;
  // Rcout << cur_alle << std::endl;
  double left = 0;
  double right;

  while(true) {
    if (cur_alle) {
      do
      {
        right = matloc[imat];
        if (patalle[ipat] == matalle[imat]) {
          tot += (right - left);
        }
        // Rcout << left << "\t" << right << std::endl;
        imat++;
        left = right;
      } while (patloc[ipat] >= matloc[imat] && imat < imatmax);
      cur_alle = 1 - cur_alle;
    }
    else {
      do
      {
        right = patloc[ipat];
        if (patalle[ipat] == matalle[imat]) {
          tot += (right - left);
        }
        // Rcout << left << "\t" << right << std::endl;
        ipat++;
        left = right;
      } while (patloc[ipat] < matloc[imat] && ipat < ipatmax);
      cur_alle = 1 - cur_alle;
    }
    if (ipat == ipatmax || imat == imatmax ) {
      break;
    }
  }
  return tot;
}

//[[Rcpp::export(".realized_IBD")]]
double realized_IBD(List ind1, List ind2) {
  int n_chrom = ind1.size();
  double len = 0;
  double tot = 0;
  for (int i = 0; i < n_chrom; i++) {
    List ind1_i = ind1[i];
    List tmpl = ind1_i["mat"];
    vec tmp = tmpl["locations"];
    len += tmp[tmp.n_elem - 1];
    List ind2_i = ind2[i];
    tot += realized_IBD_chromatids(ind1_i["mat"], ind2_i["mat"]) +
      realized_IBD_chromatids(ind1_i["mat"], ind2_i["pat"]) +
      realized_IBD_chromatids(ind1_i["pat"], ind2_i["mat"]) +
      realized_IBD_chromatids(ind1_i["pat"], ind2_i["pat"]);
  }
  // Rcout << tot << std::endl;
  // Rcout << len << std::endl;
  return tot / (4 * len);
}
