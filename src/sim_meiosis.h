#ifndef SIM_MEIOSIS_H
#define SIM_MEIOSIS_H

Rcpp::NumericVector sim_crossovers(const double L, const int m, const double p,
                                   const bool obligate_chiasma, const double Lstar);

Rcpp::List sim_meiosis(const Rcpp::List& parent,
                       const double L,
                       const int m,
                       const double p,
                       const bool obligate_chiasma,
                       const double Lstar);

Rcpp::List sim_meiosis_reliable(const Rcpp::List& parent,
                       const double L,
                       const int m,
                       const double p,
                       const bool obligate_chiasma,
                       const double Lstar);

#endif
