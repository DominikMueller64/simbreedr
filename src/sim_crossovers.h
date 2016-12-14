#ifndef SIM_CROSSOVER_H
#define SIM_CROSSOVER_H
arma::vec sim_crossovers(const double L,
                         const int m,
                         const double p,
                         const bool obligate_chiasma,
                         const double Lstar);

// Rcpp::NumericVector sim_crossovers_rcpp(const double L,
//                                         const int m,
//                                         const double p,
//                                         const bool obligate_chiasma,
//                                         const double Lstar);

#endif
