// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// chromatid_value
double chromatid_value(const List& xodat, const arma::vec& map, const arma::mat& founder, const arma::vec& eff);
RcppExport SEXP simbreed_chromatid_value(SEXP xodatSEXP, SEXP mapSEXP, SEXP founderSEXP, SEXP effSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const List& >::type xodat(xodatSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type map(mapSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type founder(founderSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type eff(effSEXP);
    rcpp_result_gen = Rcpp::wrap(chromatid_value(xodat, map, founder, eff));
    return rcpp_result_gen;
END_RCPP
}
// gamete_value
double gamete_value(const Rcpp::List& xodat, const Rcpp::List& map, const Rcpp::List& founder, const Rcpp::List& eff);
RcppExport SEXP simbreed_gamete_value(SEXP xodatSEXP, SEXP mapSEXP, SEXP founderSEXP, SEXP effSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::List& >::type xodat(xodatSEXP);
    Rcpp::traits::input_parameter< const Rcpp::List& >::type map(mapSEXP);
    Rcpp::traits::input_parameter< const Rcpp::List& >::type founder(founderSEXP);
    Rcpp::traits::input_parameter< const Rcpp::List& >::type eff(effSEXP);
    rcpp_result_gen = Rcpp::wrap(gamete_value(xodat, map, founder, eff));
    return rcpp_result_gen;
END_RCPP
}
// gamete_value2
double gamete_value2(const List& parent, const List& params, const Rcpp::List& map, const Rcpp::List& founder, const Rcpp::List& eff);
RcppExport SEXP simbreed_gamete_value2(SEXP parentSEXP, SEXP paramsSEXP, SEXP mapSEXP, SEXP founderSEXP, SEXP effSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const List& >::type parent(parentSEXP);
    Rcpp::traits::input_parameter< const List& >::type params(paramsSEXP);
    Rcpp::traits::input_parameter< const Rcpp::List& >::type map(mapSEXP);
    Rcpp::traits::input_parameter< const Rcpp::List& >::type founder(founderSEXP);
    Rcpp::traits::input_parameter< const Rcpp::List& >::type eff(effSEXP);
    rcpp_result_gen = Rcpp::wrap(gamete_value2(parent, params, map, founder, eff));
    return rcpp_result_gen;
END_RCPP
}
// bcgv
Rcpp::List bcgv(const List& parent, const int n_gam, const int n_rep, const double se_thresh, const List& params, const Rcpp::List& map, const Rcpp::List& founder, const Rcpp::List& eff);
RcppExport SEXP simbreed_bcgv(SEXP parentSEXP, SEXP n_gamSEXP, SEXP n_repSEXP, SEXP se_threshSEXP, SEXP paramsSEXP, SEXP mapSEXP, SEXP founderSEXP, SEXP effSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const List& >::type parent(parentSEXP);
    Rcpp::traits::input_parameter< const int >::type n_gam(n_gamSEXP);
    Rcpp::traits::input_parameter< const int >::type n_rep(n_repSEXP);
    Rcpp::traits::input_parameter< const double >::type se_thresh(se_threshSEXP);
    Rcpp::traits::input_parameter< const List& >::type params(paramsSEXP);
    Rcpp::traits::input_parameter< const Rcpp::List& >::type map(mapSEXP);
    Rcpp::traits::input_parameter< const Rcpp::List& >::type founder(founderSEXP);
    Rcpp::traits::input_parameter< const Rcpp::List& >::type eff(effSEXP);
    rcpp_result_gen = Rcpp::wrap(bcgv(parent, n_gam, n_rep, se_thresh, params, map, founder, eff));
    return rcpp_result_gen;
END_RCPP
}
// xo2geno_chromatid
arma::rowvec xo2geno_chromatid(const List& xodat, const arma::vec& map, const arma::mat& founder);
RcppExport SEXP simbreed_xo2geno_chromatid(SEXP xodatSEXP, SEXP mapSEXP, SEXP founderSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const List& >::type xodat(xodatSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type map(mapSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type founder(founderSEXP);
    rcpp_result_gen = Rcpp::wrap(xo2geno_chromatid(xodat, map, founder));
    return rcpp_result_gen;
END_RCPP
}
// xo2geno_chromosome
arma::mat xo2geno_chromosome(const Rcpp::List& xodat, const arma::vec& map, const arma::mat& founder, const bool homozygous);
RcppExport SEXP simbreed_xo2geno_chromosome(SEXP xodatSEXP, SEXP mapSEXP, SEXP founderSEXP, SEXP homozygousSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::List& >::type xodat(xodatSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type map(mapSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type founder(founderSEXP);
    Rcpp::traits::input_parameter< const bool >::type homozygous(homozygousSEXP);
    rcpp_result_gen = Rcpp::wrap(xo2geno_chromosome(xodat, map, founder, homozygous));
    return rcpp_result_gen;
END_RCPP
}
// xo2geno_individual
arma::mat xo2geno_individual(const Rcpp::List& xodat, const Rcpp::List& map, const Rcpp::List& founder);
RcppExport SEXP simbreed_xo2geno_individual(SEXP xodatSEXP, SEXP mapSEXP, SEXP founderSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::List& >::type xodat(xodatSEXP);
    Rcpp::traits::input_parameter< const Rcpp::List& >::type map(mapSEXP);
    Rcpp::traits::input_parameter< const Rcpp::List& >::type founder(founderSEXP);
    rcpp_result_gen = Rcpp::wrap(xo2geno_individual(xodat, map, founder));
    return rcpp_result_gen;
END_RCPP
}
// xo2geno_gamete
arma::mat xo2geno_gamete(const Rcpp::List& xodat, const Rcpp::List& map, const Rcpp::List& founder);
RcppExport SEXP simbreed_xo2geno_gamete(SEXP xodatSEXP, SEXP mapSEXP, SEXP founderSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::List& >::type xodat(xodatSEXP);
    Rcpp::traits::input_parameter< const Rcpp::List& >::type map(mapSEXP);
    Rcpp::traits::input_parameter< const Rcpp::List& >::type founder(founderSEXP);
    rcpp_result_gen = Rcpp::wrap(xo2geno_gamete(xodat, map, founder));
    return rcpp_result_gen;
END_RCPP
}
// xo2geno_population
arma::mat xo2geno_population(const Rcpp::List& xodat, const Rcpp::List& map, const Rcpp::List& founder);
RcppExport SEXP simbreed_xo2geno_population(SEXP xodatSEXP, SEXP mapSEXP, SEXP founderSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::List& >::type xodat(xodatSEXP);
    Rcpp::traits::input_parameter< const Rcpp::List& >::type map(mapSEXP);
    Rcpp::traits::input_parameter< const Rcpp::List& >::type founder(founderSEXP);
    rcpp_result_gen = Rcpp::wrap(xo2geno_population(xodat, map, founder));
    return rcpp_result_gen;
END_RCPP
}
// xo2geno
Rcpp::List xo2geno(const Rcpp::List& xodat, const Rcpp::List& map, const Rcpp::List& founder);
RcppExport SEXP simbreed_xo2geno(SEXP xodatSEXP, SEXP mapSEXP, SEXP founderSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::List& >::type xodat(xodatSEXP);
    Rcpp::traits::input_parameter< const Rcpp::List& >::type map(mapSEXP);
    Rcpp::traits::input_parameter< const Rcpp::List& >::type founder(founderSEXP);
    rcpp_result_gen = Rcpp::wrap(xo2geno(xodat, map, founder));
    return rcpp_result_gen;
END_RCPP
}
// gamete
Rcpp::List gamete(const List& parent, const List& params);
RcppExport SEXP simbreed_gamete(SEXP parentSEXP, SEXP paramsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const List& >::type parent(parentSEXP);
    Rcpp::traits::input_parameter< const List& >::type params(paramsSEXP);
    rcpp_result_gen = Rcpp::wrap(gamete(parent, params));
    return rcpp_result_gen;
END_RCPP
}
// cross
Rcpp::List cross(const List& mother, const List& father, const List& params);
RcppExport SEXP simbreed_cross(SEXP motherSEXP, SEXP fatherSEXP, SEXP paramsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const List& >::type mother(motherSEXP);
    Rcpp::traits::input_parameter< const List& >::type father(fatherSEXP);
    Rcpp::traits::input_parameter< const List& >::type params(paramsSEXP);
    rcpp_result_gen = Rcpp::wrap(cross(mother, father, params));
    return rcpp_result_gen;
END_RCPP
}
// realized_IBD_chromatids
double realized_IBD_chromatids(List chromatid1, List chromatid2);
RcppExport SEXP simbreed_realized_IBD_chromatids(SEXP chromatid1SEXP, SEXP chromatid2SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< List >::type chromatid1(chromatid1SEXP);
    Rcpp::traits::input_parameter< List >::type chromatid2(chromatid2SEXP);
    rcpp_result_gen = Rcpp::wrap(realized_IBD_chromatids(chromatid1, chromatid2));
    return rcpp_result_gen;
END_RCPP
}
// realized_IBD
double realized_IBD(List ind1, List ind2);
RcppExport SEXP simbreed_realized_IBD(SEXP ind1SEXP, SEXP ind2SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< List >::type ind1(ind1SEXP);
    Rcpp::traits::input_parameter< List >::type ind2(ind2SEXP);
    rcpp_result_gen = Rcpp::wrap(realized_IBD(ind1, ind2));
    return rcpp_result_gen;
END_RCPP
}
// sim_crossovers
arma::vec sim_crossovers(const double L, const int m, const double p, const bool obligate_chiasma, const double Lstar);
RcppExport SEXP simbreed_sim_crossovers(SEXP LSEXP, SEXP mSEXP, SEXP pSEXP, SEXP obligate_chiasmaSEXP, SEXP LstarSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const double >::type L(LSEXP);
    Rcpp::traits::input_parameter< const int >::type m(mSEXP);
    Rcpp::traits::input_parameter< const double >::type p(pSEXP);
    Rcpp::traits::input_parameter< const bool >::type obligate_chiasma(obligate_chiasmaSEXP);
    Rcpp::traits::input_parameter< const double >::type Lstar(LstarSEXP);
    rcpp_result_gen = Rcpp::wrap(sim_crossovers(L, m, p, obligate_chiasma, Lstar));
    return rcpp_result_gen;
END_RCPP
}
// sim_meiosis_xo
List sim_meiosis_xo(const List parent, arma::vec& Xlocations);
RcppExport SEXP simbreed_sim_meiosis_xo(SEXP parentSEXP, SEXP XlocationsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const List >::type parent(parentSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type Xlocations(XlocationsSEXP);
    rcpp_result_gen = Rcpp::wrap(sim_meiosis_xo(parent, Xlocations));
    return rcpp_result_gen;
END_RCPP
}
// sim_meiosis_xo_test
List sim_meiosis_xo_test(const List parent, arma::vec& Xlocations);
RcppExport SEXP simbreed_sim_meiosis_xo_test(SEXP parentSEXP, SEXP XlocationsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const List >::type parent(parentSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type Xlocations(XlocationsSEXP);
    rcpp_result_gen = Rcpp::wrap(sim_meiosis_xo_test(parent, Xlocations));
    return rcpp_result_gen;
END_RCPP
}
// sim_meiosis_reliable
Rcpp::List sim_meiosis_reliable(const Rcpp::List& parent, const double L, const int m, const double p, const bool obligate_chiasma, const double Lstar);
RcppExport SEXP simbreed_sim_meiosis_reliable(SEXP parentSEXP, SEXP LSEXP, SEXP mSEXP, SEXP pSEXP, SEXP obligate_chiasmaSEXP, SEXP LstarSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::List& >::type parent(parentSEXP);
    Rcpp::traits::input_parameter< const double >::type L(LSEXP);
    Rcpp::traits::input_parameter< const int >::type m(mSEXP);
    Rcpp::traits::input_parameter< const double >::type p(pSEXP);
    Rcpp::traits::input_parameter< const bool >::type obligate_chiasma(obligate_chiasmaSEXP);
    Rcpp::traits::input_parameter< const double >::type Lstar(LstarSEXP);
    rcpp_result_gen = Rcpp::wrap(sim_meiosis_reliable(parent, L, m, p, obligate_chiasma, Lstar));
    return rcpp_result_gen;
END_RCPP
}
// sim_meiosis
Rcpp::List sim_meiosis(const Rcpp::List& parent, const double L, const int m, const double p, const bool obligate_chiasma, const double Lstar);
RcppExport SEXP simbreed_sim_meiosis(SEXP parentSEXP, SEXP LSEXP, SEXP mSEXP, SEXP pSEXP, SEXP obligate_chiasmaSEXP, SEXP LstarSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::List& >::type parent(parentSEXP);
    Rcpp::traits::input_parameter< const double >::type L(LSEXP);
    Rcpp::traits::input_parameter< const int >::type m(mSEXP);
    Rcpp::traits::input_parameter< const double >::type p(pSEXP);
    Rcpp::traits::input_parameter< const bool >::type obligate_chiasma(obligate_chiasmaSEXP);
    Rcpp::traits::input_parameter< const double >::type Lstar(LstarSEXP);
    rcpp_result_gen = Rcpp::wrap(sim_meiosis(parent, L, m, p, obligate_chiasma, Lstar));
    return rcpp_result_gen;
END_RCPP
}
