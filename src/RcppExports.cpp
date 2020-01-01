// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// vaccC
NumericVector vaccC(NumericVector age, LogicalVector female, LogicalVector ily);
RcppExport SEXP _SC19003_vaccC(SEXP ageSEXP, SEXP femaleSEXP, SEXP ilySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type age(ageSEXP);
    Rcpp::traits::input_parameter< LogicalVector >::type female(femaleSEXP);
    Rcpp::traits::input_parameter< LogicalVector >::type ily(ilySEXP);
    rcpp_result_gen = Rcpp::wrap(vaccC(age, female, ily));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_SC19003_vaccC", (DL_FUNC) &_SC19003_vaccC, 3},
    {NULL, NULL, 0}
};

RcppExport void R_init_SC19003(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
