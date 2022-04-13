// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// cpp_GC_content
NumericMatrix cpp_GC_content(StringMatrix x);
RcppExport SEXP _FastqCleaner_cpp_GC_content(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< StringMatrix >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(cpp_GC_content(x));
    return rcpp_result_gen;
END_RCPP
}
// cpp_base_content
NumericMatrix cpp_base_content(StringMatrix x, bool relative);
RcppExport SEXP _FastqCleaner_cpp_base_content(SEXP xSEXP, SEXP relativeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< StringMatrix >::type x(xSEXP);
    Rcpp::traits::input_parameter< bool >::type relative(relativeSEXP);
    rcpp_result_gen = Rcpp::wrap(cpp_base_content(x, relative));
    return rcpp_result_gen;
END_RCPP
}
// cpp_check_quality_threshold
NumericMatrix cpp_check_quality_threshold(NumericMatrix x, NumericVector thresVector);
RcppExport SEXP _FastqCleaner_cpp_check_quality_threshold(SEXP xSEXP, SEXP thresVectorSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type x(xSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type thresVector(thresVectorSEXP);
    rcpp_result_gen = Rcpp::wrap(cpp_check_quality_threshold(x, thresVector));
    return rcpp_result_gen;
END_RCPP
}
// cpp_create_stringvec
StringVector cpp_create_stringvec(NumericVector vec, char x);
RcppExport SEXP _FastqCleaner_cpp_create_stringvec(SEXP vecSEXP, SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type vec(vecSEXP);
    Rcpp::traits::input_parameter< char >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(cpp_create_stringvec(vec, x));
    return rcpp_result_gen;
END_RCPP
}
// cpp_which_true
IntegerVector cpp_which_true(LogicalMatrix m, IntegerVector width, String origin);
RcppExport SEXP _FastqCleaner_cpp_which_true(SEXP mSEXP, SEXP widthSEXP, SEXP originSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< LogicalMatrix >::type m(mSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type width(widthSEXP);
    Rcpp::traits::input_parameter< String >::type origin(originSEXP);
    rcpp_result_gen = Rcpp::wrap(cpp_which_true(m, width, origin));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_FastqCleaner_cpp_GC_content", (DL_FUNC) &_FastqCleaner_cpp_GC_content, 1},
    {"_FastqCleaner_cpp_base_content", (DL_FUNC) &_FastqCleaner_cpp_base_content, 2},
    {"_FastqCleaner_cpp_check_quality_threshold", (DL_FUNC) &_FastqCleaner_cpp_check_quality_threshold, 2},
    {"_FastqCleaner_cpp_create_stringvec", (DL_FUNC) &_FastqCleaner_cpp_create_stringvec, 2},
    {"_FastqCleaner_cpp_which_true", (DL_FUNC) &_FastqCleaner_cpp_which_true, 3},
    {NULL, NULL, 0}
};

RcppExport void R_init_FastqCleaner(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
