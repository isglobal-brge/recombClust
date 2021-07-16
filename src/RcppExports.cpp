// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppEigen.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// BigLD
Rcpp::RObject BigLD(Rcpp::RObject corMat, Nullable<double> CLQcut, Nullable<int> clstgap, std::string CLQmode, std::string hrstType, Nullable<int> hrstParam, std::string chrN);
RcppExport SEXP _recombClust_BigLD(SEXP corMatSEXP, SEXP CLQcutSEXP, SEXP clstgapSEXP, SEXP CLQmodeSEXP, SEXP hrstTypeSEXP, SEXP hrstParamSEXP, SEXP chrNSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::RObject >::type corMat(corMatSEXP);
    Rcpp::traits::input_parameter< Nullable<double> >::type CLQcut(CLQcutSEXP);
    Rcpp::traits::input_parameter< Nullable<int> >::type clstgap(clstgapSEXP);
    Rcpp::traits::input_parameter< std::string >::type CLQmode(CLQmodeSEXP);
    Rcpp::traits::input_parameter< std::string >::type hrstType(hrstTypeSEXP);
    Rcpp::traits::input_parameter< Nullable<int> >::type hrstParam(hrstParamSEXP);
    Rcpp::traits::input_parameter< std::string >::type chrN(chrNSEXP);
    rcpp_result_gen = Rcpp::wrap(BigLD(corMat, CLQcut, clstgap, CLQmode, hrstType, hrstParam, chrN));
    return rcpp_result_gen;
END_RCPP
}
// CLQD_mod
Rcpp::RObject CLQD_mod(Rcpp::RObject CorMat, double CLQcut, int clstgap, std::string hrstType, int hrstParam, std::string CLQmode);
RcppExport SEXP _recombClust_CLQD_mod(SEXP CorMatSEXP, SEXP CLQcutSEXP, SEXP clstgapSEXP, SEXP hrstTypeSEXP, SEXP hrstParamSEXP, SEXP CLQmodeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::RObject >::type CorMat(CorMatSEXP);
    Rcpp::traits::input_parameter< double >::type CLQcut(CLQcutSEXP);
    Rcpp::traits::input_parameter< int >::type clstgap(clstgapSEXP);
    Rcpp::traits::input_parameter< std::string >::type hrstType(hrstTypeSEXP);
    Rcpp::traits::input_parameter< int >::type hrstParam(hrstParamSEXP);
    Rcpp::traits::input_parameter< std::string >::type CLQmode(CLQmodeSEXP);
    rcpp_result_gen = Rcpp::wrap(CLQD_mod(CorMat, CLQcut, clstgap, hrstType, hrstParam, CLQmode));
    return rcpp_result_gen;
END_RCPP
}
// CWriteResults
void CWriteResults(StringVector columns, NumericVector results, std::string outputfile);
RcppExport SEXP _recombClust_CWriteResults(SEXP columnsSEXP, SEXP resultsSEXP, SEXP outputfileSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< StringVector >::type columns(columnsSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type results(resultsSEXP);
    Rcpp::traits::input_parameter< std::string >::type outputfile(outputfileSEXP);
    CWriteResults(columns, results, outputfile);
    return R_NilValue;
END_RCPP
}
// get_graph_matrix_data
Rcpp::RObject get_graph_matrix_data(Rcpp::RObject OCM, int hrstParam);
RcppExport SEXP _recombClust_get_graph_matrix_data(SEXP OCMSEXP, SEXP hrstParamSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::RObject >::type OCM(OCMSEXP);
    Rcpp::traits::input_parameter< int >::type hrstParam(hrstParamSEXP);
    rcpp_result_gen = Rcpp::wrap(get_graph_matrix_data(OCM, hrstParam));
    return rcpp_result_gen;
END_RCPP
}
// LDmixtureModel
List LDmixtureModel(RObject dat, Nullable<int> maxSteps, Nullable<double> prob0, Nullable<int> blocksize);
RcppExport SEXP _recombClust_LDmixtureModel(SEXP datSEXP, SEXP maxStepsSEXP, SEXP prob0SEXP, SEXP blocksizeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< RObject >::type dat(datSEXP);
    Rcpp::traits::input_parameter< Nullable<int> >::type maxSteps(maxStepsSEXP);
    Rcpp::traits::input_parameter< Nullable<double> >::type prob0(prob0SEXP);
    Rcpp::traits::input_parameter< Nullable<int> >::type blocksize(blocksizeSEXP);
    rcpp_result_gen = Rcpp::wrap(LDmixtureModel(dat, maxSteps, prob0, blocksize));
    return rcpp_result_gen;
END_RCPP
}
// removeMatrixColsandRows
Rcpp::RObject removeMatrixColsandRows(Rcpp::RObject Mat, Rcpp::RObject vIndex);
RcppExport SEXP _recombClust_removeMatrixColsandRows(SEXP MatSEXP, SEXP vIndexSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::RObject >::type Mat(MatSEXP);
    Rcpp::traits::input_parameter< Rcpp::RObject >::type vIndex(vIndexSEXP);
    rcpp_result_gen = Rcpp::wrap(removeMatrixColsandRows(Mat, vIndex));
    return rcpp_result_gen;
END_RCPP
}
// CGetDatafromFile
std::string CGetDatafromFile(std::string file);
RcppExport SEXP _recombClust_CGetDatafromFile(SEXP fileSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string >::type file(fileSEXP);
    rcpp_result_gen = Rcpp::wrap(CGetDatafromFile(file));
    return rcpp_result_gen;
END_RCPP
}
// CgdsSNPpairMatrix
Rcpp::NumericMatrix CgdsSNPpairMatrix(Rcpp::StringMatrix filteredsnp);
RcppExport SEXP _recombClust_CgdsSNPpairMatrix(SEXP filteredsnpSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::StringMatrix >::type filteredsnp(filteredsnpSEXP);
    rcpp_result_gen = Rcpp::wrap(CgdsSNPpairMatrix(filteredsnp));
    return rcpp_result_gen;
END_RCPP
}
// CvcfSNPpairMatrix
Rcpp::NumericMatrix CvcfSNPpairMatrix(Rcpp::StringMatrix filteredsnp);
RcppExport SEXP _recombClust_CvcfSNPpairMatrix(SEXP filteredsnpSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::StringMatrix >::type filteredsnp(filteredsnpSEXP);
    rcpp_result_gen = Rcpp::wrap(CvcfSNPpairMatrix(filteredsnp));
    return rcpp_result_gen;
END_RCPP
}
// CTransformtoSampleAlleles
Rcpp::IntegerMatrix CTransformtoSampleAlleles(Rcpp::NumericVector x, Rcpp::StringVector allele, Rcpp::StringVector sample, Rcpp::StringVector variable);
RcppExport SEXP _recombClust_CTransformtoSampleAlleles(SEXP xSEXP, SEXP alleleSEXP, SEXP sampleSEXP, SEXP variableSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< Rcpp::StringVector >::type allele(alleleSEXP);
    Rcpp::traits::input_parameter< Rcpp::StringVector >::type sample(sampleSEXP);
    Rcpp::traits::input_parameter< Rcpp::StringVector >::type variable(variableSEXP);
    rcpp_result_gen = Rcpp::wrap(CTransformtoSampleAlleles(x, allele, sample, variable));
    return rcpp_result_gen;
END_RCPP
}
// getCorrelationMatrix
Rcpp::RObject getCorrelationMatrix(Rcpp::RObject mat, Rcpp::Nullable<bool> absval);
RcppExport SEXP _recombClust_getCorrelationMatrix(SEXP matSEXP, SEXP absvalSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::RObject >::type mat(matSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<bool> >::type absval(absvalSEXP);
    rcpp_result_gen = Rcpp::wrap(getCorrelationMatrix(mat, absval));
    return rcpp_result_gen;
END_RCPP
}
// getProbs
Rcpp::NumericVector getProbs(Rcpp::RObject mat, Rcpp::RObject sel);
RcppExport SEXP _recombClust_getProbs(SEXP matSEXP, SEXP selSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::RObject >::type mat(matSEXP);
    Rcpp::traits::input_parameter< Rcpp::RObject >::type sel(selSEXP);
    rcpp_result_gen = Rcpp::wrap(getProbs(mat, sel));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_recombClust_BigLD", (DL_FUNC) &_recombClust_BigLD, 7},
    {"_recombClust_CLQD_mod", (DL_FUNC) &_recombClust_CLQD_mod, 6},
    {"_recombClust_CWriteResults", (DL_FUNC) &_recombClust_CWriteResults, 3},
    {"_recombClust_get_graph_matrix_data", (DL_FUNC) &_recombClust_get_graph_matrix_data, 2},
    {"_recombClust_LDmixtureModel", (DL_FUNC) &_recombClust_LDmixtureModel, 4},
    {"_recombClust_removeMatrixColsandRows", (DL_FUNC) &_recombClust_removeMatrixColsandRows, 2},
    {"_recombClust_CGetDatafromFile", (DL_FUNC) &_recombClust_CGetDatafromFile, 1},
    {"_recombClust_CgdsSNPpairMatrix", (DL_FUNC) &_recombClust_CgdsSNPpairMatrix, 1},
    {"_recombClust_CvcfSNPpairMatrix", (DL_FUNC) &_recombClust_CvcfSNPpairMatrix, 1},
    {"_recombClust_CTransformtoSampleAlleles", (DL_FUNC) &_recombClust_CTransformtoSampleAlleles, 4},
    {"_recombClust_getCorrelationMatrix", (DL_FUNC) &_recombClust_getCorrelationMatrix, 2},
    {"_recombClust_getProbs", (DL_FUNC) &_recombClust_getProbs, 2},
    {NULL, NULL, 0}
};

RcppExport void R_init_recombClust(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
