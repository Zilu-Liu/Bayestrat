/*
 This file was autogenerated with the following R routine:
 tools::package_native_routine_registration_skeleton("/home/pperez/BLR")
*/

#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* .Call calls */
extern SEXP sample_beta(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"sample_beta", (DL_FUNC) &sample_beta, 9},
    {NULL, NULL, 0}
};

void R_init_BLR(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
