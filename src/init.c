#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* 
Routines registration obtained with 

tools::package_native_routine_registration_skeleton(".", character_only = FALSE)
 
FIXME: Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP _mixggm_conggm(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _mixggm_estepmggm(SEXP, SEXP, SEXP, SEXP);
extern SEXP _mixggm_icf(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _mixggm_profileloglik(SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"_mixggm_conggm",        (DL_FUNC) &_mixggm_conggm,        6},
    {"_mixggm_estepmggm",     (DL_FUNC) &_mixggm_estepmggm,     4},
    {"_mixggm_icf",           (DL_FUNC) &_mixggm_icf,           9},
    {"_mixggm_profileloglik", (DL_FUNC) &_mixggm_profileloglik, 3},
    {NULL, NULL, 0}
};

void R_init_mixggm(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}