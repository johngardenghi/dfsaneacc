#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>

void F77_NAME(bmqr)(int *n, int *m, int *mlim, int *lda, double *A, double *gam, int *p, int *rank);
void F77_NAME(extractQ)(int *n, int *m, int *mlim, int *lda, double *A, double *gam, int *rank, int *ldq, double *Q);
void F77_NAME(extractR)(int *m, int *mlim, int *lda, double *A, int *rank, int *ldr, double *R);
void F77_NAME(bmqrupdate)(int *n, int *m, int *rank, int *mlim, int *jtbr, int *newcolp, int *ldq, double *Q, int *ldr, double *R, int *p, double *newcol);
void F77_NAME(bmqrsolve)(int *n, int *m, int *rank, int *ldq, double *Q, int *ldr, double *R, int *p, double *b, double *x);

SEXP c_bmqr(SEXP R_n, SEXP R_m, SEXP R_mlim, SEXP R_lda, SEXP R_A, SEXP R_gam, SEXP R_p, SEXP R_rank) {
  F77_CALL(bmqr)(INTEGER(R_n), INTEGER(R_m), INTEGER(R_mlim),
		 INTEGER(R_lda), REAL(R_A), REAL(R_gam),
		 INTEGER(R_p), INTEGER(R_rank));
  return R_NilValue;
}

SEXP c_extractQ(SEXP R_n, SEXP R_m, SEXP R_mlim, SEXP R_lda, SEXP R_A, SEXP R_gam, SEXP R_rank, SEXP R_ldq, SEXP R_Q) {
  F77_CALL(extractQ)(INTEGER(R_n), INTEGER(R_m), INTEGER(R_mlim),
		     INTEGER(R_lda), REAL(R_A), REAL(R_gam),
		     INTEGER(R_rank), INTEGER(R_ldq), REAL(R_Q));
  return R_NilValue;
}

SEXP c_extractR(SEXP R_m, SEXP R_mlim, SEXP R_lda, SEXP R_A, SEXP R_rank, SEXP R_ldr, SEXP R_R) {
    F77_CALL(extractR)(INTEGER(R_m), INTEGER(R_mlim), INTEGER(R_lda),
		       REAL(R_A), INTEGER(R_rank), INTEGER(R_ldr), REAL(R_R));
  return R_NilValue;
}

SEXP c_bmqrupdate(SEXP R_n, SEXP R_m, SEXP R_rank, SEXP R_mlim, SEXP R_jtbr, SEXP R_newcolp, SEXP R_ldq, SEXP R_Q, SEXP R_ldr, SEXP R_R, SEXP R_p, SEXP R_newcol) {
  F77_CALL(bmqrupdate)(INTEGER(R_n), INTEGER(R_m), INTEGER(R_rank),
		       INTEGER(R_mlim), INTEGER(R_jtbr), INTEGER(R_newcolp), INTEGER(R_ldq),
		       REAL(R_Q), INTEGER(R_ldr), REAL(R_R), INTEGER(R_p),
		       REAL(R_newcol));
  return R_NilValue;
}

SEXP c_bmqrsolve(SEXP R_n, SEXP R_m, SEXP R_rank, SEXP R_ldq, SEXP R_Q, SEXP R_ldr, SEXP R_R, SEXP R_p, SEXP R_b, SEXP R_x) {
    F77_CALL(bmqrsolve)(INTEGER(R_n), INTEGER(R_m), INTEGER(R_rank),
			INTEGER(R_ldq), REAL(R_Q), INTEGER(R_ldr),
			REAL(R_R), INTEGER(R_p), REAL(R_b), REAL(R_x));
  return R_NilValue;
}

static const R_CallMethodDef CallEntries[] = {
  {"c_bmqr",         (DL_FUNC) &c_bmqr,         8},
  {"c_extractQ",     (DL_FUNC) &c_extractQ,     9},
  {"c_extractR",     (DL_FUNC) &c_extractR,     7},
  {"c_bmqrupdate",   (DL_FUNC) &c_bmqrupdate,  12},
  {"c_bmqrsolve",    (DL_FUNC) &c_bmqrsolve,   10},
  {NULL, NULL, 0}
};

void R_init_dfsaneacc (DllInfo *dll) {
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
  R_forceSymbols(dll, TRUE);
}
