#include <stdio.h>
#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include "cutest.h"

SEXP c_cutest_getn () {
  char *fname = "OUTSDIF.d"; /* CUTEst data file */
  integer funit = 42;        /* FORTRAN unit number for OUTSDIF.d */
  integer ierr;              /* Exit flag from OPEN and CLOSE */
  integer status;            /* Exit flag from CUTEst tools */
  integer CUTEst_ncon;       /* number of constraints */

  SEXP R_nvar = PROTECT( allocVector( INTSXP, 1 ) );

  /* Open problem description file OUTSDIF.d */
  FORTRAN_open( &funit, fname, &ierr );
  if( ierr ) {
    printf("cutest_getdims: Error opening file OUTSDIF.d.\nAborting.\n");
    exit(1);
  }

  /* Determine problem size */
  CUTEST_cdimen( &status, &funit, INTEGER(R_nvar), &CUTEst_ncon );
  if ( status ) {
    printf( "cutest_init error: cutest_cdimen - status not equal to zero. status = %d\n", status );
    exit( status );
  }

  /* Close problem description file OUTSDIF.d */
  FORTRAN_close( &funit, &ierr );

  UNPROTECT(1);

  return R_nvar;
}

SEXP c_cutest_init (SEXP R_x0, SEXP R_crhs) {
  char *fname = "OUTSDIF.d"; /* CUTEst data file */
  char pname[FSTRING_LEN+1]; /* Problem name */

  integer funit = 42;        /* FORTRAN unit number for OUTSDIF.d */
  integer iout = 6;          /* FORTRAN unit number for error output */
  integer io_buffer = 11;    /* FORTRAN unit internal input/output */
  integer ierr;              /* Exit flag from OPEN and CLOSE */
  integer status;            /* Exit flag from CUTEst tools */
  integer CUTEst_nvar;
  integer CUTEst_ncon;       /* number of constraints */

  doublereal *bl, *bu, *lambda, *cl, *crhs, f, *x0;
  integer i;
  integer e_order = 0, l_order = 0, v_order = 0;
  logical *equatn, *linear;

  /* Open problem description file OUTSDIF.d */
  FORTRAN_open( &funit, fname, &ierr );
  if( ierr ) {
    printf("cutest_init: Error opening file OUTSDIF.d.\nAborting.\n");
    exit(1);
  }

  /* Determine problem name */
  CUTEST_pname( &status, &funit, pname );
  if ( status ) {
    printf( "cutest_init error: cutest_pname - status not equal to zero. status = %d\n", status );
    exit( status );
  }

  /* Determine problem size */
  CUTEST_cdimen( &status, &funit, &CUTEst_nvar, &CUTEst_ncon );
  if ( status ) {
    printf( "cutest_init error: cutest_cdimen - status not equal to zero. status = %d\n", status );
    exit( status );
  }

  if ( CUTEst_nvar != CUTEst_ncon ) {
    printf( "cutest_init error: This program tackles nonlinear systems with n=m only!\n" );
    exit( status );
  }

  /* Reserve memory for variables, bounds, and multipliers */
  // R_x0   = PROTECT( allocVector( REALSXP, CUTEst_nvar ) );
  // R_crhs = PROTECT( allocVector( REALSXP, CUTEst_ncon+1 ) );
  // MALLOC( crhs,   CUTEst_ncon+1, doublereal );
  MALLOC( bl,     CUTEst_nvar, doublereal );
  MALLOC( bu,     CUTEst_nvar, doublereal );
  MALLOC( cl,     CUTEst_ncon+1, doublereal );
  MALLOC( equatn, CUTEst_ncon+1, logical    );
  MALLOC( lambda, CUTEst_ncon+1, doublereal );
  MALLOC( linear, CUTEst_ncon+1, logical    );
  // MALLOC( x0,     CUTEst_nvar, doublereal );

  /* Call initialization routine for CUTEst */
  CUTEST_csetup( &status, &funit, &iout, &io_buffer, &CUTEst_nvar, &CUTEst_ncon, REAL( R_x0 ), bl, bu,
		 lambda, cl, REAL( R_crhs ), equatn, linear, &e_order, &l_order, &v_order );
  if ( status ) {
    printf( "cutest_init error: cutest_csetup - status not equal to zero. status = %d\n", status );
    exit( status );
  }

  for ( i = 0; i < CUTEst_ncon; i++ ) {
    if ( ! equatn[i] ) {
      printf( "cutest_init error: The problem %s has inequalities!\n", pname );
      exit( 1 );
    }

    if ( cl[i] != REAL(R_crhs)[i] ) {
      printf( "cutest_init error: The problem %s has an equality if cl not equal to cu!\n", pname );
      exit( 1 );
    }
  }

  for ( i = 0; i < CUTEst_nvar; i++ )
    if ( bl[i] > CUTE_INF || bu[i] < CUTE_INF ) {
      printf( "cutest_init error: The problem %s has bounds!\n", pname );
      exit( 1 );
    }

  /* printf( "Problem name = %s\nx = \n", pname ); */
  /* for ( i = 0; i < CUTEst_nvar; i++ ) */
  /*   printf ( "%d %lf\n", i, x0[i] ); */

  /* Close problem description file OUTSDIF.d */
  FORTRAN_close( &funit, &ierr );

  /* Free unneeded arrays */
  FREE( bl );
  FREE( bu );
  FREE( cl );
  FREE( equatn );
  FREE( lambda );
  FREE( linear );

  return R_NilValue;
}

SEXP c_cutest_evalr (SEXP R_nvar, SEXP R_x, SEXP R_crhs) {
  const integer lcjac=1;
  integer i, indfun[1], indvar[1], nnzj, status;
  doublereal cjac[1];
  const logical grad = FALSE_;

  SEXP R_r = PROTECT( allocVector( REALSXP, INTEGER(R_nvar)[0] ) );

  CUTEST_ccfsg( &status, INTEGER(R_nvar), INTEGER(R_nvar), REAL(R_x),
		REAL(R_r), &nnzj, &lcjac, cjac, indvar, indfun,
		&grad );
  if ( status ) {
    printf( "c_cutest_evalr error: There was a nonnull flag in a call to CUTEst routine cutest_ccfsg\n" );
    exit(1);
  }

  for ( i = 0; i < INTEGER(R_nvar)[0]; i++ )
    REAL(R_r)[i] -= REAL(R_crhs)[i];

  UNPROTECT(1);

  return R_r;
}

SEXP c_cutest_end () {
  integer status;

  CUTEST_cterminate( &status );
  return R_NilValue;
}
