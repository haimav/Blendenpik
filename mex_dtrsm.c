#include "mex.h"

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <sys/timeb.h>

#include "config.h"

#ifdef BLAS_UNDERSCORE
#define DTRSM dtrsm_
#else
#define DTRSM dtrsm
#endif


double wtime()
{
  struct timeb T;

  static time_t time_diff;
  static time_t mill_diff;

  double dt;

  (void) ftime( &T );

  time_diff = T.time;
  mill_diff = T.millitm;

  dt = ((double) time_diff) + (1e-3) * ((double) mill_diff);

  return dt;
}

void mexFunction(int nargout, mxArray *argout[], int nargin, const mxArray *argin[])
{
  long m, nrhs, lda;
  double *A, *b, alpha;
  char side[2], uplo[2], trans[2], diag[2];

  A = mxGetPr(argin[0]);
  lda = mxGetM(argin[0]);
  argout[0] = mxDuplicateArray(argin[1]);
  b = mxGetPr(argout[0]);
  m = mxGetM(argout[0]);
  nrhs = mxGetN(argout[0]);
  alpha = mxGetScalar(argin[2]);
  mxGetString(argin[3], side, 2);
  mxGetString(argin[4], uplo, 2);
  mxGetString(argin[5], trans, 2);
  mxGetString(argin[6], diag, 2);

  DTRSM(side, uplo, trans, diag, &m, &nrhs, &alpha, A, &lda, b, &m);
}
