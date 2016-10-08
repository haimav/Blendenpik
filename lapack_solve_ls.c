#include "mex.h"

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <sys/timeb.h>

#include "config.h"

#ifdef BLAS_UNDERSCORE
#define DGELS dgels_
#else
#define DGELS dgels
#endif

#define MAX(x, y) ((x) > (y) ? (x) : (y))

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
  long m, n;
  double *x, *b, *A;
  long workspace_size, info;
  double *workspace;
  double wsize_d;
  long nrhs;
  long ldb;
  double dt;

  m = mxGetM(argin[0]);
  n = mxGetN(argin[0]);
  nrhs = mxGetN(argin[1]);

  argout[0] = mxCreateDoubleMatrix(n, nrhs, mxREAL);
  x = mxGetPr(argout[0]);


  /* Copy input to avoid corruption */
  A = (double *)mxMalloc(m * n * sizeof(double));
  memcpy(A, mxGetPr(argin[0]), m * n * sizeof(double));
  ldb = MAX(m, n);
  b = (double *)mxMalloc(ldb * sizeof(double));
  memcpy(b, mxGetPr(argin[1]), m * sizeof(double));

  /* Query workspace size */
  workspace_size = -1;
  workspace = &wsize_d;
  DGELS("None", &m, &n, &nrhs, A, &m, b, &ldb, workspace, &workspace_size, &info);

  /* Compute */
  workspace_size = (long)wsize_d;
  workspace = (double *)mxMalloc(sizeof(double) * workspace_size);
  dt = -wtime();
  DGELS("None", &m, &n, &nrhs, A, &m, b, &ldb, workspace, &workspace_size, &info);
  dt += wtime();
  mexPrintf("LAPACK time is %f\n", dt);
  memcpy(x, b, n * sizeof(double));
  argout[1] = mxCreateDoubleScalar(dt);
  mxFree(workspace);
  mxFree(A);
  mxFree(b);
}
