#include "mex.h"

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <sys/timeb.h>

#include "config.h"

#ifdef BLAS_UNDERSCORE
#define DORMQR dormqr_
#else
#define DORMQR dormqr
#endif

#define MIN(x, y) ((x) < (y) ? (x) : (y))

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
  long m, n, k;
  double *tau, *A, *C;
  long workspace_size, info;
  double *workspace;
  double wsize_d;
  char side[2], trans[2];

  A = mxGetPr(argin[0]);
  k = MIN(mxGetM(argin[0]), mxGetN(argin[0]));
  tau = mxGetPr(argin[1]);
  m = mxGetM(argin[2]);
  n = mxGetN(argin[2]);
  argout[0] = mxDuplicateArray(argin[2]);
  C = mxGetPr(argout[0]);
  mxGetString(argin[3], side, 2);
  mxGetString(argin[4], trans, 2);

  /* Query workspace size */
  workspace_size = -1;
  workspace = &wsize_d;
  DORMQR(side, trans, &m, &n, &k, A, &m, tau, C, &m, workspace, &workspace_size, &info);

  /* Compute */
  workspace_size = (int)wsize_d;
  workspace = (double *)mxMalloc(sizeof(double) * workspace_size);
  DORMQR(side, trans, &m, &n, &k, A, &m, tau, C, &m, workspace, &workspace_size, &info);
  mxFree(workspace);
}
