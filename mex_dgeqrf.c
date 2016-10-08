#include "mex.h"

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <sys/timeb.h>

#include "config.h"

#ifdef BLAS_UNDERSCORE
#define DGEQRF dgeqrf_
#else
#define DGEQRF dgeqrf
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
  long m, n;
  double *tau, *A;
  long workspace_size, info;
  double *workspace;
  double wsize_d;

  m = mxGetM(argin[0]);
  n = mxGetN(argin[0]);

  argout[0] = mxDuplicateArray(argin[0]);
  A = mxGetPr(argout[0]);
  if (nargout > 1) {
	  argout[1] = mxCreateDoubleMatrix(n, 1, mxREAL);
	  tau = mxGetPr(argout[1]);
  } else
	  tau = mxMalloc(n * sizeof(double));

  /* Query workspace size */
  workspace_size = -1;
  workspace = &wsize_d;
  DGEQRF(&m, &n, A, &m, tau, workspace, &workspace_size, &info);

  /* Compute */
  workspace_size = (long)wsize_d;
  workspace = (double *)mxMalloc(sizeof(double) * workspace_size);
  DGEQRF(&m, &n, A, &m, tau, workspace, &workspace_size, &info);
  mxFree(workspace);
  if (nargout < 1)
	  mxFree(tau);
}
