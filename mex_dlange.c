#include "mex.h"

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <sys/timeb.h>

#include "config.h"

#ifdef BLAS_UNDERSCORE
#define DLANGE dlange_
#else
#define DLANGE dlange
#endif

double DLANGE(char *, long *, long *, double *, long *, double *);

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
  long n, ld; 
  double norm;
  double *A, dummy;

  
  ld = mxGetM(argin[0]);
  n = mxGetN(argin[0]);
  A = mxGetPr(argin[0]);
  
  /* Compute */
  norm = DLANGE("F", &ld, &n, A, &ld, &dummy);
  argout[0] = mxCreateDoubleScalar(norm);
}
