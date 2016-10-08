#include "mex.h"

#include <sys/timeb.h>

double read_wtime()
{
  struct timeb T;
  (void) ftime( &T );

  return ((double)T.time) + (1e-3) * ((double)T.millitm);
}

/* Actual mex function */
void mexFunction(int nargout, mxArray *argout[], int nargin, const mxArray *argin[])
{
  argout[0] = mxCreateDoubleScalar(read_wtime());
}
