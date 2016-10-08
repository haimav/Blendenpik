#include "mex.h"

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <sys/timeb.h>

#include "fftw_r2r.h"


void mexFunction(int nargout, mxArray *argout[], int nargin, const mxArray *argin[]) 
{
  int n, m;
  double *X, *Y;
  unsigned kind = FFTW_R2R_NONE;
  char type[4];
  
  if (nargin == 1)
	  mexErrMsgTxt("Type not specified");

  mxGetString(argin[1], type, 4);
  if (strcmp(type, "DCT") == 0)
	kind = FFTW_R2R_DCT;
  if (strcmp(type, "DHT") == 0)
	kind = FFTW_R2R_DHT;

  if (kind == FFTW_R2R_NONE)
	  mexErrMsgTxt("Type not recognized");
      
  m = mxGetM(argin[0]);
  n = mxGetN(argin[0]);
  X = mxGetPr(argin[0]);
  
  argout[0] = mxCreateDoubleMatrix(m, n, mxREAL);
  Y = mxGetPr(argout[0]);
      
  fftw_r2r(X, Y, m, n, kind, FFTW_ESTIMATE);
}
