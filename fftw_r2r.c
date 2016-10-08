#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <sys/timeb.h>

#include "fftw_r2r.h"

int fftw_r2r(double *X, double *Y, int m, int n, unsigned type, unsigned level)
{
  fftw_plan plan;
  fftw_r2r_kind kind;
  int i, j;
  unsigned flags = FFTW_UNALIGNED;

  switch(type) {
	case FFTW_R2R_DCT:
		kind = FFTW_REDFT10;
		break;

	case FFTW_R2R_IDCT:
		kind = FFTW_REDFT01;
		break;

    case FFTW_R2R_DHT:
        kind = FFTW_DHT;
		break;

	default:
            return(0);
  }

  if (X != Y)
	  flags |= FFTW_PRESERVE_INPUT;

  plan = fftw_plan_r2r_1d(m, Y, Y, kind, flags | level);
  if (plan == NULL)
	  return(0);
  for (j = 0; j < n; j++)
	  fftw_execute_r2r(plan, X + j * m, Y + j * m);
  fftw_destroy_plan(plan);

  return(1);
}
