#include "fftw3.h"

#define FFTW_R2R_NONE  0
#define FFTW_R2R_DCT   1
#define FFTW_R2R_DHT   2
#define FFTW_R2R_IDCT  3

int fftw_r2r(double *X, double *Y, int m, int n, unsigned type, unsigned level);