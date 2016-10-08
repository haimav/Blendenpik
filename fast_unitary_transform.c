#include "mex.h"

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <sys/timeb.h>

#include "config.h"

#ifdef USE_FFTW
#include "fftw_r2r.h"
#endif

#ifdef USE_SWHT
#include "spiral_wht.h"
#endif

#ifdef USE_FFTW
int dct_sizes[FFTW_TIMES];
int dht_sizes[FFTW_TIMES];

void load_fftw_wisdom()
{
  FILE *wisdom_file;
  char *wisdom_string;
  long size, start;

  wisdom_file = fopen(FFTW_WISDOM_FILE, "r");
  if (wisdom_file == NULL)
  	mexErrMsgTxt("Could not find wisdom file");
  fread(dht_sizes, sizeof(int), FFTW_TIMES, wisdom_file);
  fread(dct_sizes, sizeof(int), FFTW_TIMES, wisdom_file);
  start = ftell(wisdom_file);
  fseek(wisdom_file, 0, SEEK_END);
  size = ftell(wisdom_file) - start;
  fseek(wisdom_file, start, SEEK_SET);
  wisdom_string = malloc(sizeof(char) * size);
  fread(wisdom_string, 1, size, wisdom_file);
  fclose(wisdom_file);

  fftw_import_wisdom_from_string(wisdom_string);
}
#endif

void mexFunction(int nargout, mxArray *argout[], int nargin, const mxArray *argin[])
{
  int n, m0, mm, m;
  double *X, *Y, *D;
  int i, j;
  unsigned kind;
  char type[5];
  double scale = 0.0;
  int wht_tree_size = 0;  /* only for SPIRAL WHT */
  int inverse = 0;
  #ifdef USE_FFTW
  int wis_level = FFTW_ESTIMATE;
  #endif

  m0 = mxGetM(argin[0]);
  n = mxGetN(argin[0]);
  X = mxGetPr(argin[0]);
  D = mxGetPr(argin[1]);

  if (mxGetM(argin[1]) * mxGetN(argin[1]) < m0)
	  mexErrMsgTxt("Scale vector too small!");

  mxGetString(argin[2], type, 5);

#ifdef USE_FFTW
  if (strcmp(type, "DCT") == 0) {
	inverse = 0;
    if (nargin < 4) {
      load_fftw_wisdom();
      wis_level = FFTW_WISDOM_FLAG;
      mm = (m0 % FFTW_QUANT == 0) ? m0 : m0 + (FFTW_QUANT - (m0 % FFTW_QUANT));
      m = dct_sizes[mm / FFTW_QUANT - 1];
    } else
      m = m0;
    kind = FFTW_R2R_DCT;
    scale = 1 / sqrt(2 * m);

  }

  if (strcmp(type, "DHT") == 0) {
	inverse = 0;
    if (nargin < 4) {
      load_fftw_wisdom();
      wis_level = FFTW_WISDOM_FLAG;
      mm = (m0 % FFTW_QUANT == 0) ? m0 : m0 + (FFTW_QUANT - (m0 % FFTW_QUANT));
      m = dht_sizes[mm / FFTW_QUANT - 1];
    } else
      m = m0;
    kind = FFTW_R2R_DHT;
    scale = 1 / sqrt(m);
  }

    if (strcmp(type, "IDCT") == 0) {
      inverse = 1;
      if (nargin < 4) {
        load_fftw_wisdom();
        wis_level = FFTW_WISDOM_FLAG;
        mm = (m0 % FFTW_QUANT == 0) ? m0 : m0 + (FFTW_QUANT - (m0 % FFTW_QUANT));
        m = dct_sizes[mm / FFTW_QUANT - 1];
      } else
        m = m0;
      kind = FFTW_R2R_IDCT;
      scale = 1 / sqrt(2 * m);

    }

    if (strcmp(type, "IDHT") == 0) {
	  inverse = 1;
      if (nargin < 4) {
        load_fftw_wisdom();
        wis_level = FFTW_WISDOM_FLAG;
        mm = (m0 % FFTW_QUANT == 0) ? m0 : m0 + (FFTW_QUANT - (m0 % FFTW_QUANT));
        m = dht_sizes[mm / FFTW_QUANT - 1];
      } else
        m = m0;
      kind = FFTW_R2R_DHT;
      scale = 1 / sqrt(m);
  }
#endif

  /* if (strcmp(type, "WHS") == 0) {
    wht_tree_size = 0;
    kind = FFTW_R2R_NONE;
    m = (int)pow(2, ceil(log(m0) / log(2)));
    scale = 1 / sqrt(m);
  }*/

#ifdef USE_SWHT
  if (strcmp(type, "WHT") == 0) {
	inverse = 0;
    wht_tree_size = ceil(log(m0) / log(2)); /* --> SPIRAL WHT */
    kind = FFTW_R2R_NONE;
    m = (int)pow(2, ceil(log(m0) / log(2)));
    scale = 1 / sqrt(m);
  }

  if (strcmp(type, "IWHT") == 0) {
	  inverse = 1;
      wht_tree_size = ceil(log(m0) / log(2)); /* --> SPIRAL WHT */
      kind = FFTW_R2R_NONE;
      m = (int)pow(2, ceil(log(m0) / log(2)));
      scale = 1 / sqrt(m);
  }
#endif

  if (scale == 0.0)
    mexErrMsgTxt("Unrecognized transformation type!");

  argout[0] = mxCreateDoubleMatrix(m, n, mxREAL);
  Y = mxGetPr(argout[0]);
  if (!inverse) {
	  for(i = 0; i < n; i++) {
		for(j = 0; j < m0; j++)
		  Y[m * i + j] = scale * D[j] * X[m0 * i + j];
		if (m0 != m)
		  memset(Y + m * i + m0, 0, (m - m0) * sizeof(double));
	  }
  } else {
	   for(i = 0; i < n; i++) {
		  memcpy(Y + m * i, X + m0 * i, m0 * sizeof(double));
	      if (m0 != m)
		  	memset(Y + m * i + m0, 0, (m - m0) * sizeof(double));
       }
  }

#ifdef USE_FFTW
    /* For IDCT need to rescale first row */
    if (kind == FFTW_R2R_IDCT)
      for (i = 0; i < n; i++)
        Y[i * m] *= sqrt(2);
#endif

#ifdef USE_FFTW
  if (kind != FFTW_R2R_NONE) {
    int r;

    r = fftw_r2r(Y, Y, m, n, kind, wis_level);
    if (!r)
      mexErrMsgTxt("FFTW wisdom not loaded");
  } else {
	  #endif
	  #ifdef USE_SWHT
    /*if (wht_tree_size == 0)
      WHT(Y, m, n);
      else { */
    Wht *wht_tree = wht_get_tree(wht_tree_size);
    for(i = 0; i < n; i++)
      wht_apply(wht_tree, 1, Y + i * m);
    wht_delete(wht_tree);
    /*}*/
    #endif
  }

#ifdef USE_FFTW
  /* For DCT need to rescale first row */
  if (kind == FFTW_R2R_DCT)
    for (i = 0; i < n; i++)
      Y[i * m] /= sqrt(2);
#endif

  if (inverse) {
	  for(i = 0; i < n; i++)
		for(j = 0; j < m; j++)
		  Y[m * i + j] = scale * D[j] * Y[m * i + j];
  }
}
