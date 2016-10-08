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
  int m, m0, mm;
  char type[5];
  #ifdef USE_FFTW
  int wis_level = FFTW_ESTIMATE;
  #endif


  m0 = (int)mxGetScalar(argin[0]);
  mxGetString(argin[1], type, 5);

#ifdef USE_FFTW
  if (strcmp(type, "DCT") == 0) {
      load_fftw_wisdom();
      mm = (m0 % FFTW_QUANT == 0) ? m0 : m0 + (FFTW_QUANT - (m0 % FFTW_QUANT));
      m = dct_sizes[mm / FFTW_QUANT - 1];
  }

  if (strcmp(type, "DHT") == 0) {
      load_fftw_wisdom();
      mm = (m0 % FFTW_QUANT == 0) ? m0 : m0 + (FFTW_QUANT - (m0 % FFTW_QUANT));
      m = dht_sizes[mm / FFTW_QUANT - 1];
  }

    if (strcmp(type, "IDCT") == 0) {
        load_fftw_wisdom();
        mm = (m0 % FFTW_QUANT == 0) ? m0 : m0 + (FFTW_QUANT - (m0 % FFTW_QUANT));
        m = dct_sizes[mm / FFTW_QUANT - 1];
    }

    if (strcmp(type, "IDHT") == 0) {
        load_fftw_wisdom();
        mm = (m0 % FFTW_QUANT == 0) ? m0 : m0 + (FFTW_QUANT - (m0 % FFTW_QUANT));
        m = dht_sizes[mm / FFTW_QUANT - 1];
  }
#endif

  /* if (strcmp(type, "WHS") == 0) {
    wht_tree_size = 0;
    kind = FFTW_R2R_NONE;
    m = (int)pow(2, ceil(log(m0) / log(2)));
    scale = 1 / sqrt(m);
  }*/

	argout[0] = mxCreateDoubleScalar(m);
}
