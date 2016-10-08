#include "mex.h"

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <sys/timeb.h>

#include "config.h"
#include "fftw3.h"


double read_wtime() 
{
  struct timeb T;
  (void) ftime( &T );
  
  return ((double)T.time) + (1e-3) * ((double)T.millitm);
}
void mexFunction(int nargout, mxArray *argout[], int nargin, const mxArray *argin[]) 
{
  int n, m, i;
  mxArray *XX;
  double *X;
  fftw_plan plan;
  char type[4];
  double t0, tt0, elp, t_all0;
  char wisdom_file_name[255];
  FILE *wisdom_file;
  unsigned level_flag;
  int level;
  int tm;

  int dct_sizes[FFTW_TIMES];
  double dct_times[FFTW_TIMES];
  int dht_sizes[FFTW_TIMES];
  double dht_times[FFTW_TIMES];

  t_all0 = read_wtime();

  level = (int)(mxGetScalar(argin[0]));
  switch(level) {
  case 0:
    level_flag = FFTW_ESTIMATE;
    break;

  case 1:
    level_flag = FFTW_MEASURE;
    break;

  case 2:
    level_flag = FFTW_PATIENT;
    break;
			
  case 3:
    level_flag = FFTW_EXHAUSTIVE;
    break;
  }

  n = 10;


  for(m = FFTW_TIMES * FFTW_QUANT; m >= FFTW_QUANT; m -= FFTW_QUANT) {
    tm = m / FFTW_QUANT - 1;
    t0 = read_wtime();
    XX = mxCreateDoubleMatrix(m, n, mxREAL);
    X = mxGetPr(XX);
    memset(X, 1, m * n * sizeof(double));
    plan = fftw_plan_r2r_1d(m, X, X, FFTW_DHT, FFTW_UNALIGNED | level_flag);
    fftw_execute_r2r(plan, X, X); /* execute one time for estimate */
    fftw_destroy_plan(plan);
    plan = fftw_plan_r2r_1d(m, X, X, FFTW_DHT, FFTW_UNALIGNED | level_flag); /* replan for estimate */
    tt0 = read_wtime();
    for(i = 0; i < n; i++)
      fftw_execute_r2r(plan, X + i * m, X + i * m);
    elp = read_wtime() - tt0;
    dht_times[tm] = elp;
    dht_sizes[tm] = m;
    fftw_destroy_plan(plan);
    mxDestroyArray(XX);

    if (tm < FFTW_TIMES-1 && dht_times[tm] > dht_times[tm + 1]) {
      dht_times[tm] = dht_times[tm + 1];
      dht_sizes[tm] = dht_sizes[tm + 1];
    } 
    mexPrintf("DHT size %d->%d run time %.2e planning time is %.2e sec\n", m, dht_sizes[tm], elp, read_wtime() - t0);
    mexEvalString("drawnow;");
  }

  for(m = FFTW_TIMES * FFTW_QUANT; m >= FFTW_QUANT; m -= FFTW_QUANT) {
    tm = m / FFTW_QUANT - 1;
    t0 = read_wtime();
    XX = mxCreateDoubleMatrix(m, n, mxREAL);
    X = mxGetPr(XX);
    memset(X, 1, m * n * sizeof(double));
    plan = fftw_plan_r2r_1d(m, X, X, FFTW_REDFT10, FFTW_UNALIGNED | level_flag);
    fftw_execute_r2r(plan, X, X); /* execute one time for estimate */
    fftw_destroy_plan(plan);
    plan = fftw_plan_r2r_1d(m, X, X, FFTW_REDFT10, FFTW_UNALIGNED | level_flag);
    tt0 = read_wtime();
    for(i = 0; i < n; i++)
      fftw_execute_r2r(plan, X + i * m, X + i * m);
    elp = read_wtime() - tt0;
    dct_times[tm] = elp;
    dct_sizes[tm] = m;
    fftw_destroy_plan(plan);
    mxDestroyArray(XX);

    if (tm < FFTW_TIMES-1 && dct_times[tm] > dct_times[tm + 1]) {
      dct_times[tm] = dct_times[tm + 1];
      dct_sizes[tm] = dct_sizes[tm + 1];
    } 
    mexPrintf("DCT size %d->%d run time %.2e planning time is %.2e sec\n", m, dct_sizes[tm], elp, read_wtime() - t0);
    mexEvalString("drawnow;");
  }

  /* Write wisdom to file. FFTW's function doesn't work for some reason... */
  mxGetString(argin[1], wisdom_file_name, 255);	
  wisdom_file = fopen(wisdom_file_name, "w+");
  if (wisdom_file == NULL)
    mexErrMsgTxt("Cannot create wisdom file"); 
  fwrite(dht_sizes, sizeof(int), FFTW_TIMES, wisdom_file);
  fwrite(dct_sizes, sizeof(int), FFTW_TIMES, wisdom_file);
  fputs(fftw_export_wisdom_to_string(), wisdom_file);
  fclose(wisdom_file);

  mexPrintf("TOTAL time is %.2e sec\n", read_wtime() - t_all0); 
}
