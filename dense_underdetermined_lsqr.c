#include "mex.h"

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

#include "config.h"

#ifdef BLAS_UNDERSCORE
#define DGEMV dgemv_
#define DTRSV dtrsv_
#define DNRM2 dnrm2_
#else
#define DGEMV dgemv
#define DTRSV dtrsv
#define DNRM2 dnrm2
#endif


double DNRM2(long *, double *, long *);


void scale(double *x, int n, double s)
{
  int i;

  for(i = 0; i < n; i++)
    x[i] *= s;
}

void mexFunction(int nargout, mxArray *argout[], int nargin, const mxArray *argin[])
{
  double tol;
  double *A, *b, *L, *x, *u, *ut, *v, *vt, *d, *z;
  double beta, alpha, normr;
  double c, s, phibar, phi, nn;
  double thet, rhot, rho;
  double mbeta;

  double *resvec, *xvec, *Atr, *r;

  long m, n;
  long maxit, it, i;

  long int_zero = 0;
  long int_one = 1;
  double dbl_one = 1.0;
  double dbl_mone = -1.0;
  double dbl_zero = 0.0;

  A = mxGetPr(argin[0]);
  b = mxGetPr(argin[1]);
  L = mxIsEmpty(argin[2]) ? NULL : mxGetPr(argin[2]);

  m = mxGetM(argin[0]);
  n = mxGetN(argin[0]);

  tol = mxGetScalar(argin[3]) * DNRM2(&m, b, &int_one);
  maxit = (int)mxGetScalar(argin[4]);

  u = malloc(m * sizeof(double));
  ut = malloc(m * sizeof(double));
  v = malloc(n * sizeof(double));
  vt = malloc(n * sizeof(double));
  d = malloc(n * sizeof(double));
  z = malloc(n * sizeof(double));

  argout[0] = mxCreateDoubleMatrix(n, 1, mxREAL);
  x = mxGetPr(argout[0]);

  if (nargout > 2) {
    argout[2] = mxCreateDoubleMatrix(maxit+1, 1, mxREAL);
    resvec = mxGetPr(argout[2]);
    argout[3] = mxCreateDoubleMatrix(maxit+1, 1, mxREAL);
    xvec = mxGetPr(argout[3]);

    r = malloc(m * sizeof(double));
    memcpy(r, b, m * sizeof(double));

    resvec[0] = DNRM2(&m, r, &int_one);
    xvec[0] = DNRM2(&n, x, &int_one);
  }

  memset(x, 0, n * sizeof(double));
  memset(d, 0, n * sizeof(double));

  memcpy(u, b, m * sizeof(double));
  if (L != NULL)
  	DTRSV("L", "N", "Not Unit", &m, L, &m, u, &int_one);

  beta = DNRM2(&m, u, &int_one);
  normr = beta;
  scale(u, m, 1/beta);
  c = 1; s = 0; phibar = beta;

  memcpy(z, u, m * sizeof(double));
  if (L != NULL)
    DTRSV("L", "T", "Not Unit", &m, L, &m, z, &int_one); 
  DGEMV("T", &m, &n, &dbl_one, A, &m, z, &int_one, &dbl_zero, v, &int_one);

  alpha = DNRM2(&n, v, &int_one);
  scale(v, n, 1/alpha);

  it = 0;
  while (it < maxit) {

    DGEMV("N", &m, &n, &dbl_one, A, &m, v, &int_one, &dbl_zero, ut, &int_one);
    if (L != NULL)
      DTRSV("L", "N", "Not Unit", &m, L, &m, ut, &int_one);
    for (i = 0; i < m; i++)
      u[i] = ut[i] - alpha * u[i];

    beta = DNRM2(&m, u, &int_one);
    scale(u, m, 1/beta);

    thet = - s * alpha;
    rhot = c * alpha;
    rho = sqrt(rhot * rhot + beta * beta);
    c = rhot / rho;
    s = - beta / rho;
    phi = c * phibar;
    phibar = s * phibar;
		
    for (i = 0; i < n; i++) {      
      d[i] = (v[i] - thet * d[i]) / rho;
      x[i] = x[i] + phi * d[i];
    }
    it++;

    if (nargout > 2) {
      memcpy(r, b, m * sizeof(double));
      DGEMV("N", &m, &n,  &dbl_mone, A, &m, x, &int_one, &dbl_one, r, &int_one);
      resvec[it] = DNRM2(&m, r, &int_one);
      xvec[it] = DNRM2(&n, x, &int_one);
    }

    normr = fabs(s) * normr;
    if (normr < tol)
      break;

    mbeta = -beta;
    memcpy(z, u, m * sizeof(double));
    if (L != NULL)
      DTRSV("L", "T", "Not Unit", &m, L, &m, z, &int_one);
    DGEMV("T", &m, &n, &dbl_one, A, &m, z, &int_one, &mbeta, v, &int_one);


    alpha = DNRM2(&n, v, &int_one);
    scale(v, n, 1/alpha);	
  }

  if (nargout > 2){
    mxSetM(argout[2] , it + 1);
    mxSetM(argout[3] , it + 1);
  }

  if (nargout > 1)
    argout[1] = mxCreateDoubleScalar(it);
  if (nn > tol)
    mexPrintf("dense_lsqr: did not converge\n");
  else
    mexPrintf("dense_lsqr: converged at iteration %d\n", it);

  free(u);
  free(ut);
  free(v);
  free(d);
  free(z);
  if (nargout > 2)
    free(r);
}
