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
#define DDOT  ddot_
#define DAXPY daxpy_
#else
#define DGEMV dgemv
#define DTRSV dtrsv
#define DNRM2 dnrm2
#define DDOT  ddot
#define DAXPY daxpy
#endif


double DNRM2(long *, double *, long *);
double DDOT(long *, double *, long *, double *, long *);

#define MAX_FULLIT 80

void scale(double *x, int n, double s)
{
  int i;

  for(i = 0; i < n; i++)
    x[i] *= s;
}

void mexFunction(int nargout, mxArray *argout[], int nargin, const mxArray *argin[]) 
{
  double tol;
  double *A, *b, *R, *x, *U, *u, *V, *v, *vt, *d, *z;
  double beta, alpha, normr, normar, norma;
  double c, s, phibar, phi, nn;
  double thet, rhot, rho;

  double *lsvec, *resvec, *Atr, *r;

  long m, n;
  long maxit, it, i, k;

  long int_zero = 0;
  long int_one = 1;
  double dbl_one = 1.0;
  double dbl_mone = -1.0;
  double dbl_zero = 0.0;


  A = mxGetPr(argin[0]);
  b = mxGetPr(argin[1]);
  R = mxIsEmpty(argin[2]) ? NULL : mxGetPr(argin[2]);

  tol = mxGetScalar(argin[3]);
  maxit = (int)mxGetScalar(argin[4]);

  m = mxGetM(argin[0]);
  n = mxGetN(argin[0]);

  U = malloc(MAX_FULLIT * m * sizeof(double));
  u = malloc(m * sizeof(double));
  V = malloc(MAX_FULLIT * n * sizeof(double));
  v = malloc(n * sizeof(double));
  vt = malloc(n * sizeof(double));
  d = malloc(n * sizeof(double));
  z = malloc(n * sizeof(double));

  argout[0] = mxCreateDoubleMatrix(n, 1, mxREAL);
  x = mxGetPr(argout[0]);

  if (nargout > 2) {
    argout[2] = mxCreateDoubleMatrix(maxit+1, 1, mxREAL);
    lsvec = mxGetPr(argout[2]);
    argout[3] = mxCreateDoubleMatrix(maxit+1, 1, mxREAL);
    resvec = mxGetPr(argout[3]);

    Atr = malloc(n * sizeof(double));
    r = malloc(m * sizeof(double));
    memcpy(r, b, m * sizeof(double));

    DGEMV("T", &m, &n, &dbl_one, A, &m, r, &int_one, &dbl_zero, Atr, &int_one);
    resvec[0] = DNRM2(&m, r, &int_one);
    lsvec[0] = DNRM2(&n, Atr, &int_one);
  }

  memset(x, 0, n * sizeof(double));
  memset(d, 0, n * sizeof(double));
  memcpy(u, b, m * sizeof(double));
  beta = DNRM2(&m, u, &int_one);
  normr = beta;
  scale(u, m, 1/beta);
  memcpy(U, u, m * sizeof(double));

  c = 1; s = 0; phibar = beta;
  DGEMV("T", &m, &n, &dbl_one, A, &m, u, &int_one, &dbl_zero, v, &int_one);
  if (R != NULL)
    DTRSV("U", "T", "Not Unit", &n, R, &n, v, &int_one);

  alpha = DNRM2(&n, v, &int_one);
  scale(v, n, 1/alpha);
  memcpy(V, v, n * sizeof(double));
	
  normar = alpha * beta;
  norma = 0;

  it = 0;
  while (it < maxit) {
    double malpha = -alpha;
    memcpy(z, v, n * sizeof(double));
    if (R != NULL)
      DTRSV("U", "N", "Not Unit", &n, R, &n, z, &int_one);
    DGEMV("N", &m, &n, &dbl_one, A, &m, z, &int_one, &malpha, u, &int_one);

    if (it < MAX_FULLIT)
      for(k = 0; k <= it; k++) {
        double eta = -DDOT(&m, u, &int_one, U + k * m, &int_one);
        DAXPY(&m, &eta, U + k * m, &int_one, u, &int_one);
      }
	
    
    beta = DNRM2(&m, u, &int_one);
    scale(u, m, 1/beta);
    if (it < MAX_FULLIT-1)
      memcpy(U + (it + 1) * m, u, m * sizeof(double));

    norma = sqrt(norma * norma + alpha * alpha + beta * beta);

    thet = - s * alpha;
    rhot = c * alpha;
    rho = sqrt(rhot * rhot + beta * beta);
    c = rhot / rho;
    s = - beta / rho;
    phi = c * phibar;
    phibar = s * phibar;
		
    for (i = 0; i < n; i++) {      
      d[i] = (z[i] - thet * d[i]) / rho;
      x[i] = x[i] + phi * d[i];
    }
    it++;

    if (nargout > 2) {
      memcpy(r, b, m * sizeof(double));
      DGEMV("N", &m, &n,  &dbl_mone, A, &m, x, &int_one, &dbl_one, r, &int_one);
      DGEMV("T", &m, &n, &dbl_one, A, &m, r, &int_one, &dbl_zero, Atr, &int_one);
      lsvec[it] = DNRM2(&n, Atr, &int_one);
      resvec[it] = DNRM2(&m, r, &int_one);
    }

    normr = fabs(s) * normr;
    nn = normar / (normr * norma);

    if (nn < tol)
      break;

    DGEMV("T", &m, &n, &dbl_one, A, &m, u, &int_one, &dbl_zero, vt, &int_one);
    if (R != NULL)
      DTRSV("U", "T", "Not Unit", &n, R, &n, vt, &int_one);
    for (i = 0; i < n; i++)
      v[i] = vt[i] - beta * v[i];

    if (it < MAX_FULLIT)
      for(k = 0; k < it; k++) {
        double eta = -DDOT(&n, v, &int_one, V + k * n, &int_one);
        DAXPY(&n, &eta, V + k * n, &int_one, v, &int_one);
      }
    alpha = DNRM2(&n, v, &int_one);
    scale(v, n, 1/alpha);	
    if (it < MAX_FULLIT)
      memcpy(V + it * n, v, n * sizeof(double));

    normar = alpha * fabs(s * phi);
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
  free(v);
  free(vt);
  free(d);
  free(z);
  if (nargout > 2) {
    free(Atr);
    free(r);
  }
}
