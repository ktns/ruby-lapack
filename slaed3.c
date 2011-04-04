#include "rb_lapack.h"

extern VOID slaed3_(integer *k, integer *n, integer *n1, real *d, real *q, integer *ldq, real *rho, real *dlamda, real *q2, integer *indx, integer *ctot, real *w, real *s, integer *info);

static VALUE
rb_slaed3(int argc, VALUE *argv, VALUE self){
  VALUE rb_n1;
  integer n1; 
  VALUE rb_rho;
  real rho; 
  VALUE rb_dlamda;
  real *dlamda; 
  VALUE rb_q2;
  real *q2; 
  VALUE rb_indx;
  integer *indx; 
  VALUE rb_ctot;
  integer *ctot; 
  VALUE rb_w;
  real *w; 
  VALUE rb_d;
  real *d; 
  VALUE rb_q;
  real *q; 
  VALUE rb_info;
  integer info; 
  VALUE rb_dlamda_out__;
  real *dlamda_out__;
  VALUE rb_w_out__;
  real *w_out__;
  real *s;

  integer k;
  integer n;
  integer ldq;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  d, q, info, dlamda, w = NumRu::Lapack.slaed3( n1, rho, dlamda, q2, indx, ctot, w)\n    or\n  NumRu::Lapack.slaed3  # print help\n\n\nFORTRAN MANUAL\n      SUBROUTINE SLAED3( K, N, N1, D, Q, LDQ, RHO, DLAMDA, Q2, INDX, CTOT, W, S, INFO )\n\n*  Purpose\n*  =======\n*\n*  SLAED3 finds the roots of the secular equation, as defined by the\n*  values in D, W, and RHO, between 1 and K.  It makes the\n*  appropriate calls to SLAED4 and then updates the eigenvectors by\n*  multiplying the matrix of eigenvectors of the pair of eigensystems\n*  being combined by the matrix of eigenvectors of the K-by-K system\n*  which is solved here.\n*\n*  This code makes very mild assumptions about floating point\n*  arithmetic. It will work on machines with a guard digit in\n*  add/subtract, or on those binary machines without guard digits\n*  which subtract like the Cray X-MP, Cray Y-MP, Cray C-90, or Cray-2.\n*  It could conceivably fail on hexadecimal or decimal machines\n*  without guard digits, but we know of none.\n*\n\n*  Arguments\n*  =========\n*\n*  K       (input) INTEGER\n*          The number of terms in the rational function to be solved by\n*          SLAED4.  K >= 0.\n*\n*  N       (input) INTEGER\n*          The number of rows and columns in the Q matrix.\n*          N >= K (deflation may result in N>K).\n*\n*  N1      (input) INTEGER\n*          The location of the last eigenvalue in the leading submatrix.\n*          min(1,N) <= N1 <= N/2.\n*\n*  D       (output) REAL array, dimension (N)\n*          D(I) contains the updated eigenvalues for\n*          1 <= I <= K.\n*\n*  Q       (output) REAL array, dimension (LDQ,N)\n*          Initially the first K columns are used as workspace.\n*          On output the columns 1 to K contain\n*          the updated eigenvectors.\n*\n*  LDQ     (input) INTEGER\n*          The leading dimension of the array Q.  LDQ >= max(1,N).\n*\n*  RHO     (input) REAL\n*          The value of the parameter in the rank one update equation.\n*          RHO >= 0 required.\n*\n*  DLAMDA  (input/output) REAL array, dimension (K)\n*          The first K elements of this array contain the old roots\n*          of the deflated updating problem.  These are the poles\n*          of the secular equation. May be changed on output by\n*          having lowest order bit set to zero on Cray X-MP, Cray Y-MP,\n*          Cray-2, or Cray C-90, as described above.\n*\n*  Q2      (input) REAL array, dimension (LDQ2, N)\n*          The first K columns of this matrix contain the non-deflated\n*          eigenvectors for the split problem.\n*\n*  INDX    (input) INTEGER array, dimension (N)\n*          The permutation used to arrange the columns of the deflated\n*          Q matrix into three groups (see SLAED2).\n*          The rows of the eigenvectors found by SLAED4 must be likewise\n*          permuted before the matrix multiply can take place.\n*\n*  CTOT    (input) INTEGER array, dimension (4)\n*          A count of the total number of the various types of columns\n*          in Q, as described in INDX.  The fourth column type is any\n*          column which has been deflated.\n*\n*  W       (input/output) REAL array, dimension (K)\n*          The first K elements of this array contain the components\n*          of the deflation-adjusted updating vector. Destroyed on\n*          output.\n*\n*  S       (workspace) REAL array, dimension (N1 + 1)*K\n*          Will contain the eigenvectors of the repaired matrix which\n*          will be multiplied by the previously accumulated eigenvectors\n*          to update the system.\n*\n*  LDS     (input) INTEGER\n*          The leading dimension of S.  LDS >= max(1,K).\n*\n*  INFO    (output) INTEGER\n*          = 0:  successful exit.\n*          < 0:  if INFO = -i, the i-th argument had an illegal value.\n*          > 0:  if INFO = 1, an eigenvalue did not converge\n*\n\n*  Further Details\n*  ===============\n*\n*  Based on contributions by\n*     Jeff Rutter, Computer Science Division, University of California\n*     at Berkeley, USA\n*  Modified by Francoise Tisseur, University of Tennessee.\n*\n*  =====================================================================\n*\n\n");
    return Qnil;
  }
  if (argc != 7)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 7)", argc);
  rb_n1 = argv[0];
  rb_rho = argv[1];
  rb_dlamda = argv[2];
  rb_q2 = argv[3];
  rb_indx = argv[4];
  rb_ctot = argv[5];
  rb_w = argv[6];

  if (!NA_IsNArray(rb_ctot))
    rb_raise(rb_eArgError, "ctot (6th argument) must be NArray");
  if (NA_RANK(rb_ctot) != 1)
    rb_raise(rb_eArgError, "rank of ctot (6th argument) must be %d", 1);
  if (NA_SHAPE0(rb_ctot) != (4))
    rb_raise(rb_eRuntimeError, "shape 0 of ctot must be %d", 4);
  if (NA_TYPE(rb_ctot) != NA_LINT)
    rb_ctot = na_change_type(rb_ctot, NA_LINT);
  ctot = NA_PTR_TYPE(rb_ctot, integer*);
  if (!NA_IsNArray(rb_w))
    rb_raise(rb_eArgError, "w (7th argument) must be NArray");
  if (NA_RANK(rb_w) != 1)
    rb_raise(rb_eArgError, "rank of w (7th argument) must be %d", 1);
  k = NA_SHAPE0(rb_w);
  if (NA_TYPE(rb_w) != NA_SFLOAT)
    rb_w = na_change_type(rb_w, NA_SFLOAT);
  w = NA_PTR_TYPE(rb_w, real*);
  if (!NA_IsNArray(rb_q2))
    rb_raise(rb_eArgError, "q2 (4th argument) must be NArray");
  if (NA_RANK(rb_q2) != 2)
    rb_raise(rb_eArgError, "rank of q2 (4th argument) must be %d", 2);
  n = NA_SHAPE1(rb_q2);
  if (NA_SHAPE0(rb_q2) != n)
    rb_raise(rb_eRuntimeError, "shape 0 of q2 must be the same as shape 1 of q2");
  if (NA_TYPE(rb_q2) != NA_SFLOAT)
    rb_q2 = na_change_type(rb_q2, NA_SFLOAT);
  q2 = NA_PTR_TYPE(rb_q2, real*);
  if (!NA_IsNArray(rb_dlamda))
    rb_raise(rb_eArgError, "dlamda (3th argument) must be NArray");
  if (NA_RANK(rb_dlamda) != 1)
    rb_raise(rb_eArgError, "rank of dlamda (3th argument) must be %d", 1);
  if (NA_SHAPE0(rb_dlamda) != k)
    rb_raise(rb_eRuntimeError, "shape 0 of dlamda must be the same as shape 0 of w");
  if (NA_TYPE(rb_dlamda) != NA_SFLOAT)
    rb_dlamda = na_change_type(rb_dlamda, NA_SFLOAT);
  dlamda = NA_PTR_TYPE(rb_dlamda, real*);
  if (!NA_IsNArray(rb_indx))
    rb_raise(rb_eArgError, "indx (5th argument) must be NArray");
  if (NA_RANK(rb_indx) != 1)
    rb_raise(rb_eArgError, "rank of indx (5th argument) must be %d", 1);
  if (NA_SHAPE0(rb_indx) != n)
    rb_raise(rb_eRuntimeError, "shape 0 of indx must be the same as shape 1 of q2");
  if (NA_TYPE(rb_indx) != NA_LINT)
    rb_indx = na_change_type(rb_indx, NA_LINT);
  indx = NA_PTR_TYPE(rb_indx, integer*);
  n1 = NUM2INT(rb_n1);
  rho = (real)NUM2DBL(rb_rho);
  ldq = MAX(1,n);
  {
    int shape[1];
    shape[0] = n;
    rb_d = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  d = NA_PTR_TYPE(rb_d, real*);
  {
    int shape[2];
    shape[0] = ldq;
    shape[1] = n;
    rb_q = na_make_object(NA_SFLOAT, 2, shape, cNArray);
  }
  q = NA_PTR_TYPE(rb_q, real*);
  {
    int shape[1];
    shape[0] = k;
    rb_dlamda_out__ = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  dlamda_out__ = NA_PTR_TYPE(rb_dlamda_out__, real*);
  MEMCPY(dlamda_out__, dlamda, real, NA_TOTAL(rb_dlamda));
  rb_dlamda = rb_dlamda_out__;
  dlamda = dlamda_out__;
  {
    int shape[1];
    shape[0] = k;
    rb_w_out__ = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  w_out__ = NA_PTR_TYPE(rb_w_out__, real*);
  MEMCPY(w_out__, w, real, NA_TOTAL(rb_w));
  rb_w = rb_w_out__;
  w = w_out__;
  s = ALLOC_N(real, (MAX(1,k))*(n1 + 1));

  slaed3_(&k, &n, &n1, d, q, &ldq, &rho, dlamda, q2, indx, ctot, w, s, &info);

  free(s);
  rb_info = INT2NUM(info);
  return rb_ary_new3(5, rb_d, rb_q, rb_info, rb_dlamda, rb_w);
}

void
init_lapack_slaed3(VALUE mLapack){
  rb_define_module_function(mLapack, "slaed3", rb_slaed3, -1);
}
