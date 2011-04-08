#include "rb_lapack.h"

extern VOID dlaed3_(integer *k, integer *n, integer *n1, doublereal *d, doublereal *q, integer *ldq, doublereal *rho, doublereal *dlamda, doublereal *q2, integer *indx, integer *ctot, doublereal *w, doublereal *s, integer *info);

static VALUE sHelp, sUsage;

static VALUE
rblapack_dlaed3(int argc, VALUE *argv, VALUE self){
  VALUE rblapack_n1;
  integer n1; 
  VALUE rblapack_rho;
  doublereal rho; 
  VALUE rblapack_dlamda;
  doublereal *dlamda; 
  VALUE rblapack_q2;
  doublereal *q2; 
  VALUE rblapack_indx;
  integer *indx; 
  VALUE rblapack_ctot;
  integer *ctot; 
  VALUE rblapack_w;
  doublereal *w; 
  VALUE rblapack_d;
  doublereal *d; 
  VALUE rblapack_q;
  doublereal *q; 
  VALUE rblapack_info;
  integer info; 
  VALUE rblapack_dlamda_out__;
  doublereal *dlamda_out__;
  VALUE rblapack_w_out__;
  doublereal *w_out__;
  doublereal *s;

  integer k;
  integer n;
  integer ldq;

  VALUE rb_options;
  if (argc > 0 && TYPE(argv[argc-1]) == T_HASH) {
    argc--;
    rb_options = argv[argc];
    if (rb_hash_aref(rb_options, sHelp) == Qtrue) {
      printf("%s\n", "USAGE:\n  d, q, info, dlamda, w = NumRu::Lapack.dlaed3( n1, rho, dlamda, q2, indx, ctot, w, [:usage => usage, :help => help])\n\n\nFORTRAN MANUAL\n      SUBROUTINE DLAED3( K, N, N1, D, Q, LDQ, RHO, DLAMDA, Q2, INDX, CTOT, W, S, INFO )\n\n*  Purpose\n*  =======\n*\n*  DLAED3 finds the roots of the secular equation, as defined by the\n*  values in D, W, and RHO, between 1 and K.  It makes the\n*  appropriate calls to DLAED4 and then updates the eigenvectors by\n*  multiplying the matrix of eigenvectors of the pair of eigensystems\n*  being combined by the matrix of eigenvectors of the K-by-K system\n*  which is solved here.\n*\n*  This code makes very mild assumptions about floating point\n*  arithmetic. It will work on machines with a guard digit in\n*  add/subtract, or on those binary machines without guard digits\n*  which subtract like the Cray X-MP, Cray Y-MP, Cray C-90, or Cray-2.\n*  It could conceivably fail on hexadecimal or decimal machines\n*  without guard digits, but we know of none.\n*\n\n*  Arguments\n*  =========\n*\n*  K       (input) INTEGER\n*          The number of terms in the rational function to be solved by\n*          DLAED4.  K >= 0.\n*\n*  N       (input) INTEGER\n*          The number of rows and columns in the Q matrix.\n*          N >= K (deflation may result in N>K).\n*\n*  N1      (input) INTEGER\n*          The location of the last eigenvalue in the leading submatrix.\n*          min(1,N) <= N1 <= N/2.\n*\n*  D       (output) DOUBLE PRECISION array, dimension (N)\n*          D(I) contains the updated eigenvalues for\n*          1 <= I <= K.\n*\n*  Q       (output) DOUBLE PRECISION array, dimension (LDQ,N)\n*          Initially the first K columns are used as workspace.\n*          On output the columns 1 to K contain\n*          the updated eigenvectors.\n*\n*  LDQ     (input) INTEGER\n*          The leading dimension of the array Q.  LDQ >= max(1,N).\n*\n*  RHO     (input) DOUBLE PRECISION\n*          The value of the parameter in the rank one update equation.\n*          RHO >= 0 required.\n*\n*  DLAMDA  (input/output) DOUBLE PRECISION array, dimension (K)\n*          The first K elements of this array contain the old roots\n*          of the deflated updating problem.  These are the poles\n*          of the secular equation. May be changed on output by\n*          having lowest order bit set to zero on Cray X-MP, Cray Y-MP,\n*          Cray-2, or Cray C-90, as described above.\n*\n*  Q2      (input) DOUBLE PRECISION array, dimension (LDQ2, N)\n*          The first K columns of this matrix contain the non-deflated\n*          eigenvectors for the split problem.\n*\n*  INDX    (input) INTEGER array, dimension (N)\n*          The permutation used to arrange the columns of the deflated\n*          Q matrix into three groups (see DLAED2).\n*          The rows of the eigenvectors found by DLAED4 must be likewise\n*          permuted before the matrix multiply can take place.\n*\n*  CTOT    (input) INTEGER array, dimension (4)\n*          A count of the total number of the various types of columns\n*          in Q, as described in INDX.  The fourth column type is any\n*          column which has been deflated.\n*\n*  W       (input/output) DOUBLE PRECISION array, dimension (K)\n*          The first K elements of this array contain the components\n*          of the deflation-adjusted updating vector. Destroyed on\n*          output.\n*\n*  S       (workspace) DOUBLE PRECISION array, dimension (N1 + 1)*K\n*          Will contain the eigenvectors of the repaired matrix which\n*          will be multiplied by the previously accumulated eigenvectors\n*          to update the system.\n*\n*  LDS     (input) INTEGER\n*          The leading dimension of S.  LDS >= max(1,K).\n*\n*  INFO    (output) INTEGER\n*          = 0:  successful exit.\n*          < 0:  if INFO = -i, the i-th argument had an illegal value.\n*          > 0:  if INFO = 1, an eigenvalue did not converge\n*\n\n*  Further Details\n*  ===============\n*\n*  Based on contributions by\n*     Jeff Rutter, Computer Science Division, University of California\n*     at Berkeley, USA\n*  Modified by Francoise Tisseur, University of Tennessee.\n*\n*  =====================================================================\n*\n\n");
      rb_exit(0);
    }
    if (rb_hash_aref(rb_options, sUsage) == Qtrue) {
      printf("%s\n", "USAGE:\n  d, q, info, dlamda, w = NumRu::Lapack.dlaed3( n1, rho, dlamda, q2, indx, ctot, w, [:usage => usage, :help => help])\n");
      rb_exit(0);
    } 
  } else
    rb_options = Qnil;
  if (argc != 7)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 7)", argc);
  rblapack_n1 = argv[0];
  rblapack_rho = argv[1];
  rblapack_dlamda = argv[2];
  rblapack_q2 = argv[3];
  rblapack_indx = argv[4];
  rblapack_ctot = argv[5];
  rblapack_w = argv[6];
  if (rb_options != Qnil) {
  }

  if (!NA_IsNArray(rblapack_ctot))
    rb_raise(rb_eArgError, "ctot (6th argument) must be NArray");
  if (NA_RANK(rblapack_ctot) != 1)
    rb_raise(rb_eArgError, "rank of ctot (6th argument) must be %d", 1);
  if (NA_SHAPE0(rblapack_ctot) != (4))
    rb_raise(rb_eRuntimeError, "shape 0 of ctot must be %d", 4);
  if (NA_TYPE(rblapack_ctot) != NA_LINT)
    rblapack_ctot = na_change_type(rblapack_ctot, NA_LINT);
  ctot = NA_PTR_TYPE(rblapack_ctot, integer*);
  if (!NA_IsNArray(rblapack_w))
    rb_raise(rb_eArgError, "w (7th argument) must be NArray");
  if (NA_RANK(rblapack_w) != 1)
    rb_raise(rb_eArgError, "rank of w (7th argument) must be %d", 1);
  k = NA_SHAPE0(rblapack_w);
  if (NA_TYPE(rblapack_w) != NA_DFLOAT)
    rblapack_w = na_change_type(rblapack_w, NA_DFLOAT);
  w = NA_PTR_TYPE(rblapack_w, doublereal*);
  if (!NA_IsNArray(rblapack_q2))
    rb_raise(rb_eArgError, "q2 (4th argument) must be NArray");
  if (NA_RANK(rblapack_q2) != 2)
    rb_raise(rb_eArgError, "rank of q2 (4th argument) must be %d", 2);
  n = NA_SHAPE1(rblapack_q2);
  if (NA_SHAPE0(rblapack_q2) != n)
    rb_raise(rb_eRuntimeError, "shape 0 of q2 must be the same as shape 1 of q2");
  if (NA_TYPE(rblapack_q2) != NA_DFLOAT)
    rblapack_q2 = na_change_type(rblapack_q2, NA_DFLOAT);
  q2 = NA_PTR_TYPE(rblapack_q2, doublereal*);
  if (!NA_IsNArray(rblapack_dlamda))
    rb_raise(rb_eArgError, "dlamda (3th argument) must be NArray");
  if (NA_RANK(rblapack_dlamda) != 1)
    rb_raise(rb_eArgError, "rank of dlamda (3th argument) must be %d", 1);
  if (NA_SHAPE0(rblapack_dlamda) != k)
    rb_raise(rb_eRuntimeError, "shape 0 of dlamda must be the same as shape 0 of w");
  if (NA_TYPE(rblapack_dlamda) != NA_DFLOAT)
    rblapack_dlamda = na_change_type(rblapack_dlamda, NA_DFLOAT);
  dlamda = NA_PTR_TYPE(rblapack_dlamda, doublereal*);
  if (!NA_IsNArray(rblapack_indx))
    rb_raise(rb_eArgError, "indx (5th argument) must be NArray");
  if (NA_RANK(rblapack_indx) != 1)
    rb_raise(rb_eArgError, "rank of indx (5th argument) must be %d", 1);
  if (NA_SHAPE0(rblapack_indx) != n)
    rb_raise(rb_eRuntimeError, "shape 0 of indx must be the same as shape 1 of q2");
  if (NA_TYPE(rblapack_indx) != NA_LINT)
    rblapack_indx = na_change_type(rblapack_indx, NA_LINT);
  indx = NA_PTR_TYPE(rblapack_indx, integer*);
  n1 = NUM2INT(rblapack_n1);
  rho = NUM2DBL(rblapack_rho);
  ldq = MAX(1,n);
  {
    int shape[1];
    shape[0] = n;
    rblapack_d = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  d = NA_PTR_TYPE(rblapack_d, doublereal*);
  {
    int shape[2];
    shape[0] = ldq;
    shape[1] = n;
    rblapack_q = na_make_object(NA_DFLOAT, 2, shape, cNArray);
  }
  q = NA_PTR_TYPE(rblapack_q, doublereal*);
  {
    int shape[1];
    shape[0] = k;
    rblapack_dlamda_out__ = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  dlamda_out__ = NA_PTR_TYPE(rblapack_dlamda_out__, doublereal*);
  MEMCPY(dlamda_out__, dlamda, doublereal, NA_TOTAL(rblapack_dlamda));
  rblapack_dlamda = rblapack_dlamda_out__;
  dlamda = dlamda_out__;
  {
    int shape[1];
    shape[0] = k;
    rblapack_w_out__ = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  w_out__ = NA_PTR_TYPE(rblapack_w_out__, doublereal*);
  MEMCPY(w_out__, w, doublereal, NA_TOTAL(rblapack_w));
  rblapack_w = rblapack_w_out__;
  w = w_out__;
  s = ALLOC_N(doublereal, (MAX(1,k))*((n1 + 1)));

  dlaed3_(&k, &n, &n1, d, q, &ldq, &rho, dlamda, q2, indx, ctot, w, s, &info);

  free(s);
  rblapack_info = INT2NUM(info);
  return rb_ary_new3(5, rblapack_d, rblapack_q, rblapack_info, rblapack_dlamda, rblapack_w);
}

void
init_lapack_dlaed3(VALUE mLapack){
  rb_define_module_function(mLapack, "dlaed3", rblapack_dlaed3, -1);
  sHelp = ID2SYM(rb_intern("help"));
  sUsage = ID2SYM(rb_intern("usage"));
}
