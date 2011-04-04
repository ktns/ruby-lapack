#include "rb_lapack.h"

extern VOID slaed9_(integer *k, integer *kstart, integer *kstop, integer *n, real *d, real *q, integer *ldq, real *rho, real *dlamda, real *w, real *s, integer *lds, integer *info);

static VALUE
rb_slaed9(int argc, VALUE *argv, VALUE self){
  VALUE rb_kstart;
  integer kstart; 
  VALUE rb_kstop;
  integer kstop; 
  VALUE rb_n;
  integer n; 
  VALUE rb_rho;
  real rho; 
  VALUE rb_dlamda;
  real *dlamda; 
  VALUE rb_w;
  real *w; 
  VALUE rb_d;
  real *d; 
  VALUE rb_s;
  real *s; 
  VALUE rb_info;
  integer info; 
  real *q;

  integer k;
  integer lds;
  integer ldq;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  d, s, info = NumRu::Lapack.slaed9( kstart, kstop, n, rho, dlamda, w)\n    or\n  NumRu::Lapack.slaed9  # print help\n\n\nFORTRAN MANUAL\n      SUBROUTINE SLAED9( K, KSTART, KSTOP, N, D, Q, LDQ, RHO, DLAMDA, W, S, LDS, INFO )\n\n*  Purpose\n*  =======\n*\n*  SLAED9 finds the roots of the secular equation, as defined by the\n*  values in D, Z, and RHO, between KSTART and KSTOP.  It makes the\n*  appropriate calls to SLAED4 and then stores the new matrix of\n*  eigenvectors for use in calculating the next level of Z vectors.\n*\n\n*  Arguments\n*  =========\n*\n*  K       (input) INTEGER\n*          The number of terms in the rational function to be solved by\n*          SLAED4.  K >= 0.\n*\n*  KSTART  (input) INTEGER\n*  KSTOP   (input) INTEGER\n*          The updated eigenvalues Lambda(I), KSTART <= I <= KSTOP\n*          are to be computed.  1 <= KSTART <= KSTOP <= K.\n*\n*  N       (input) INTEGER\n*          The number of rows and columns in the Q matrix.\n*          N >= K (delation may result in N > K).\n*\n*  D       (output) REAL array, dimension (N)\n*          D(I) contains the updated eigenvalues\n*          for KSTART <= I <= KSTOP.\n*\n*  Q       (workspace) REAL array, dimension (LDQ,N)\n*\n*  LDQ     (input) INTEGER\n*          The leading dimension of the array Q.  LDQ >= max( 1, N ).\n*\n*  RHO     (input) REAL\n*          The value of the parameter in the rank one update equation.\n*          RHO >= 0 required.\n*\n*  DLAMDA  (input) REAL array, dimension (K)\n*          The first K elements of this array contain the old roots\n*          of the deflated updating problem.  These are the poles\n*          of the secular equation.\n*\n*  W       (input) REAL array, dimension (K)\n*          The first K elements of this array contain the components\n*          of the deflation-adjusted updating vector.\n*\n*  S       (output) REAL array, dimension (LDS, K)\n*          Will contain the eigenvectors of the repaired matrix which\n*          will be stored for subsequent Z vector calculation and\n*          multiplied by the previously accumulated eigenvectors\n*          to update the system.\n*\n*  LDS     (input) INTEGER\n*          The leading dimension of S.  LDS >= max( 1, K ).\n*\n*  INFO    (output) INTEGER\n*          = 0:  successful exit.\n*          < 0:  if INFO = -i, the i-th argument had an illegal value.\n*          > 0:  if INFO = 1, an eigenvalue did not converge\n*\n\n*  Further Details\n*  ===============\n*\n*  Based on contributions by\n*     Jeff Rutter, Computer Science Division, University of California\n*     at Berkeley, USA\n*\n*  =====================================================================\n*\n*     .. Local Scalars ..\n      INTEGER            I, J\n      REAL               TEMP\n*     ..\n*     .. External Functions ..\n      REAL               SLAMC3, SNRM2\n      EXTERNAL           SLAMC3, SNRM2\n*     ..\n*     .. External Subroutines ..\n      EXTERNAL           SCOPY, SLAED4, XERBLA\n*     ..\n*     .. Intrinsic Functions ..\n      INTRINSIC          MAX, SIGN, SQRT\n*     ..\n\n");
    return Qnil;
  }
  if (argc != 6)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 6)", argc);
  rb_kstart = argv[0];
  rb_kstop = argv[1];
  rb_n = argv[2];
  rb_rho = argv[3];
  rb_dlamda = argv[4];
  rb_w = argv[5];

  if (!NA_IsNArray(rb_w))
    rb_raise(rb_eArgError, "w (6th argument) must be NArray");
  if (NA_RANK(rb_w) != 1)
    rb_raise(rb_eArgError, "rank of w (6th argument) must be %d", 1);
  k = NA_SHAPE0(rb_w);
  if (NA_TYPE(rb_w) != NA_SFLOAT)
    rb_w = na_change_type(rb_w, NA_SFLOAT);
  w = NA_PTR_TYPE(rb_w, real*);
  kstart = NUM2INT(rb_kstart);
  if (!NA_IsNArray(rb_dlamda))
    rb_raise(rb_eArgError, "dlamda (5th argument) must be NArray");
  if (NA_RANK(rb_dlamda) != 1)
    rb_raise(rb_eArgError, "rank of dlamda (5th argument) must be %d", 1);
  if (NA_SHAPE0(rb_dlamda) != k)
    rb_raise(rb_eRuntimeError, "shape 0 of dlamda must be the same as shape 0 of w");
  if (NA_TYPE(rb_dlamda) != NA_SFLOAT)
    rb_dlamda = na_change_type(rb_dlamda, NA_SFLOAT);
  dlamda = NA_PTR_TYPE(rb_dlamda, real*);
  rho = (real)NUM2DBL(rb_rho);
  kstop = NUM2INT(rb_kstop);
  n = NUM2INT(rb_n);
  lds = MAX( 1, k );
  ldq = MAX( 1, n );
  {
    int shape[1];
    shape[0] = MAX(1,n);
    rb_d = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  d = NA_PTR_TYPE(rb_d, real*);
  {
    int shape[2];
    shape[0] = lds;
    shape[1] = k;
    rb_s = na_make_object(NA_SFLOAT, 2, shape, cNArray);
  }
  s = NA_PTR_TYPE(rb_s, real*);
  q = ALLOC_N(real, (ldq)*(MAX(1,n)));

  slaed9_(&k, &kstart, &kstop, &n, d, q, &ldq, &rho, dlamda, w, s, &lds, &info);

  free(q);
  rb_info = INT2NUM(info);
  return rb_ary_new3(3, rb_d, rb_s, rb_info);
}

void
init_lapack_slaed9(VALUE mLapack){
  rb_define_module_function(mLapack, "slaed9", rb_slaed9, -1);
}
