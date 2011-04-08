#include "rb_lapack.h"

extern VOID dlaed9_(integer *k, integer *kstart, integer *kstop, integer *n, doublereal *d, doublereal *q, integer *ldq, doublereal *rho, doublereal *dlamda, doublereal *w, doublereal *s, integer *lds, integer *info);

static VALUE sHelp, sUsage;

static VALUE
rblapack_dlaed9(int argc, VALUE *argv, VALUE self){
  VALUE rblapack_kstart;
  integer kstart; 
  VALUE rblapack_kstop;
  integer kstop; 
  VALUE rblapack_n;
  integer n; 
  VALUE rblapack_rho;
  doublereal rho; 
  VALUE rblapack_dlamda;
  doublereal *dlamda; 
  VALUE rblapack_w;
  doublereal *w; 
  VALUE rblapack_d;
  doublereal *d; 
  VALUE rblapack_s;
  doublereal *s; 
  VALUE rblapack_info;
  integer info; 
  doublereal *q;

  integer k;
  integer lds;
  integer ldq;

  VALUE rb_options;
  if (argc > 0 && TYPE(argv[argc-1]) == T_HASH) {
    argc--;
    rb_options = argv[argc];
    if (rb_hash_aref(rb_options, sHelp) == Qtrue) {
      printf("%s\n", "USAGE:\n  d, s, info = NumRu::Lapack.dlaed9( kstart, kstop, n, rho, dlamda, w, [:usage => usage, :help => help])\n\n\nFORTRAN MANUAL\n      SUBROUTINE DLAED9( K, KSTART, KSTOP, N, D, Q, LDQ, RHO, DLAMDA, W, S, LDS, INFO )\n\n*  Purpose\n*  =======\n*\n*  DLAED9 finds the roots of the secular equation, as defined by the\n*  values in D, Z, and RHO, between KSTART and KSTOP.  It makes the\n*  appropriate calls to DLAED4 and then stores the new matrix of\n*  eigenvectors for use in calculating the next level of Z vectors.\n*\n\n*  Arguments\n*  =========\n*\n*  K       (input) INTEGER\n*          The number of terms in the rational function to be solved by\n*          DLAED4.  K >= 0.\n*\n*  KSTART  (input) INTEGER\n*  KSTOP   (input) INTEGER\n*          The updated eigenvalues Lambda(I), KSTART <= I <= KSTOP\n*          are to be computed.  1 <= KSTART <= KSTOP <= K.\n*\n*  N       (input) INTEGER\n*          The number of rows and columns in the Q matrix.\n*          N >= K (delation may result in N > K).\n*\n*  D       (output) DOUBLE PRECISION array, dimension (N)\n*          D(I) contains the updated eigenvalues\n*          for KSTART <= I <= KSTOP.\n*\n*  Q       (workspace) DOUBLE PRECISION array, dimension (LDQ,N)\n*\n*  LDQ     (input) INTEGER\n*          The leading dimension of the array Q.  LDQ >= max( 1, N ).\n*\n*  RHO     (input) DOUBLE PRECISION\n*          The value of the parameter in the rank one update equation.\n*          RHO >= 0 required.\n*\n*  DLAMDA  (input) DOUBLE PRECISION array, dimension (K)\n*          The first K elements of this array contain the old roots\n*          of the deflated updating problem.  These are the poles\n*          of the secular equation.\n*\n*  W       (input) DOUBLE PRECISION array, dimension (K)\n*          The first K elements of this array contain the components\n*          of the deflation-adjusted updating vector.\n*\n*  S       (output) DOUBLE PRECISION array, dimension (LDS, K)\n*          Will contain the eigenvectors of the repaired matrix which\n*          will be stored for subsequent Z vector calculation and\n*          multiplied by the previously accumulated eigenvectors\n*          to update the system.\n*\n*  LDS     (input) INTEGER\n*          The leading dimension of S.  LDS >= max( 1, K ).\n*\n*  INFO    (output) INTEGER\n*          = 0:  successful exit.\n*          < 0:  if INFO = -i, the i-th argument had an illegal value.\n*          > 0:  if INFO = 1, an eigenvalue did not converge\n*\n\n*  Further Details\n*  ===============\n*\n*  Based on contributions by\n*     Jeff Rutter, Computer Science Division, University of California\n*     at Berkeley, USA\n*\n*  =====================================================================\n*\n*     .. Local Scalars ..\n      INTEGER            I, J\n      DOUBLE PRECISION   TEMP\n*     ..\n*     .. External Functions ..\n      DOUBLE PRECISION   DLAMC3, DNRM2\n      EXTERNAL           DLAMC3, DNRM2\n*     ..\n*     .. External Subroutines ..\n      EXTERNAL           DCOPY, DLAED4, XERBLA\n*     ..\n*     .. Intrinsic Functions ..\n      INTRINSIC          MAX, SIGN, SQRT\n*     ..\n\n");
      rb_exit(0);
    }
    if (rb_hash_aref(rb_options, sUsage) == Qtrue) {
      printf("%s\n", "USAGE:\n  d, s, info = NumRu::Lapack.dlaed9( kstart, kstop, n, rho, dlamda, w, [:usage => usage, :help => help])\n");
      rb_exit(0);
    } 
  } else
    rb_options = Qnil;
  if (argc != 6)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 6)", argc);
  rblapack_kstart = argv[0];
  rblapack_kstop = argv[1];
  rblapack_n = argv[2];
  rblapack_rho = argv[3];
  rblapack_dlamda = argv[4];
  rblapack_w = argv[5];
  if (rb_options != Qnil) {
  }

  if (!NA_IsNArray(rblapack_w))
    rb_raise(rb_eArgError, "w (6th argument) must be NArray");
  if (NA_RANK(rblapack_w) != 1)
    rb_raise(rb_eArgError, "rank of w (6th argument) must be %d", 1);
  k = NA_SHAPE0(rblapack_w);
  if (NA_TYPE(rblapack_w) != NA_DFLOAT)
    rblapack_w = na_change_type(rblapack_w, NA_DFLOAT);
  w = NA_PTR_TYPE(rblapack_w, doublereal*);
  kstart = NUM2INT(rblapack_kstart);
  if (!NA_IsNArray(rblapack_dlamda))
    rb_raise(rb_eArgError, "dlamda (5th argument) must be NArray");
  if (NA_RANK(rblapack_dlamda) != 1)
    rb_raise(rb_eArgError, "rank of dlamda (5th argument) must be %d", 1);
  if (NA_SHAPE0(rblapack_dlamda) != k)
    rb_raise(rb_eRuntimeError, "shape 0 of dlamda must be the same as shape 0 of w");
  if (NA_TYPE(rblapack_dlamda) != NA_DFLOAT)
    rblapack_dlamda = na_change_type(rblapack_dlamda, NA_DFLOAT);
  dlamda = NA_PTR_TYPE(rblapack_dlamda, doublereal*);
  rho = NUM2DBL(rblapack_rho);
  kstop = NUM2INT(rblapack_kstop);
  n = NUM2INT(rblapack_n);
  lds = MAX( 1, k );
  ldq = MAX( 1, n );
  {
    int shape[1];
    shape[0] = MAX(1,n);
    rblapack_d = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  d = NA_PTR_TYPE(rblapack_d, doublereal*);
  {
    int shape[2];
    shape[0] = lds;
    shape[1] = k;
    rblapack_s = na_make_object(NA_DFLOAT, 2, shape, cNArray);
  }
  s = NA_PTR_TYPE(rblapack_s, doublereal*);
  q = ALLOC_N(doublereal, (ldq)*(MAX(1,n)));

  dlaed9_(&k, &kstart, &kstop, &n, d, q, &ldq, &rho, dlamda, w, s, &lds, &info);

  free(q);
  rblapack_info = INT2NUM(info);
  return rb_ary_new3(3, rblapack_d, rblapack_s, rblapack_info);
}

void
init_lapack_dlaed9(VALUE mLapack){
  rb_define_module_function(mLapack, "dlaed9", rblapack_dlaed9, -1);
  sHelp = ID2SYM(rb_intern("help"));
  sUsage = ID2SYM(rb_intern("usage"));
}
