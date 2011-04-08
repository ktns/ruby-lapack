#include "rb_lapack.h"

extern VOID slaein_(logical *rightv, logical *noinit, integer *n, real *h, integer *ldh, real *wr, real *wi, real *vr, real *vi, real *b, integer *ldb, real *work, real *eps3, real *smlnum, real *bignum, integer *info);

static VALUE sHelp, sUsage;

static VALUE
rblapack_slaein(int argc, VALUE *argv, VALUE self){
  VALUE rblapack_rightv;
  logical rightv; 
  VALUE rblapack_noinit;
  logical noinit; 
  VALUE rblapack_h;
  real *h; 
  VALUE rblapack_wr;
  real wr; 
  VALUE rblapack_wi;
  real wi; 
  VALUE rblapack_vr;
  real *vr; 
  VALUE rblapack_vi;
  real *vi; 
  VALUE rblapack_eps3;
  real eps3; 
  VALUE rblapack_smlnum;
  real smlnum; 
  VALUE rblapack_bignum;
  real bignum; 
  VALUE rblapack_info;
  integer info; 
  VALUE rblapack_vr_out__;
  real *vr_out__;
  VALUE rblapack_vi_out__;
  real *vi_out__;
  real *b;
  real *work;

  integer ldh;
  integer n;
  integer ldb;

  VALUE rb_options;
  if (argc > 0 && TYPE(argv[argc-1]) == T_HASH) {
    argc--;
    rb_options = argv[argc];
    if (rb_hash_aref(rb_options, sHelp) == Qtrue) {
      printf("%s\n", "USAGE:\n  info, vr, vi = NumRu::Lapack.slaein( rightv, noinit, h, wr, wi, vr, vi, eps3, smlnum, bignum, [:usage => usage, :help => help])\n\n\nFORTRAN MANUAL\n      SUBROUTINE SLAEIN( RIGHTV, NOINIT, N, H, LDH, WR, WI, VR, VI, B, LDB, WORK, EPS3, SMLNUM, BIGNUM, INFO )\n\n*  Purpose\n*  =======\n*\n*  SLAEIN uses inverse iteration to find a right or left eigenvector\n*  corresponding to the eigenvalue (WR,WI) of a real upper Hessenberg\n*  matrix H.\n*\n\n*  Arguments\n*  =========\n*\n*  RIGHTV   (input) LOGICAL\n*          = .TRUE. : compute right eigenvector;\n*          = .FALSE.: compute left eigenvector.\n*\n*  NOINIT   (input) LOGICAL\n*          = .TRUE. : no initial vector supplied in (VR,VI).\n*          = .FALSE.: initial vector supplied in (VR,VI).\n*\n*  N       (input) INTEGER\n*          The order of the matrix H.  N >= 0.\n*\n*  H       (input) REAL array, dimension (LDH,N)\n*          The upper Hessenberg matrix H.\n*\n*  LDH     (input) INTEGER\n*          The leading dimension of the array H.  LDH >= max(1,N).\n*\n*  WR      (input) REAL\n*  WI      (input) REAL\n*          The real and imaginary parts of the eigenvalue of H whose\n*          corresponding right or left eigenvector is to be computed.\n*\n*  VR      (input/output) REAL array, dimension (N)\n*  VI      (input/output) REAL array, dimension (N)\n*          On entry, if NOINIT = .FALSE. and WI = 0.0, VR must contain\n*          a real starting vector for inverse iteration using the real\n*          eigenvalue WR; if NOINIT = .FALSE. and WI.ne.0.0, VR and VI\n*          must contain the real and imaginary parts of a complex\n*          starting vector for inverse iteration using the complex\n*          eigenvalue (WR,WI); otherwise VR and VI need not be set.\n*          On exit, if WI = 0.0 (real eigenvalue), VR contains the\n*          computed real eigenvector; if WI.ne.0.0 (complex eigenvalue),\n*          VR and VI contain the real and imaginary parts of the\n*          computed complex eigenvector. The eigenvector is normalized\n*          so that the component of largest magnitude has magnitude 1;\n*          here the magnitude of a complex number (x,y) is taken to be\n*          |x| + |y|.\n*          VI is not referenced if WI = 0.0.\n*\n*  B       (workspace) REAL array, dimension (LDB,N)\n*\n*  LDB     (input) INTEGER\n*          The leading dimension of the array B.  LDB >= N+1.\n*\n*  WORK   (workspace) REAL array, dimension (N)\n*\n*  EPS3    (input) REAL\n*          A small machine-dependent value which is used to perturb\n*          close eigenvalues, and to replace zero pivots.\n*\n*  SMLNUM  (input) REAL\n*          A machine-dependent value close to the underflow threshold.\n*\n*  BIGNUM  (input) REAL\n*          A machine-dependent value close to the overflow threshold.\n*\n*  INFO    (output) INTEGER\n*          = 0:  successful exit\n*          = 1:  inverse iteration did not converge; VR is set to the\n*                last iterate, and so is VI if WI.ne.0.0.\n*\n\n*  =====================================================================\n*\n\n");
      rb_exit(0);
    }
    if (rb_hash_aref(rb_options, sUsage) == Qtrue) {
      printf("%s\n", "USAGE:\n  info, vr, vi = NumRu::Lapack.slaein( rightv, noinit, h, wr, wi, vr, vi, eps3, smlnum, bignum, [:usage => usage, :help => help])\n");
      rb_exit(0);
    } 
  } else
    rb_options = Qnil;
  if (argc != 10)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 10)", argc);
  rblapack_rightv = argv[0];
  rblapack_noinit = argv[1];
  rblapack_h = argv[2];
  rblapack_wr = argv[3];
  rblapack_wi = argv[4];
  rblapack_vr = argv[5];
  rblapack_vi = argv[6];
  rblapack_eps3 = argv[7];
  rblapack_smlnum = argv[8];
  rblapack_bignum = argv[9];
  if (rb_options != Qnil) {
  }

  smlnum = (real)NUM2DBL(rblapack_smlnum);
  eps3 = (real)NUM2DBL(rblapack_eps3);
  wr = (real)NUM2DBL(rblapack_wr);
  if (!NA_IsNArray(rblapack_vr))
    rb_raise(rb_eArgError, "vr (6th argument) must be NArray");
  if (NA_RANK(rblapack_vr) != 1)
    rb_raise(rb_eArgError, "rank of vr (6th argument) must be %d", 1);
  n = NA_SHAPE0(rblapack_vr);
  if (NA_TYPE(rblapack_vr) != NA_SFLOAT)
    rblapack_vr = na_change_type(rblapack_vr, NA_SFLOAT);
  vr = NA_PTR_TYPE(rblapack_vr, real*);
  rightv = (rblapack_rightv == Qtrue);
  noinit = (rblapack_noinit == Qtrue);
  bignum = (real)NUM2DBL(rblapack_bignum);
  if (!NA_IsNArray(rblapack_vi))
    rb_raise(rb_eArgError, "vi (7th argument) must be NArray");
  if (NA_RANK(rblapack_vi) != 1)
    rb_raise(rb_eArgError, "rank of vi (7th argument) must be %d", 1);
  if (NA_SHAPE0(rblapack_vi) != n)
    rb_raise(rb_eRuntimeError, "shape 0 of vi must be the same as shape 0 of vr");
  if (NA_TYPE(rblapack_vi) != NA_SFLOAT)
    rblapack_vi = na_change_type(rblapack_vi, NA_SFLOAT);
  vi = NA_PTR_TYPE(rblapack_vi, real*);
  if (!NA_IsNArray(rblapack_h))
    rb_raise(rb_eArgError, "h (3th argument) must be NArray");
  if (NA_RANK(rblapack_h) != 2)
    rb_raise(rb_eArgError, "rank of h (3th argument) must be %d", 2);
  if (NA_SHAPE1(rblapack_h) != n)
    rb_raise(rb_eRuntimeError, "shape 1 of h must be the same as shape 0 of vr");
  ldh = NA_SHAPE0(rblapack_h);
  if (NA_TYPE(rblapack_h) != NA_SFLOAT)
    rblapack_h = na_change_type(rblapack_h, NA_SFLOAT);
  h = NA_PTR_TYPE(rblapack_h, real*);
  wi = (real)NUM2DBL(rblapack_wi);
  ldb = n+1;
  {
    int shape[1];
    shape[0] = n;
    rblapack_vr_out__ = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  vr_out__ = NA_PTR_TYPE(rblapack_vr_out__, real*);
  MEMCPY(vr_out__, vr, real, NA_TOTAL(rblapack_vr));
  rblapack_vr = rblapack_vr_out__;
  vr = vr_out__;
  {
    int shape[1];
    shape[0] = n;
    rblapack_vi_out__ = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  vi_out__ = NA_PTR_TYPE(rblapack_vi_out__, real*);
  MEMCPY(vi_out__, vi, real, NA_TOTAL(rblapack_vi));
  rblapack_vi = rblapack_vi_out__;
  vi = vi_out__;
  b = ALLOC_N(real, (ldb)*(n));
  work = ALLOC_N(real, (n));

  slaein_(&rightv, &noinit, &n, h, &ldh, &wr, &wi, vr, vi, b, &ldb, work, &eps3, &smlnum, &bignum, &info);

  free(b);
  free(work);
  rblapack_info = INT2NUM(info);
  return rb_ary_new3(3, rblapack_info, rblapack_vr, rblapack_vi);
}

void
init_lapack_slaein(VALUE mLapack){
  rb_define_module_function(mLapack, "slaein", rblapack_slaein, -1);
  sHelp = ID2SYM(rb_intern("help"));
  sUsage = ID2SYM(rb_intern("usage"));
}
