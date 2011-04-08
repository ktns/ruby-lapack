#include "rb_lapack.h"

extern VOID sggbak_(char *job, char *side, integer *n, integer *ilo, integer *ihi, real *lscale, real *rscale, integer *m, real *v, integer *ldv, integer *info);

static VALUE sHelp, sUsage;

static VALUE
rblapack_sggbak(int argc, VALUE *argv, VALUE self){
  VALUE rblapack_job;
  char job; 
  VALUE rblapack_side;
  char side; 
  VALUE rblapack_ilo;
  integer ilo; 
  VALUE rblapack_ihi;
  integer ihi; 
  VALUE rblapack_lscale;
  real *lscale; 
  VALUE rblapack_rscale;
  real *rscale; 
  VALUE rblapack_v;
  real *v; 
  VALUE rblapack_info;
  integer info; 
  VALUE rblapack_v_out__;
  real *v_out__;

  integer n;
  integer ldv;
  integer m;

  VALUE rb_options;
  if (argc > 0 && TYPE(argv[argc-1]) == T_HASH) {
    argc--;
    rb_options = argv[argc];
    if (rb_hash_aref(rb_options, sHelp) == Qtrue) {
      printf("%s\n", "USAGE:\n  info, v = NumRu::Lapack.sggbak( job, side, ilo, ihi, lscale, rscale, v, [:usage => usage, :help => help])\n\n\nFORTRAN MANUAL\n      SUBROUTINE SGGBAK( JOB, SIDE, N, ILO, IHI, LSCALE, RSCALE, M, V, LDV, INFO )\n\n*  Purpose\n*  =======\n*\n*  SGGBAK forms the right or left eigenvectors of a real generalized\n*  eigenvalue problem A*x = lambda*B*x, by backward transformation on\n*  the computed eigenvectors of the balanced pair of matrices output by\n*  SGGBAL.\n*\n\n*  Arguments\n*  =========\n*\n*  JOB     (input) CHARACTER*1\n*          Specifies the type of backward transformation required:\n*          = 'N':  do nothing, return immediately;\n*          = 'P':  do backward transformation for permutation only;\n*          = 'S':  do backward transformation for scaling only;\n*          = 'B':  do backward transformations for both permutation and\n*                  scaling.\n*          JOB must be the same as the argument JOB supplied to SGGBAL.\n*\n*  SIDE    (input) CHARACTER*1\n*          = 'R':  V contains right eigenvectors;\n*          = 'L':  V contains left eigenvectors.\n*\n*  N       (input) INTEGER\n*          The number of rows of the matrix V.  N >= 0.\n*\n*  ILO     (input) INTEGER\n*  IHI     (input) INTEGER\n*          The integers ILO and IHI determined by SGGBAL.\n*          1 <= ILO <= IHI <= N, if N > 0; ILO=1 and IHI=0, if N=0.\n*\n*  LSCALE  (input) REAL array, dimension (N)\n*          Details of the permutations and/or scaling factors applied\n*          to the left side of A and B, as returned by SGGBAL.\n*\n*  RSCALE  (input) REAL array, dimension (N)\n*          Details of the permutations and/or scaling factors applied\n*          to the right side of A and B, as returned by SGGBAL.\n*\n*  M       (input) INTEGER\n*          The number of columns of the matrix V.  M >= 0.\n*\n*  V       (input/output) REAL array, dimension (LDV,M)\n*          On entry, the matrix of right or left eigenvectors to be\n*          transformed, as returned by STGEVC.\n*          On exit, V is overwritten by the transformed eigenvectors.\n*\n*  LDV     (input) INTEGER\n*          The leading dimension of the matrix V. LDV >= max(1,N).\n*\n*  INFO    (output) INTEGER\n*          = 0:  successful exit.\n*          < 0:  if INFO = -i, the i-th argument had an illegal value.\n*\n\n*  Further Details\n*  ===============\n*\n*  See R.C. Ward, Balancing the generalized eigenvalue problem,\n*                 SIAM J. Sci. Stat. Comp. 2 (1981), 141-152.\n*\n*  =====================================================================\n*\n*     .. Local Scalars ..\n      LOGICAL            LEFTV, RIGHTV\n      INTEGER            I, K\n*     ..\n*     .. External Functions ..\n      LOGICAL            LSAME\n      EXTERNAL           LSAME\n*     ..\n*     .. External Subroutines ..\n      EXTERNAL           SSCAL, SSWAP, XERBLA\n*     ..\n*     .. Intrinsic Functions ..\n      INTRINSIC          MAX\n*     ..\n\n");
      rb_exit(0);
    }
    if (rb_hash_aref(rb_options, sUsage) == Qtrue) {
      printf("%s\n", "USAGE:\n  info, v = NumRu::Lapack.sggbak( job, side, ilo, ihi, lscale, rscale, v, [:usage => usage, :help => help])\n");
      rb_exit(0);
    } 
  } else
    rb_options = Qnil;
  if (argc != 7)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 7)", argc);
  rblapack_job = argv[0];
  rblapack_side = argv[1];
  rblapack_ilo = argv[2];
  rblapack_ihi = argv[3];
  rblapack_lscale = argv[4];
  rblapack_rscale = argv[5];
  rblapack_v = argv[6];
  if (rb_options != Qnil) {
  }

  if (!NA_IsNArray(rblapack_v))
    rb_raise(rb_eArgError, "v (7th argument) must be NArray");
  if (NA_RANK(rblapack_v) != 2)
    rb_raise(rb_eArgError, "rank of v (7th argument) must be %d", 2);
  m = NA_SHAPE1(rblapack_v);
  ldv = NA_SHAPE0(rblapack_v);
  if (NA_TYPE(rblapack_v) != NA_SFLOAT)
    rblapack_v = na_change_type(rblapack_v, NA_SFLOAT);
  v = NA_PTR_TYPE(rblapack_v, real*);
  ilo = NUM2INT(rblapack_ilo);
  side = StringValueCStr(rblapack_side)[0];
  if (!NA_IsNArray(rblapack_rscale))
    rb_raise(rb_eArgError, "rscale (6th argument) must be NArray");
  if (NA_RANK(rblapack_rscale) != 1)
    rb_raise(rb_eArgError, "rank of rscale (6th argument) must be %d", 1);
  n = NA_SHAPE0(rblapack_rscale);
  if (NA_TYPE(rblapack_rscale) != NA_SFLOAT)
    rblapack_rscale = na_change_type(rblapack_rscale, NA_SFLOAT);
  rscale = NA_PTR_TYPE(rblapack_rscale, real*);
  job = StringValueCStr(rblapack_job)[0];
  if (!NA_IsNArray(rblapack_lscale))
    rb_raise(rb_eArgError, "lscale (5th argument) must be NArray");
  if (NA_RANK(rblapack_lscale) != 1)
    rb_raise(rb_eArgError, "rank of lscale (5th argument) must be %d", 1);
  if (NA_SHAPE0(rblapack_lscale) != n)
    rb_raise(rb_eRuntimeError, "shape 0 of lscale must be the same as shape 0 of rscale");
  if (NA_TYPE(rblapack_lscale) != NA_SFLOAT)
    rblapack_lscale = na_change_type(rblapack_lscale, NA_SFLOAT);
  lscale = NA_PTR_TYPE(rblapack_lscale, real*);
  ihi = NUM2INT(rblapack_ihi);
  {
    int shape[2];
    shape[0] = ldv;
    shape[1] = m;
    rblapack_v_out__ = na_make_object(NA_SFLOAT, 2, shape, cNArray);
  }
  v_out__ = NA_PTR_TYPE(rblapack_v_out__, real*);
  MEMCPY(v_out__, v, real, NA_TOTAL(rblapack_v));
  rblapack_v = rblapack_v_out__;
  v = v_out__;

  sggbak_(&job, &side, &n, &ilo, &ihi, lscale, rscale, &m, v, &ldv, &info);

  rblapack_info = INT2NUM(info);
  return rb_ary_new3(2, rblapack_info, rblapack_v);
}

void
init_lapack_sggbak(VALUE mLapack){
  rb_define_module_function(mLapack, "sggbak", rblapack_sggbak, -1);
  sHelp = ID2SYM(rb_intern("help"));
  sUsage = ID2SYM(rb_intern("usage"));
}
