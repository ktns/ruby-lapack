#include "rb_lapack.h"

extern VOID sgebak_(char *job, char *side, integer *n, integer *ilo, integer *ihi, real *scale, integer *m, real *v, integer *ldv, integer *info);

static VALUE sHelp, sUsage;

static VALUE
rblapack_sgebak(int argc, VALUE *argv, VALUE self){
  VALUE rblapack_job;
  char job; 
  VALUE rblapack_side;
  char side; 
  VALUE rblapack_ilo;
  integer ilo; 
  VALUE rblapack_ihi;
  integer ihi; 
  VALUE rblapack_scale;
  real *scale; 
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
      printf("%s\n", "USAGE:\n  info, v = NumRu::Lapack.sgebak( job, side, ilo, ihi, scale, v, [:usage => usage, :help => help])\n\n\nFORTRAN MANUAL\n      SUBROUTINE SGEBAK( JOB, SIDE, N, ILO, IHI, SCALE, M, V, LDV, INFO )\n\n*  Purpose\n*  =======\n*\n*  SGEBAK forms the right or left eigenvectors of a real general matrix\n*  by backward transformation on the computed eigenvectors of the\n*  balanced matrix output by SGEBAL.\n*\n\n*  Arguments\n*  =========\n*\n*  JOB     (input) CHARACTER*1\n*          Specifies the type of backward transformation required:\n*          = 'N', do nothing, return immediately;\n*          = 'P', do backward transformation for permutation only;\n*          = 'S', do backward transformation for scaling only;\n*          = 'B', do backward transformations for both permutation and\n*                 scaling.\n*          JOB must be the same as the argument JOB supplied to SGEBAL.\n*\n*  SIDE    (input) CHARACTER*1\n*          = 'R':  V contains right eigenvectors;\n*          = 'L':  V contains left eigenvectors.\n*\n*  N       (input) INTEGER\n*          The number of rows of the matrix V.  N >= 0.\n*\n*  ILO     (input) INTEGER\n*  IHI     (input) INTEGER\n*          The integers ILO and IHI determined by SGEBAL.\n*          1 <= ILO <= IHI <= N, if N > 0; ILO=1 and IHI=0, if N=0.\n*\n*  SCALE   (input) REAL array, dimension (N)\n*          Details of the permutation and scaling factors, as returned\n*          by SGEBAL.\n*\n*  M       (input) INTEGER\n*          The number of columns of the matrix V.  M >= 0.\n*\n*  V       (input/output) REAL array, dimension (LDV,M)\n*          On entry, the matrix of right or left eigenvectors to be\n*          transformed, as returned by SHSEIN or STREVC.\n*          On exit, V is overwritten by the transformed eigenvectors.\n*\n*  LDV     (input) INTEGER\n*          The leading dimension of the array V. LDV >= max(1,N).\n*\n*  INFO    (output) INTEGER\n*          = 0:  successful exit\n*          < 0:  if INFO = -i, the i-th argument had an illegal value.\n*\n\n*  =====================================================================\n*\n\n");
      rb_exit(0);
    }
    if (rb_hash_aref(rb_options, sUsage) == Qtrue) {
      printf("%s\n", "USAGE:\n  info, v = NumRu::Lapack.sgebak( job, side, ilo, ihi, scale, v, [:usage => usage, :help => help])\n");
      rb_exit(0);
    } 
  } else
    rb_options = Qnil;
  if (argc != 6)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 6)", argc);
  rblapack_job = argv[0];
  rblapack_side = argv[1];
  rblapack_ilo = argv[2];
  rblapack_ihi = argv[3];
  rblapack_scale = argv[4];
  rblapack_v = argv[5];
  if (rb_options != Qnil) {
  }

  if (!NA_IsNArray(rblapack_v))
    rb_raise(rb_eArgError, "v (6th argument) must be NArray");
  if (NA_RANK(rblapack_v) != 2)
    rb_raise(rb_eArgError, "rank of v (6th argument) must be %d", 2);
  m = NA_SHAPE1(rblapack_v);
  ldv = NA_SHAPE0(rblapack_v);
  if (NA_TYPE(rblapack_v) != NA_SFLOAT)
    rblapack_v = na_change_type(rblapack_v, NA_SFLOAT);
  v = NA_PTR_TYPE(rblapack_v, real*);
  if (!NA_IsNArray(rblapack_scale))
    rb_raise(rb_eArgError, "scale (5th argument) must be NArray");
  if (NA_RANK(rblapack_scale) != 1)
    rb_raise(rb_eArgError, "rank of scale (5th argument) must be %d", 1);
  n = NA_SHAPE0(rblapack_scale);
  if (NA_TYPE(rblapack_scale) != NA_SFLOAT)
    rblapack_scale = na_change_type(rblapack_scale, NA_SFLOAT);
  scale = NA_PTR_TYPE(rblapack_scale, real*);
  ilo = NUM2INT(rblapack_ilo);
  side = StringValueCStr(rblapack_side)[0];
  job = StringValueCStr(rblapack_job)[0];
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

  sgebak_(&job, &side, &n, &ilo, &ihi, scale, &m, v, &ldv, &info);

  rblapack_info = INT2NUM(info);
  return rb_ary_new3(2, rblapack_info, rblapack_v);
}

void
init_lapack_sgebak(VALUE mLapack){
  rb_define_module_function(mLapack, "sgebak", rblapack_sgebak, -1);
  sHelp = ID2SYM(rb_intern("help"));
  sUsage = ID2SYM(rb_intern("usage"));
}
