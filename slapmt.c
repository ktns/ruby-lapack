#include "rb_lapack.h"

extern VOID slapmt_(logical *forwrd, integer *m, integer *n, real *x, integer *ldx, integer *k);

static VALUE sHelp, sUsage;

static VALUE
rblapack_slapmt(int argc, VALUE *argv, VALUE self){
  VALUE rblapack_forwrd;
  logical forwrd; 
  VALUE rblapack_m;
  integer m; 
  VALUE rblapack_x;
  real *x; 
  VALUE rblapack_k;
  integer *k; 
  VALUE rblapack_x_out__;
  real *x_out__;
  VALUE rblapack_k_out__;
  integer *k_out__;

  integer ldx;
  integer n;

  VALUE rb_options;
  if (argc > 0 && TYPE(argv[argc-1]) == T_HASH) {
    argc--;
    rb_options = argv[argc];
    if (rb_hash_aref(rb_options, sHelp) == Qtrue) {
      printf("%s\n", "USAGE:\n  x, k = NumRu::Lapack.slapmt( forwrd, m, x, k, [:usage => usage, :help => help])\n\n\nFORTRAN MANUAL\n      SUBROUTINE SLAPMT( FORWRD, M, N, X, LDX, K )\n\n*  Purpose\n*  =======\n*\n*  SLAPMT rearranges the columns of the M by N matrix X as specified\n*  by the permutation K(1),K(2),...,K(N) of the integers 1,...,N.\n*  If FORWRD = .TRUE.,  forward permutation:\n*\n*       X(*,K(J)) is moved X(*,J) for J = 1,2,...,N.\n*\n*  If FORWRD = .FALSE., backward permutation:\n*\n*       X(*,J) is moved to X(*,K(J)) for J = 1,2,...,N.\n*\n\n*  Arguments\n*  =========\n*\n*  FORWRD  (input) LOGICAL\n*          = .TRUE., forward permutation\n*          = .FALSE., backward permutation\n*\n*  M       (input) INTEGER\n*          The number of rows of the matrix X. M >= 0.\n*\n*  N       (input) INTEGER\n*          The number of columns of the matrix X. N >= 0.\n*\n*  X       (input/output) REAL array, dimension (LDX,N)\n*          On entry, the M by N matrix X.\n*          On exit, X contains the permuted matrix X.\n*\n*  LDX     (input) INTEGER\n*          The leading dimension of the array X, LDX >= MAX(1,M).\n*\n*  K       (input/output) INTEGER array, dimension (N)\n*          On entry, K contains the permutation vector. K is used as\n*          internal workspace, but reset to its original value on\n*          output.\n*\n\n*  =====================================================================\n*\n*     .. Local Scalars ..\n      INTEGER            I, II, J, IN\n      REAL               TEMP\n*     ..\n\n");
      rb_exit(0);
    }
    if (rb_hash_aref(rb_options, sUsage) == Qtrue) {
      printf("%s\n", "USAGE:\n  x, k = NumRu::Lapack.slapmt( forwrd, m, x, k, [:usage => usage, :help => help])\n");
      rb_exit(0);
    } 
  } else
    rb_options = Qnil;
  if (argc != 4)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 4)", argc);
  rblapack_forwrd = argv[0];
  rblapack_m = argv[1];
  rblapack_x = argv[2];
  rblapack_k = argv[3];
  if (rb_options != Qnil) {
  }

  if (!NA_IsNArray(rblapack_k))
    rb_raise(rb_eArgError, "k (4th argument) must be NArray");
  if (NA_RANK(rblapack_k) != 1)
    rb_raise(rb_eArgError, "rank of k (4th argument) must be %d", 1);
  n = NA_SHAPE0(rblapack_k);
  if (NA_TYPE(rblapack_k) != NA_LINT)
    rblapack_k = na_change_type(rblapack_k, NA_LINT);
  k = NA_PTR_TYPE(rblapack_k, integer*);
  forwrd = (rblapack_forwrd == Qtrue);
  if (!NA_IsNArray(rblapack_x))
    rb_raise(rb_eArgError, "x (3th argument) must be NArray");
  if (NA_RANK(rblapack_x) != 2)
    rb_raise(rb_eArgError, "rank of x (3th argument) must be %d", 2);
  if (NA_SHAPE1(rblapack_x) != n)
    rb_raise(rb_eRuntimeError, "shape 1 of x must be the same as shape 0 of k");
  ldx = NA_SHAPE0(rblapack_x);
  if (NA_TYPE(rblapack_x) != NA_SFLOAT)
    rblapack_x = na_change_type(rblapack_x, NA_SFLOAT);
  x = NA_PTR_TYPE(rblapack_x, real*);
  m = NUM2INT(rblapack_m);
  {
    int shape[2];
    shape[0] = ldx;
    shape[1] = n;
    rblapack_x_out__ = na_make_object(NA_SFLOAT, 2, shape, cNArray);
  }
  x_out__ = NA_PTR_TYPE(rblapack_x_out__, real*);
  MEMCPY(x_out__, x, real, NA_TOTAL(rblapack_x));
  rblapack_x = rblapack_x_out__;
  x = x_out__;
  {
    int shape[1];
    shape[0] = n;
    rblapack_k_out__ = na_make_object(NA_LINT, 1, shape, cNArray);
  }
  k_out__ = NA_PTR_TYPE(rblapack_k_out__, integer*);
  MEMCPY(k_out__, k, integer, NA_TOTAL(rblapack_k));
  rblapack_k = rblapack_k_out__;
  k = k_out__;

  slapmt_(&forwrd, &m, &n, x, &ldx, k);

  return rb_ary_new3(2, rblapack_x, rblapack_k);
}

void
init_lapack_slapmt(VALUE mLapack){
  rb_define_module_function(mLapack, "slapmt", rblapack_slapmt, -1);
  sHelp = ID2SYM(rb_intern("help"));
  sUsage = ID2SYM(rb_intern("usage"));
}
