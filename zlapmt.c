#include "rb_lapack.h"

extern VOID zlapmt_(logical *forwrd, integer *m, integer *n, doublecomplex *x, integer *ldx, integer *k);

static VALUE
rb_zlapmt(int argc, VALUE *argv, VALUE self){
  VALUE rb_forwrd;
  logical forwrd; 
  VALUE rb_m;
  integer m; 
  VALUE rb_x;
  doublecomplex *x; 
  VALUE rb_k;
  integer *k; 
  VALUE rb_x_out__;
  doublecomplex *x_out__;
  VALUE rb_k_out__;
  integer *k_out__;

  integer ldx;
  integer n;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  x, k = NumRu::Lapack.zlapmt( forwrd, m, x, k)\n    or\n  NumRu::Lapack.zlapmt  # print help\n\n\nFORTRAN MANUAL\n      SUBROUTINE ZLAPMT( FORWRD, M, N, X, LDX, K )\n\n*  Purpose\n*  =======\n*\n*  ZLAPMT rearranges the columns of the M by N matrix X as specified\n*  by the permutation K(1),K(2),...,K(N) of the integers 1,...,N.\n*  If FORWRD = .TRUE.,  forward permutation:\n*\n*       X(*,K(J)) is moved X(*,J) for J = 1,2,...,N.\n*\n*  If FORWRD = .FALSE., backward permutation:\n*\n*       X(*,J) is moved to X(*,K(J)) for J = 1,2,...,N.\n*\n\n*  Arguments\n*  =========\n*\n*  FORWRD  (input) LOGICAL\n*          = .TRUE., forward permutation\n*          = .FALSE., backward permutation\n*\n*  M       (input) INTEGER\n*          The number of rows of the matrix X. M >= 0.\n*\n*  N       (input) INTEGER\n*          The number of columns of the matrix X. N >= 0.\n*\n*  X       (input/output) COMPLEX*16 array, dimension (LDX,N)\n*          On entry, the M by N matrix X.\n*          On exit, X contains the permuted matrix X.\n*\n*  LDX     (input) INTEGER\n*          The leading dimension of the array X, LDX >= MAX(1,M).\n*\n*  K       (input/output) INTEGER array, dimension (N)\n*          On entry, K contains the permutation vector. K is used as\n*          internal workspace, but reset to its original value on\n*          output.\n*\n\n*  =====================================================================\n*\n*     .. Local Scalars ..\n      INTEGER            I, II, IN, J\n      COMPLEX*16         TEMP\n*     ..\n\n");
    return Qnil;
  }
  if (argc != 4)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 4)", argc);
  rb_forwrd = argv[0];
  rb_m = argv[1];
  rb_x = argv[2];
  rb_k = argv[3];

  if (!NA_IsNArray(rb_k))
    rb_raise(rb_eArgError, "k (4th argument) must be NArray");
  if (NA_RANK(rb_k) != 1)
    rb_raise(rb_eArgError, "rank of k (4th argument) must be %d", 1);
  n = NA_SHAPE0(rb_k);
  if (NA_TYPE(rb_k) != NA_LINT)
    rb_k = na_change_type(rb_k, NA_LINT);
  k = NA_PTR_TYPE(rb_k, integer*);
  forwrd = (rb_forwrd == Qtrue);
  if (!NA_IsNArray(rb_x))
    rb_raise(rb_eArgError, "x (3th argument) must be NArray");
  if (NA_RANK(rb_x) != 2)
    rb_raise(rb_eArgError, "rank of x (3th argument) must be %d", 2);
  if (NA_SHAPE1(rb_x) != n)
    rb_raise(rb_eRuntimeError, "shape 1 of x must be the same as shape 0 of k");
  ldx = NA_SHAPE0(rb_x);
  if (NA_TYPE(rb_x) != NA_DCOMPLEX)
    rb_x = na_change_type(rb_x, NA_DCOMPLEX);
  x = NA_PTR_TYPE(rb_x, doublecomplex*);
  m = NUM2INT(rb_m);
  {
    int shape[2];
    shape[0] = ldx;
    shape[1] = n;
    rb_x_out__ = na_make_object(NA_DCOMPLEX, 2, shape, cNArray);
  }
  x_out__ = NA_PTR_TYPE(rb_x_out__, doublecomplex*);
  MEMCPY(x_out__, x, doublecomplex, NA_TOTAL(rb_x));
  rb_x = rb_x_out__;
  x = x_out__;
  {
    int shape[1];
    shape[0] = n;
    rb_k_out__ = na_make_object(NA_LINT, 1, shape, cNArray);
  }
  k_out__ = NA_PTR_TYPE(rb_k_out__, integer*);
  MEMCPY(k_out__, k, integer, NA_TOTAL(rb_k));
  rb_k = rb_k_out__;
  k = k_out__;

  zlapmt_(&forwrd, &m, &n, x, &ldx, k);

  return rb_ary_new3(2, rb_x, rb_k);
}

void
init_lapack_zlapmt(VALUE mLapack){
  rb_define_module_function(mLapack, "zlapmt", rb_zlapmt, -1);
}
