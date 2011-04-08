#include "rb_lapack.h"

extern VOID slagts_(integer *job, integer *n, real *a, real *b, real *c, real *d, integer *in, real *y, real *tol, integer *info);

static VALUE sHelp, sUsage;

static VALUE
rblapack_slagts(int argc, VALUE *argv, VALUE self){
  VALUE rblapack_job;
  integer job; 
  VALUE rblapack_a;
  real *a; 
  VALUE rblapack_b;
  real *b; 
  VALUE rblapack_c;
  real *c; 
  VALUE rblapack_d;
  real *d; 
  VALUE rblapack_in;
  integer *in; 
  VALUE rblapack_y;
  real *y; 
  VALUE rblapack_tol;
  real tol; 
  VALUE rblapack_info;
  integer info; 
  VALUE rblapack_y_out__;
  real *y_out__;

  integer n;

  VALUE rb_options;
  if (argc > 0 && TYPE(argv[argc-1]) == T_HASH) {
    argc--;
    rb_options = argv[argc];
    if (rb_hash_aref(rb_options, sHelp) == Qtrue) {
      printf("%s\n", "USAGE:\n  info, y, tol = NumRu::Lapack.slagts( job, a, b, c, d, in, y, tol, [:usage => usage, :help => help])\n\n\nFORTRAN MANUAL\n      SUBROUTINE SLAGTS( JOB, N, A, B, C, D, IN, Y, TOL, INFO )\n\n*  Purpose\n*  =======\n*\n*  SLAGTS may be used to solve one of the systems of equations\n*\n*     (T - lambda*I)*x = y   or   (T - lambda*I)'*x = y,\n*\n*  where T is an n by n tridiagonal matrix, for x, following the\n*  factorization of (T - lambda*I) as\n*\n*     (T - lambda*I) = P*L*U ,\n*\n*  by routine SLAGTF. The choice of equation to be solved is\n*  controlled by the argument JOB, and in each case there is an option\n*  to perturb zero or very small diagonal elements of U, this option\n*  being intended for use in applications such as inverse iteration.\n*\n\n*  Arguments\n*  =========\n*\n*  JOB     (input) INTEGER\n*          Specifies the job to be performed by SLAGTS as follows:\n*          =  1: The equations  (T - lambda*I)x = y  are to be solved,\n*                but diagonal elements of U are not to be perturbed.\n*          = -1: The equations  (T - lambda*I)x = y  are to be solved\n*                and, if overflow would otherwise occur, the diagonal\n*                elements of U are to be perturbed. See argument TOL\n*                below.\n*          =  2: The equations  (T - lambda*I)'x = y  are to be solved,\n*                but diagonal elements of U are not to be perturbed.\n*          = -2: The equations  (T - lambda*I)'x = y  are to be solved\n*                and, if overflow would otherwise occur, the diagonal\n*                elements of U are to be perturbed. See argument TOL\n*                below.\n*\n*  N       (input) INTEGER\n*          The order of the matrix T.\n*\n*  A       (input) REAL array, dimension (N)\n*          On entry, A must contain the diagonal elements of U as\n*          returned from SLAGTF.\n*\n*  B       (input) REAL array, dimension (N-1)\n*          On entry, B must contain the first super-diagonal elements of\n*          U as returned from SLAGTF.\n*\n*  C       (input) REAL array, dimension (N-1)\n*          On entry, C must contain the sub-diagonal elements of L as\n*          returned from SLAGTF.\n*\n*  D       (input) REAL array, dimension (N-2)\n*          On entry, D must contain the second super-diagonal elements\n*          of U as returned from SLAGTF.\n*\n*  IN      (input) INTEGER array, dimension (N)\n*          On entry, IN must contain details of the matrix P as returned\n*          from SLAGTF.\n*\n*  Y       (input/output) REAL array, dimension (N)\n*          On entry, the right hand side vector y.\n*          On exit, Y is overwritten by the solution vector x.\n*\n*  TOL     (input/output) REAL\n*          On entry, with  JOB .lt. 0, TOL should be the minimum\n*          perturbation to be made to very small diagonal elements of U.\n*          TOL should normally be chosen as about eps*norm(U), where eps\n*          is the relative machine precision, but if TOL is supplied as\n*          non-positive, then it is reset to eps*max( abs( u(i,j) ) ).\n*          If  JOB .gt. 0  then TOL is not referenced.\n*\n*          On exit, TOL is changed as described above, only if TOL is\n*          non-positive on entry. Otherwise TOL is unchanged.\n*\n*  INFO    (output) INTEGER\n*          = 0   : successful exit\n*          .lt. 0: if INFO = -i, the i-th argument had an illegal value\n*          .gt. 0: overflow would occur when computing the INFO(th)\n*                  element of the solution vector x. This can only occur\n*                  when JOB is supplied as positive and either means\n*                  that a diagonal element of U is very small, or that\n*                  the elements of the right-hand side vector y are very\n*                  large.\n*\n\n*  =====================================================================\n*\n\n");
      rb_exit(0);
    }
    if (rb_hash_aref(rb_options, sUsage) == Qtrue) {
      printf("%s\n", "USAGE:\n  info, y, tol = NumRu::Lapack.slagts( job, a, b, c, d, in, y, tol, [:usage => usage, :help => help])\n");
      rb_exit(0);
    } 
  } else
    rb_options = Qnil;
  if (argc != 8)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 8)", argc);
  rblapack_job = argv[0];
  rblapack_a = argv[1];
  rblapack_b = argv[2];
  rblapack_c = argv[3];
  rblapack_d = argv[4];
  rblapack_in = argv[5];
  rblapack_y = argv[6];
  rblapack_tol = argv[7];
  if (rb_options != Qnil) {
  }

  if (!NA_IsNArray(rblapack_a))
    rb_raise(rb_eArgError, "a (2th argument) must be NArray");
  if (NA_RANK(rblapack_a) != 1)
    rb_raise(rb_eArgError, "rank of a (2th argument) must be %d", 1);
  n = NA_SHAPE0(rblapack_a);
  if (NA_TYPE(rblapack_a) != NA_SFLOAT)
    rblapack_a = na_change_type(rblapack_a, NA_SFLOAT);
  a = NA_PTR_TYPE(rblapack_a, real*);
  tol = (real)NUM2DBL(rblapack_tol);
  if (!NA_IsNArray(rblapack_y))
    rb_raise(rb_eArgError, "y (7th argument) must be NArray");
  if (NA_RANK(rblapack_y) != 1)
    rb_raise(rb_eArgError, "rank of y (7th argument) must be %d", 1);
  if (NA_SHAPE0(rblapack_y) != n)
    rb_raise(rb_eRuntimeError, "shape 0 of y must be the same as shape 0 of a");
  if (NA_TYPE(rblapack_y) != NA_SFLOAT)
    rblapack_y = na_change_type(rblapack_y, NA_SFLOAT);
  y = NA_PTR_TYPE(rblapack_y, real*);
  if (!NA_IsNArray(rblapack_in))
    rb_raise(rb_eArgError, "in (6th argument) must be NArray");
  if (NA_RANK(rblapack_in) != 1)
    rb_raise(rb_eArgError, "rank of in (6th argument) must be %d", 1);
  if (NA_SHAPE0(rblapack_in) != n)
    rb_raise(rb_eRuntimeError, "shape 0 of in must be the same as shape 0 of a");
  if (NA_TYPE(rblapack_in) != NA_LINT)
    rblapack_in = na_change_type(rblapack_in, NA_LINT);
  in = NA_PTR_TYPE(rblapack_in, integer*);
  job = NUM2INT(rblapack_job);
  if (!NA_IsNArray(rblapack_c))
    rb_raise(rb_eArgError, "c (4th argument) must be NArray");
  if (NA_RANK(rblapack_c) != 1)
    rb_raise(rb_eArgError, "rank of c (4th argument) must be %d", 1);
  if (NA_SHAPE0(rblapack_c) != (n-1))
    rb_raise(rb_eRuntimeError, "shape 0 of c must be %d", n-1);
  if (NA_TYPE(rblapack_c) != NA_SFLOAT)
    rblapack_c = na_change_type(rblapack_c, NA_SFLOAT);
  c = NA_PTR_TYPE(rblapack_c, real*);
  if (!NA_IsNArray(rblapack_b))
    rb_raise(rb_eArgError, "b (3th argument) must be NArray");
  if (NA_RANK(rblapack_b) != 1)
    rb_raise(rb_eArgError, "rank of b (3th argument) must be %d", 1);
  if (NA_SHAPE0(rblapack_b) != (n-1))
    rb_raise(rb_eRuntimeError, "shape 0 of b must be %d", n-1);
  if (NA_TYPE(rblapack_b) != NA_SFLOAT)
    rblapack_b = na_change_type(rblapack_b, NA_SFLOAT);
  b = NA_PTR_TYPE(rblapack_b, real*);
  if (!NA_IsNArray(rblapack_d))
    rb_raise(rb_eArgError, "d (5th argument) must be NArray");
  if (NA_RANK(rblapack_d) != 1)
    rb_raise(rb_eArgError, "rank of d (5th argument) must be %d", 1);
  if (NA_SHAPE0(rblapack_d) != (n-2))
    rb_raise(rb_eRuntimeError, "shape 0 of d must be %d", n-2);
  if (NA_TYPE(rblapack_d) != NA_SFLOAT)
    rblapack_d = na_change_type(rblapack_d, NA_SFLOAT);
  d = NA_PTR_TYPE(rblapack_d, real*);
  {
    int shape[1];
    shape[0] = n;
    rblapack_y_out__ = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  y_out__ = NA_PTR_TYPE(rblapack_y_out__, real*);
  MEMCPY(y_out__, y, real, NA_TOTAL(rblapack_y));
  rblapack_y = rblapack_y_out__;
  y = y_out__;

  slagts_(&job, &n, a, b, c, d, in, y, &tol, &info);

  rblapack_info = INT2NUM(info);
  rblapack_tol = rb_float_new((double)tol);
  return rb_ary_new3(3, rblapack_info, rblapack_y, rblapack_tol);
}

void
init_lapack_slagts(VALUE mLapack){
  rb_define_module_function(mLapack, "slagts", rblapack_slagts, -1);
  sHelp = ID2SYM(rb_intern("help"));
  sUsage = ID2SYM(rb_intern("usage"));
}
