#include "rb_lapack.h"

extern VOID slasq1_(integer *n, real *d, real *e, real *work, integer *info);

static VALUE sHelp, sUsage;

static VALUE
rblapack_slasq1(int argc, VALUE *argv, VALUE self){
  VALUE rblapack_d;
  real *d; 
  VALUE rblapack_e;
  real *e; 
  VALUE rblapack_info;
  integer info; 
  VALUE rblapack_d_out__;
  real *d_out__;
  VALUE rblapack_e_out__;
  real *e_out__;
  real *work;

  integer n;

  VALUE rb_options;
  if (argc > 0 && TYPE(argv[argc-1]) == T_HASH) {
    argc--;
    rb_options = argv[argc];
    if (rb_hash_aref(rb_options, sHelp) == Qtrue) {
      printf("%s\n", "USAGE:\n  info, d, e = NumRu::Lapack.slasq1( d, e, [:usage => usage, :help => help])\n\n\nFORTRAN MANUAL\n      SUBROUTINE SLASQ1( N, D, E, WORK, INFO )\n\n*  Purpose\n*  =======\n*\n*  SLASQ1 computes the singular values of a real N-by-N bidiagonal\n*  matrix with diagonal D and off-diagonal E. The singular values\n*  are computed to high relative accuracy, in the absence of\n*  denormalization, underflow and overflow. The algorithm was first\n*  presented in\n*\n*  \"Accurate singular values and differential qd algorithms\" by K. V.\n*  Fernando and B. N. Parlett, Numer. Math., Vol-67, No. 2, pp. 191-230,\n*  1994,\n*\n*  and the present implementation is described in \"An implementation of\n*  the dqds Algorithm (Positive Case)\", LAPACK Working Note.\n*\n\n*  Arguments\n*  =========\n*\n*  N     (input) INTEGER\n*        The number of rows and columns in the matrix. N >= 0.\n*\n*  D     (input/output) REAL array, dimension (N)\n*        On entry, D contains the diagonal elements of the\n*        bidiagonal matrix whose SVD is desired. On normal exit,\n*        D contains the singular values in decreasing order.\n*\n*  E     (input/output) REAL array, dimension (N)\n*        On entry, elements E(1:N-1) contain the off-diagonal elements\n*        of the bidiagonal matrix whose SVD is desired.\n*        On exit, E is overwritten.\n*\n*  WORK  (workspace) REAL array, dimension (4*N)\n*\n*  INFO  (output) INTEGER\n*        = 0: successful exit\n*        < 0: if INFO = -i, the i-th argument had an illegal value\n*        > 0: the algorithm failed\n*             = 1, a split was marked by a positive value in E\n*             = 2, current block of Z not diagonalized after 30*N\n*                  iterations (in inner while loop)\n*             = 3, termination criterion of outer while loop not met \n*                  (program created more than N unreduced blocks)\n*\n\n*  =====================================================================\n*\n\n");
      rb_exit(0);
    }
    if (rb_hash_aref(rb_options, sUsage) == Qtrue) {
      printf("%s\n", "USAGE:\n  info, d, e = NumRu::Lapack.slasq1( d, e, [:usage => usage, :help => help])\n");
      rb_exit(0);
    } 
  } else
    rb_options = Qnil;
  if (argc != 2)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 2)", argc);
  rblapack_d = argv[0];
  rblapack_e = argv[1];
  if (rb_options != Qnil) {
  }

  if (!NA_IsNArray(rblapack_d))
    rb_raise(rb_eArgError, "d (1th argument) must be NArray");
  if (NA_RANK(rblapack_d) != 1)
    rb_raise(rb_eArgError, "rank of d (1th argument) must be %d", 1);
  n = NA_SHAPE0(rblapack_d);
  if (NA_TYPE(rblapack_d) != NA_SFLOAT)
    rblapack_d = na_change_type(rblapack_d, NA_SFLOAT);
  d = NA_PTR_TYPE(rblapack_d, real*);
  if (!NA_IsNArray(rblapack_e))
    rb_raise(rb_eArgError, "e (2th argument) must be NArray");
  if (NA_RANK(rblapack_e) != 1)
    rb_raise(rb_eArgError, "rank of e (2th argument) must be %d", 1);
  if (NA_SHAPE0(rblapack_e) != n)
    rb_raise(rb_eRuntimeError, "shape 0 of e must be the same as shape 0 of d");
  if (NA_TYPE(rblapack_e) != NA_SFLOAT)
    rblapack_e = na_change_type(rblapack_e, NA_SFLOAT);
  e = NA_PTR_TYPE(rblapack_e, real*);
  {
    int shape[1];
    shape[0] = n;
    rblapack_d_out__ = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  d_out__ = NA_PTR_TYPE(rblapack_d_out__, real*);
  MEMCPY(d_out__, d, real, NA_TOTAL(rblapack_d));
  rblapack_d = rblapack_d_out__;
  d = d_out__;
  {
    int shape[1];
    shape[0] = n;
    rblapack_e_out__ = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  e_out__ = NA_PTR_TYPE(rblapack_e_out__, real*);
  MEMCPY(e_out__, e, real, NA_TOTAL(rblapack_e));
  rblapack_e = rblapack_e_out__;
  e = e_out__;
  work = ALLOC_N(real, (4*n));

  slasq1_(&n, d, e, work, &info);

  free(work);
  rblapack_info = INT2NUM(info);
  return rb_ary_new3(3, rblapack_info, rblapack_d, rblapack_e);
}

void
init_lapack_slasq1(VALUE mLapack){
  rb_define_module_function(mLapack, "slasq1", rblapack_slasq1, -1);
  sHelp = ID2SYM(rb_intern("help"));
  sUsage = ID2SYM(rb_intern("usage"));
}
