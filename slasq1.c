#include "rb_lapack.h"

static VALUE
rb_slasq1(int argc, VALUE *argv, VALUE self){
  VALUE rb_d;
  real *d; 
  VALUE rb_e;
  real *e; 
  VALUE rb_info;
  integer info; 
  VALUE rb_d_out__;
  real *d_out__;
  VALUE rb_e_out__;
  real *e_out__;
  real *work;

  integer n;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  info, d, e = NumRu::Lapack.slasq1( d, e)\n    or\n  NumRu::Lapack.slasq1  # print help\n\n\nFORTRAN MANUAL\n      SUBROUTINE SLASQ1( N, D, E, WORK, INFO )\n\n*  Purpose\n*  =======\n*\n*  SLASQ1 computes the singular values of a real N-by-N bidiagonal\n*  matrix with diagonal D and off-diagonal E. The singular values\n*  are computed to high relative accuracy, in the absence of\n*  denormalization, underflow and overflow. The algorithm was first\n*  presented in\n*\n*  \"Accurate singular values and differential qd algorithms\" by K. V.\n*  Fernando and B. N. Parlett, Numer. Math., Vol-67, No. 2, pp. 191-230,\n*  1994,\n*\n*  and the present implementation is described in \"An implementation of\n*  the dqds Algorithm (Positive Case)\", LAPACK Working Note.\n*\n\n*  Arguments\n*  =========\n*\n*  N     (input) INTEGER\n*        The number of rows and columns in the matrix. N >= 0.\n*\n*  D     (input/output) REAL array, dimension (N)\n*        On entry, D contains the diagonal elements of the\n*        bidiagonal matrix whose SVD is desired. On normal exit,\n*        D contains the singular values in decreasing order.\n*\n*  E     (input/output) REAL array, dimension (N)\n*        On entry, elements E(1:N-1) contain the off-diagonal elements\n*        of the bidiagonal matrix whose SVD is desired.\n*        On exit, E is overwritten.\n*\n*  WORK  (workspace) REAL array, dimension (4*N)\n*\n*  INFO  (output) INTEGER\n*        = 0: successful exit\n*        < 0: if INFO = -i, the i-th argument had an illegal value\n*        > 0: the algorithm failed\n*             = 1, a split was marked by a positive value in E\n*             = 2, current block of Z not diagonalized after 30*N\n*                  iterations (in inner while loop)\n*             = 3, termination criterion of outer while loop not met \n*                  (program created more than N unreduced blocks)\n*\n\n*  =====================================================================\n*\n\n");
    return Qnil;
  }
  if (argc != 2)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 2)", argc);
  rb_d = argv[0];
  rb_e = argv[1];

  if (!NA_IsNArray(rb_d))
    rb_raise(rb_eArgError, "d (1th argument) must be NArray");
  if (NA_RANK(rb_d) != 1)
    rb_raise(rb_eArgError, "rank of d (1th argument) must be %d", 1);
  n = NA_SHAPE0(rb_d);
  if (NA_TYPE(rb_d) != NA_SFLOAT)
    rb_d = na_change_type(rb_d, NA_SFLOAT);
  d = NA_PTR_TYPE(rb_d, real*);
  if (!NA_IsNArray(rb_e))
    rb_raise(rb_eArgError, "e (2th argument) must be NArray");
  if (NA_RANK(rb_e) != 1)
    rb_raise(rb_eArgError, "rank of e (2th argument) must be %d", 1);
  if (NA_SHAPE0(rb_e) != n)
    rb_raise(rb_eRuntimeError, "shape 0 of e must be the same as shape 0 of d");
  if (NA_TYPE(rb_e) != NA_SFLOAT)
    rb_e = na_change_type(rb_e, NA_SFLOAT);
  e = NA_PTR_TYPE(rb_e, real*);
  {
    int shape[1];
    shape[0] = n;
    rb_d_out__ = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  d_out__ = NA_PTR_TYPE(rb_d_out__, real*);
  MEMCPY(d_out__, d, real, NA_TOTAL(rb_d));
  rb_d = rb_d_out__;
  d = d_out__;
  {
    int shape[1];
    shape[0] = n;
    rb_e_out__ = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  e_out__ = NA_PTR_TYPE(rb_e_out__, real*);
  MEMCPY(e_out__, e, real, NA_TOTAL(rb_e));
  rb_e = rb_e_out__;
  e = e_out__;
  work = ALLOC_N(real, (4*n));

  slasq1_(&n, d, e, work, &info);

  free(work);
  rb_info = INT2NUM(info);
  return rb_ary_new3(3, rb_info, rb_d, rb_e);
}

void
init_lapack_slasq1(VALUE mLapack){
  rb_define_module_function(mLapack, "slasq1", rb_slasq1, -1);
}
