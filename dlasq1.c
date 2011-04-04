#include "rb_lapack.h"

extern VOID dlasq1_(integer *n, doublereal *d, doublereal *e, doublereal *work, integer *info);

static VALUE
rb_dlasq1(int argc, VALUE *argv, VALUE self){
  VALUE rb_d;
  doublereal *d; 
  VALUE rb_e;
  doublereal *e; 
  VALUE rb_info;
  integer info; 
  VALUE rb_d_out__;
  doublereal *d_out__;
  VALUE rb_e_out__;
  doublereal *e_out__;
  doublereal *work;

  integer n;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  info, d, e = NumRu::Lapack.dlasq1( d, e)\n    or\n  NumRu::Lapack.dlasq1  # print help\n\n\nFORTRAN MANUAL\n      SUBROUTINE DLASQ1( N, D, E, WORK, INFO )\n\n*  Purpose\n*  =======\n*\n*  DLASQ1 computes the singular values of a real N-by-N bidiagonal\n*  matrix with diagonal D and off-diagonal E. The singular values\n*  are computed to high relative accuracy, in the absence of\n*  denormalization, underflow and overflow. The algorithm was first\n*  presented in\n*\n*  \"Accurate singular values and differential qd algorithms\" by K. V.\n*  Fernando and B. N. Parlett, Numer. Math., Vol-67, No. 2, pp. 191-230,\n*  1994,\n*\n*  and the present implementation is described in \"An implementation of\n*  the dqds Algorithm (Positive Case)\", LAPACK Working Note.\n*\n\n*  Arguments\n*  =========\n*\n*  N     (input) INTEGER\n*        The number of rows and columns in the matrix. N >= 0.\n*\n*  D     (input/output) DOUBLE PRECISION array, dimension (N)\n*        On entry, D contains the diagonal elements of the\n*        bidiagonal matrix whose SVD is desired. On normal exit,\n*        D contains the singular values in decreasing order.\n*\n*  E     (input/output) DOUBLE PRECISION array, dimension (N)\n*        On entry, elements E(1:N-1) contain the off-diagonal elements\n*        of the bidiagonal matrix whose SVD is desired.\n*        On exit, E is overwritten.\n*\n*  WORK  (workspace) DOUBLE PRECISION array, dimension (4*N)\n*\n*  INFO  (output) INTEGER\n*        = 0: successful exit\n*        < 0: if INFO = -i, the i-th argument had an illegal value\n*        > 0: the algorithm failed\n*             = 1, a split was marked by a positive value in E\n*             = 2, current block of Z not diagonalized after 30*N\n*                  iterations (in inner while loop)\n*             = 3, termination criterion of outer while loop not met \n*                  (program created more than N unreduced blocks)\n*\n\n*  =====================================================================\n*\n\n");
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
  if (NA_TYPE(rb_d) != NA_DFLOAT)
    rb_d = na_change_type(rb_d, NA_DFLOAT);
  d = NA_PTR_TYPE(rb_d, doublereal*);
  if (!NA_IsNArray(rb_e))
    rb_raise(rb_eArgError, "e (2th argument) must be NArray");
  if (NA_RANK(rb_e) != 1)
    rb_raise(rb_eArgError, "rank of e (2th argument) must be %d", 1);
  if (NA_SHAPE0(rb_e) != n)
    rb_raise(rb_eRuntimeError, "shape 0 of e must be the same as shape 0 of d");
  if (NA_TYPE(rb_e) != NA_DFLOAT)
    rb_e = na_change_type(rb_e, NA_DFLOAT);
  e = NA_PTR_TYPE(rb_e, doublereal*);
  {
    int shape[1];
    shape[0] = n;
    rb_d_out__ = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  d_out__ = NA_PTR_TYPE(rb_d_out__, doublereal*);
  MEMCPY(d_out__, d, doublereal, NA_TOTAL(rb_d));
  rb_d = rb_d_out__;
  d = d_out__;
  {
    int shape[1];
    shape[0] = n;
    rb_e_out__ = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  e_out__ = NA_PTR_TYPE(rb_e_out__, doublereal*);
  MEMCPY(e_out__, e, doublereal, NA_TOTAL(rb_e));
  rb_e = rb_e_out__;
  e = e_out__;
  work = ALLOC_N(doublereal, (4*n));

  dlasq1_(&n, d, e, work, &info);

  free(work);
  rb_info = INT2NUM(info);
  return rb_ary_new3(3, rb_info, rb_d, rb_e);
}

void
init_lapack_dlasq1(VALUE mLapack){
  rb_define_module_function(mLapack, "dlasq1", rb_dlasq1, -1);
}
