#include "rb_lapack.h"

extern VOID dlarrr_(integer *n, doublereal *d, doublereal *e, integer *info);

static VALUE
rb_dlarrr(int argc, VALUE *argv, VALUE self){
  VALUE rb_d;
  doublereal *d; 
  VALUE rb_e;
  doublereal *e; 
  VALUE rb_info;
  integer info; 
  VALUE rb_e_out__;
  doublereal *e_out__;

  integer n;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  info, e = NumRu::Lapack.dlarrr( d, e)\n    or\n  NumRu::Lapack.dlarrr  # print help\n\n\nFORTRAN MANUAL\n      SUBROUTINE DLARRR( N, D, E, INFO )\n\n*  Purpose\n*  =======\n*\n*  Perform tests to decide whether the symmetric tridiagonal matrix T\n*  warrants expensive computations which guarantee high relative accuracy\n*  in the eigenvalues.\n*\n\n*  Arguments\n*  =========\n*\n*  N       (input) INTEGER\n*          The order of the matrix. N > 0.\n*\n*  D       (input) DOUBLE PRECISION array, dimension (N)\n*          The N diagonal elements of the tridiagonal matrix T.\n*\n*  E       (input/output) DOUBLE PRECISION array, dimension (N)\n*          On entry, the first (N-1) entries contain the subdiagonal\n*          elements of the tridiagonal matrix T; E(N) is set to ZERO.\n*\n*  INFO    (output) INTEGER\n*          INFO = 0(default) : the matrix warrants computations preserving\n*                              relative accuracy.\n*          INFO = 1          : the matrix warrants computations guaranteeing\n*                              only absolute accuracy.\n*\n\n*  Further Details\n*  ===============\n*\n*  Based on contributions by\n*     Beresford Parlett, University of California, Berkeley, USA\n*     Jim Demmel, University of California, Berkeley, USA\n*     Inderjit Dhillon, University of Texas, Austin, USA\n*     Osni Marques, LBNL/NERSC, USA\n*     Christof Voemel, University of California, Berkeley, USA\n*\n*  =====================================================================\n*\n\n");
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
    rb_e_out__ = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  e_out__ = NA_PTR_TYPE(rb_e_out__, doublereal*);
  MEMCPY(e_out__, e, doublereal, NA_TOTAL(rb_e));
  rb_e = rb_e_out__;
  e = e_out__;

  dlarrr_(&n, d, e, &info);

  rb_info = INT2NUM(info);
  return rb_ary_new3(2, rb_info, rb_e);
}

void
init_lapack_dlarrr(VALUE mLapack){
  rb_define_module_function(mLapack, "dlarrr", rb_dlarrr, -1);
}
