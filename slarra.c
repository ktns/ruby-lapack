#include "rb_lapack.h"

static VALUE
rb_slarra(int argc, VALUE *argv, VALUE self){
  VALUE rb_d;
  real *d; 
  VALUE rb_e;
  real *e; 
  VALUE rb_e2;
  real *e2; 
  VALUE rb_spltol;
  real spltol; 
  VALUE rb_tnrm;
  real tnrm; 
  VALUE rb_nsplit;
  integer nsplit; 
  VALUE rb_isplit;
  integer *isplit; 
  VALUE rb_info;
  integer info; 
  VALUE rb_e_out__;
  real *e_out__;
  VALUE rb_e2_out__;
  real *e2_out__;

  integer n;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  nsplit, isplit, info, e, e2 = NumRu::Lapack.slarra( d, e, e2, spltol, tnrm)\n    or\n  NumRu::Lapack.slarra  # print help\n\n\nFORTRAN MANUAL\n      SUBROUTINE SLARRA( N, D, E, E2, SPLTOL, TNRM, NSPLIT, ISPLIT, INFO )\n\n*  Purpose\n*  =======\n*\n*  Compute the splitting points with threshold SPLTOL.\n*  SLARRA sets any \"small\" off-diagonal elements to zero.\n*\n\n*  Arguments\n*  =========\n*\n*  N       (input) INTEGER\n*          The order of the matrix. N > 0.\n*\n*  D       (input) REAL             array, dimension (N)\n*          On entry, the N diagonal elements of the tridiagonal\n*          matrix T.\n*\n*  E       (input/output) REAL             array, dimension (N)\n*          On entry, the first (N-1) entries contain the subdiagonal\n*          elements of the tridiagonal matrix T; E(N) need not be set.\n*          On exit, the entries E( ISPLIT( I ) ), 1 <= I <= NSPLIT,\n*          are set to zero, the other entries of E are untouched.\n*\n*  E2      (input/output) REAL             array, dimension (N)\n*          On entry, the first (N-1) entries contain the SQUARES of the\n*          subdiagonal elements of the tridiagonal matrix T;\n*          E2(N) need not be set.\n*          On exit, the entries E2( ISPLIT( I ) ),\n*          1 <= I <= NSPLIT, have been set to zero\n*\n*  SPLTOL (input) REAL            \n*          The threshold for splitting. Two criteria can be used:\n*          SPLTOL<0 : criterion based on absolute off-diagonal value\n*          SPLTOL>0 : criterion that preserves relative accuracy\n*\n*  TNRM (input) REAL            \n*          The norm of the matrix.\n*\n*  NSPLIT  (output) INTEGER\n*          The number of blocks T splits into. 1 <= NSPLIT <= N.\n*\n*  ISPLIT  (output) INTEGER array, dimension (N)\n*          The splitting points, at which T breaks up into blocks.\n*          The first block consists of rows/columns 1 to ISPLIT(1),\n*          the second of rows/columns ISPLIT(1)+1 through ISPLIT(2),\n*          etc., and the NSPLIT-th consists of rows/columns\n*          ISPLIT(NSPLIT-1)+1 through ISPLIT(NSPLIT)=N.\n*\n*\n*  INFO    (output) INTEGER\n*          = 0:  successful exit\n*\n\n*  Further Details\n*  ===============\n*\n*  Based on contributions by\n*     Beresford Parlett, University of California, Berkeley, USA\n*     Jim Demmel, University of California, Berkeley, USA\n*     Inderjit Dhillon, University of Texas, Austin, USA\n*     Osni Marques, LBNL/NERSC, USA\n*     Christof Voemel, University of California, Berkeley, USA\n*\n*  =====================================================================\n*\n\n");
    return Qnil;
  }
  if (argc != 5)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 5)", argc);
  rb_d = argv[0];
  rb_e = argv[1];
  rb_e2 = argv[2];
  rb_spltol = argv[3];
  rb_tnrm = argv[4];

  spltol = (real)NUM2DBL(rb_spltol);
  tnrm = (real)NUM2DBL(rb_tnrm);
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
  if (!NA_IsNArray(rb_e2))
    rb_raise(rb_eArgError, "e2 (3th argument) must be NArray");
  if (NA_RANK(rb_e2) != 1)
    rb_raise(rb_eArgError, "rank of e2 (3th argument) must be %d", 1);
  if (NA_SHAPE0(rb_e2) != n)
    rb_raise(rb_eRuntimeError, "shape 0 of e2 must be the same as shape 0 of d");
  if (NA_TYPE(rb_e2) != NA_SFLOAT)
    rb_e2 = na_change_type(rb_e2, NA_SFLOAT);
  e2 = NA_PTR_TYPE(rb_e2, real*);
  {
    int shape[1];
    shape[0] = DIM_LEN(n);
    rb_isplit = na_make_object(NA_LINT, 1, shape, cNArray);
  }
  isplit = NA_PTR_TYPE(rb_isplit, integer*);
  {
    int shape[1];
    shape[0] = n;
    rb_e_out__ = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  e_out__ = NA_PTR_TYPE(rb_e_out__, real*);
  MEMCPY(e_out__, e, real, NA_TOTAL(rb_e));
  rb_e = rb_e_out__;
  e = e_out__;
  {
    int shape[1];
    shape[0] = n;
    rb_e2_out__ = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  e2_out__ = NA_PTR_TYPE(rb_e2_out__, real*);
  MEMCPY(e2_out__, e2, real, NA_TOTAL(rb_e2));
  rb_e2 = rb_e2_out__;
  e2 = e2_out__;

  slarra_(&n, d, e, e2, &spltol, &tnrm, &nsplit, isplit, &info);

  rb_nsplit = INT2NUM(nsplit);
  rb_info = INT2NUM(info);
  return rb_ary_new3(5, rb_nsplit, rb_isplit, rb_info, rb_e, rb_e2);
}

void
init_lapack_slarra(VALUE mLapack){
  rb_define_module_function(mLapack, "slarra", rb_slarra, -1);
}
