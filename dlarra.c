#include "rb_lapack.h"

static VALUE
rb_dlarra(int argc, VALUE *argv, VALUE self){
  VALUE rb_d;
  doublereal *d; 
  VALUE rb_e;
  doublereal *e; 
  VALUE rb_e2;
  doublereal *e2; 
  VALUE rb_spltol;
  doublereal spltol; 
  VALUE rb_tnrm;
  doublereal tnrm; 
  VALUE rb_nsplit;
  integer nsplit; 
  VALUE rb_isplit;
  integer *isplit; 
  VALUE rb_info;
  integer info; 
  VALUE rb_e_out__;
  doublereal *e_out__;
  VALUE rb_e2_out__;
  doublereal *e2_out__;

  integer n;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  nsplit, isplit, info, e, e2 = NumRu::Lapack.dlarra( d, e, e2, spltol, tnrm)\n    or\n  NumRu::Lapack.dlarra  # print help\n\n\nFORTRAN MANUAL\n      SUBROUTINE DLARRA( N, D, E, E2, SPLTOL, TNRM, NSPLIT, ISPLIT, INFO )\n\n*  Purpose\n*  =======\n*\n*  Compute the splitting points with threshold SPLTOL.\n*  DLARRA sets any \"small\" off-diagonal elements to zero.\n*\n\n*  Arguments\n*  =========\n*\n*  N       (input) INTEGER\n*          The order of the matrix. N > 0.\n*\n*  D       (input) DOUBLE PRECISION array, dimension (N)\n*          On entry, the N diagonal elements of the tridiagonal\n*          matrix T.\n*\n*  E       (input/output) DOUBLE PRECISION array, dimension (N)\n*          On entry, the first (N-1) entries contain the subdiagonal\n*          elements of the tridiagonal matrix T; E(N) need not be set.\n*          On exit, the entries E( ISPLIT( I ) ), 1 <= I <= NSPLIT,\n*          are set to zero, the other entries of E are untouched.\n*\n*  E2      (input/output) DOUBLE PRECISION array, dimension (N)\n*          On entry, the first (N-1) entries contain the SQUARES of the\n*          subdiagonal elements of the tridiagonal matrix T;\n*          E2(N) need not be set.\n*          On exit, the entries E2( ISPLIT( I ) ),\n*          1 <= I <= NSPLIT, have been set to zero\n*\n*  SPLTOL (input) DOUBLE PRECISION\n*          The threshold for splitting. Two criteria can be used:\n*          SPLTOL<0 : criterion based on absolute off-diagonal value\n*          SPLTOL>0 : criterion that preserves relative accuracy\n*\n*  TNRM (input) DOUBLE PRECISION\n*          The norm of the matrix.\n*\n*  NSPLIT  (output) INTEGER\n*          The number of blocks T splits into. 1 <= NSPLIT <= N.\n*\n*  ISPLIT  (output) INTEGER array, dimension (N)\n*          The splitting points, at which T breaks up into blocks.\n*          The first block consists of rows/columns 1 to ISPLIT(1),\n*          the second of rows/columns ISPLIT(1)+1 through ISPLIT(2),\n*          etc., and the NSPLIT-th consists of rows/columns\n*          ISPLIT(NSPLIT-1)+1 through ISPLIT(NSPLIT)=N.\n*\n*\n*  INFO    (output) INTEGER\n*          = 0:  successful exit\n*\n\n*  Further Details\n*  ===============\n*\n*  Based on contributions by\n*     Beresford Parlett, University of California, Berkeley, USA\n*     Jim Demmel, University of California, Berkeley, USA\n*     Inderjit Dhillon, University of Texas, Austin, USA\n*     Osni Marques, LBNL/NERSC, USA\n*     Christof Voemel, University of California, Berkeley, USA\n*\n*  =====================================================================\n*\n\n");
    return Qnil;
  }
  if (argc != 5)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 5)", argc);
  rb_d = argv[0];
  rb_e = argv[1];
  rb_e2 = argv[2];
  rb_spltol = argv[3];
  rb_tnrm = argv[4];

  spltol = NUM2DBL(rb_spltol);
  tnrm = NUM2DBL(rb_tnrm);
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
  if (!NA_IsNArray(rb_e2))
    rb_raise(rb_eArgError, "e2 (3th argument) must be NArray");
  if (NA_RANK(rb_e2) != 1)
    rb_raise(rb_eArgError, "rank of e2 (3th argument) must be %d", 1);
  if (NA_SHAPE0(rb_e2) != n)
    rb_raise(rb_eRuntimeError, "shape 0 of e2 must be the same as shape 0 of d");
  if (NA_TYPE(rb_e2) != NA_DFLOAT)
    rb_e2 = na_change_type(rb_e2, NA_DFLOAT);
  e2 = NA_PTR_TYPE(rb_e2, doublereal*);
  {
    int shape[1];
    shape[0] = DIM_LEN(n);
    rb_isplit = na_make_object(NA_LINT, 1, shape, cNArray);
  }
  isplit = NA_PTR_TYPE(rb_isplit, integer*);
  {
    int shape[1];
    shape[0] = n;
    rb_e_out__ = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  e_out__ = NA_PTR_TYPE(rb_e_out__, doublereal*);
  MEMCPY(e_out__, e, doublereal, NA_TOTAL(rb_e));
  rb_e = rb_e_out__;
  e = e_out__;
  {
    int shape[1];
    shape[0] = n;
    rb_e2_out__ = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  e2_out__ = NA_PTR_TYPE(rb_e2_out__, doublereal*);
  MEMCPY(e2_out__, e2, doublereal, NA_TOTAL(rb_e2));
  rb_e2 = rb_e2_out__;
  e2 = e2_out__;

  dlarra_(&n, d, e, e2, &spltol, &tnrm, &nsplit, isplit, &info);

  rb_nsplit = INT2NUM(nsplit);
  rb_info = INT2NUM(info);
  return rb_ary_new3(5, rb_nsplit, rb_isplit, rb_info, rb_e, rb_e2);
}

void
init_lapack_dlarra(VALUE mLapack){
  rb_define_module_function(mLapack, "dlarra", rb_dlarra, -1);
}
