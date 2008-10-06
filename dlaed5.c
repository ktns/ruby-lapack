#include "rb_lapack.h"

static VALUE
rb_dlaed5(int argc, VALUE *argv, VALUE self){
  VALUE rb_i;
  integer i; 
  VALUE rb_d;
  doublereal *d; 
  VALUE rb_z;
  doublereal *z; 
  VALUE rb_rho;
  doublereal rho; 
  VALUE rb_delta;
  doublereal *delta; 
  VALUE rb_dlam;
  doublereal dlam; 


  if (argc == 0) {
    printf("%s\n", "USAGE:\n  delta, dlam = NumRu::Lapack.dlaed5( i, d, z, rho)\n    or\n  NumRu::Lapack.dlaed5  # print help\n\n\nFORTRAN MANUAL\n      SUBROUTINE DLAED5( I, D, Z, DELTA, RHO, DLAM )\n\n*  Purpose\n*  =======\n*\n*  This subroutine computes the I-th eigenvalue of a symmetric rank-one\n*  modification of a 2-by-2 diagonal matrix\n*\n*             diag( D )  +  RHO *  Z * transpose(Z) .\n*\n*  The diagonal elements in the array D are assumed to satisfy\n*\n*             D(i) < D(j)  for  i < j .\n*\n*  We also assume RHO > 0 and that the Euclidean norm of the vector\n*  Z is one.\n*\n\n*  Arguments\n*  =========\n*\n*  I      (input) INTEGER\n*         The index of the eigenvalue to be computed.  I = 1 or I = 2.\n*\n*  D      (input) DOUBLE PRECISION array, dimension (2)\n*         The original eigenvalues.  We assume D(1) < D(2).\n*\n*  Z      (input) DOUBLE PRECISION array, dimension (2)\n*         The components of the updating vector.\n*\n*  DELTA  (output) DOUBLE PRECISION array, dimension (2)\n*         The vector DELTA contains the information necessary\n*         to construct the eigenvectors.\n*\n*  RHO    (input) DOUBLE PRECISION\n*         The scalar in the symmetric updating formula.\n*\n*  DLAM   (output) DOUBLE PRECISION\n*         The computed lambda_I, the I-th updated eigenvalue.\n*\n\n*  Further Details\n*  ===============\n*\n*  Based on contributions by\n*     Ren-Cang Li, Computer Science Division, University of California\n*     at Berkeley, USA\n*\n*  =====================================================================\n*\n\n");
    return Qnil;
  }
  if (argc != 4)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 4)", argc);
  rb_i = argv[0];
  rb_d = argv[1];
  rb_z = argv[2];
  rb_rho = argv[3];

  i = NUM2INT(rb_i);
  rho = NUM2DBL(rb_rho);
  if (!NA_IsNArray(rb_d))
    rb_raise(rb_eArgError, "d (2th argument) must be NArray");
  if (NA_RANK(rb_d) != 1)
    rb_raise(rb_eArgError, "rank of d (2th argument) must be %d", 1);
  if (NA_SHAPE0(rb_d) != (2))
    rb_raise(rb_eRuntimeError, "shape 0 of d must be %d", 2);
  if (NA_TYPE(rb_d) != NA_DFLOAT)
    rb_d = na_change_type(rb_d, NA_DFLOAT);
  d = NA_PTR_TYPE(rb_d, doublereal*);
  if (!NA_IsNArray(rb_z))
    rb_raise(rb_eArgError, "z (3th argument) must be NArray");
  if (NA_RANK(rb_z) != 1)
    rb_raise(rb_eArgError, "rank of z (3th argument) must be %d", 1);
  if (NA_SHAPE0(rb_z) != (2))
    rb_raise(rb_eRuntimeError, "shape 0 of z must be %d", 2);
  if (NA_TYPE(rb_z) != NA_DFLOAT)
    rb_z = na_change_type(rb_z, NA_DFLOAT);
  z = NA_PTR_TYPE(rb_z, doublereal*);
  {
    int shape[1];
    shape[0] = 2;
    rb_delta = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  delta = NA_PTR_TYPE(rb_delta, doublereal*);

  dlaed5_(&i, d, z, delta, &rho, &dlam);

  rb_dlam = rb_float_new((double)dlam);
  return rb_ary_new3(2, rb_delta, rb_dlam);
}

void
init_lapack_dlaed5(VALUE mLapack){
  rb_define_module_function(mLapack, "dlaed5", rb_dlaed5, -1);
}
