#include "rb_lapack.h"

extern VOID dlaed5_(integer *i, doublereal *d, doublereal *z, doublereal *delta, doublereal *rho, doublereal *dlam);

static VALUE sHelp, sUsage;

static VALUE
rblapack_dlaed5(int argc, VALUE *argv, VALUE self){
  VALUE rblapack_i;
  integer i; 
  VALUE rblapack_d;
  doublereal *d; 
  VALUE rblapack_z;
  doublereal *z; 
  VALUE rblapack_rho;
  doublereal rho; 
  VALUE rblapack_delta;
  doublereal *delta; 
  VALUE rblapack_dlam;
  doublereal dlam; 


  VALUE rb_options;
  if (argc > 0 && TYPE(argv[argc-1]) == T_HASH) {
    argc--;
    rb_options = argv[argc];
    if (rb_hash_aref(rb_options, sHelp) == Qtrue) {
      printf("%s\n", "USAGE:\n  delta, dlam = NumRu::Lapack.dlaed5( i, d, z, rho, [:usage => usage, :help => help])\n\n\nFORTRAN MANUAL\n      SUBROUTINE DLAED5( I, D, Z, DELTA, RHO, DLAM )\n\n*  Purpose\n*  =======\n*\n*  This subroutine computes the I-th eigenvalue of a symmetric rank-one\n*  modification of a 2-by-2 diagonal matrix\n*\n*             diag( D )  +  RHO *  Z * transpose(Z) .\n*\n*  The diagonal elements in the array D are assumed to satisfy\n*\n*             D(i) < D(j)  for  i < j .\n*\n*  We also assume RHO > 0 and that the Euclidean norm of the vector\n*  Z is one.\n*\n\n*  Arguments\n*  =========\n*\n*  I      (input) INTEGER\n*         The index of the eigenvalue to be computed.  I = 1 or I = 2.\n*\n*  D      (input) DOUBLE PRECISION array, dimension (2)\n*         The original eigenvalues.  We assume D(1) < D(2).\n*\n*  Z      (input) DOUBLE PRECISION array, dimension (2)\n*         The components of the updating vector.\n*\n*  DELTA  (output) DOUBLE PRECISION array, dimension (2)\n*         The vector DELTA contains the information necessary\n*         to construct the eigenvectors.\n*\n*  RHO    (input) DOUBLE PRECISION\n*         The scalar in the symmetric updating formula.\n*\n*  DLAM   (output) DOUBLE PRECISION\n*         The computed lambda_I, the I-th updated eigenvalue.\n*\n\n*  Further Details\n*  ===============\n*\n*  Based on contributions by\n*     Ren-Cang Li, Computer Science Division, University of California\n*     at Berkeley, USA\n*\n*  =====================================================================\n*\n\n");
      rb_exit(0);
    }
    if (rb_hash_aref(rb_options, sUsage) == Qtrue) {
      printf("%s\n", "USAGE:\n  delta, dlam = NumRu::Lapack.dlaed5( i, d, z, rho, [:usage => usage, :help => help])\n");
      rb_exit(0);
    } 
  } else
    rb_options = Qnil;
  if (argc != 4)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 4)", argc);
  rblapack_i = argv[0];
  rblapack_d = argv[1];
  rblapack_z = argv[2];
  rblapack_rho = argv[3];
  if (rb_options != Qnil) {
  }

  rho = NUM2DBL(rblapack_rho);
  if (!NA_IsNArray(rblapack_z))
    rb_raise(rb_eArgError, "z (3th argument) must be NArray");
  if (NA_RANK(rblapack_z) != 1)
    rb_raise(rb_eArgError, "rank of z (3th argument) must be %d", 1);
  if (NA_SHAPE0(rblapack_z) != (2))
    rb_raise(rb_eRuntimeError, "shape 0 of z must be %d", 2);
  if (NA_TYPE(rblapack_z) != NA_DFLOAT)
    rblapack_z = na_change_type(rblapack_z, NA_DFLOAT);
  z = NA_PTR_TYPE(rblapack_z, doublereal*);
  if (!NA_IsNArray(rblapack_d))
    rb_raise(rb_eArgError, "d (2th argument) must be NArray");
  if (NA_RANK(rblapack_d) != 1)
    rb_raise(rb_eArgError, "rank of d (2th argument) must be %d", 1);
  if (NA_SHAPE0(rblapack_d) != (2))
    rb_raise(rb_eRuntimeError, "shape 0 of d must be %d", 2);
  if (NA_TYPE(rblapack_d) != NA_DFLOAT)
    rblapack_d = na_change_type(rblapack_d, NA_DFLOAT);
  d = NA_PTR_TYPE(rblapack_d, doublereal*);
  i = NUM2INT(rblapack_i);
  {
    int shape[1];
    shape[0] = 2;
    rblapack_delta = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  delta = NA_PTR_TYPE(rblapack_delta, doublereal*);

  dlaed5_(&i, d, z, delta, &rho, &dlam);

  rblapack_dlam = rb_float_new((double)dlam);
  return rb_ary_new3(2, rblapack_delta, rblapack_dlam);
}

void
init_lapack_dlaed5(VALUE mLapack){
  rb_define_module_function(mLapack, "dlaed5", rblapack_dlaed5, -1);
  sHelp = ID2SYM(rb_intern("help"));
  sUsage = ID2SYM(rb_intern("usage"));
}
