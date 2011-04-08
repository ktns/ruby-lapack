#include "rb_lapack.h"

extern VOID dsterf_(integer *n, doublereal *d, doublereal *e, integer *info);

static VALUE sHelp, sUsage;

static VALUE
rblapack_dsterf(int argc, VALUE *argv, VALUE self){
  VALUE rblapack_d;
  doublereal *d; 
  VALUE rblapack_e;
  doublereal *e; 
  VALUE rblapack_info;
  integer info; 
  VALUE rblapack_d_out__;
  doublereal *d_out__;
  VALUE rblapack_e_out__;
  doublereal *e_out__;

  integer n;

  VALUE rb_options;
  if (argc > 0 && TYPE(argv[argc-1]) == T_HASH) {
    argc--;
    rb_options = argv[argc];
    if (rb_hash_aref(rb_options, sHelp) == Qtrue) {
      printf("%s\n", "USAGE:\n  info, d, e = NumRu::Lapack.dsterf( d, e, [:usage => usage, :help => help])\n\n\nFORTRAN MANUAL\n      SUBROUTINE DSTERF( N, D, E, INFO )\n\n*  Purpose\n*  =======\n*\n*  DSTERF computes all eigenvalues of a symmetric tridiagonal matrix\n*  using the Pal-Walker-Kahan variant of the QL or QR algorithm.\n*\n\n*  Arguments\n*  =========\n*\n*  N       (input) INTEGER\n*          The order of the matrix.  N >= 0.\n*\n*  D       (input/output) DOUBLE PRECISION array, dimension (N)\n*          On entry, the n diagonal elements of the tridiagonal matrix.\n*          On exit, if INFO = 0, the eigenvalues in ascending order.\n*\n*  E       (input/output) DOUBLE PRECISION array, dimension (N-1)\n*          On entry, the (n-1) subdiagonal elements of the tridiagonal\n*          matrix.\n*          On exit, E has been destroyed.\n*\n*  INFO    (output) INTEGER\n*          = 0:  successful exit\n*          < 0:  if INFO = -i, the i-th argument had an illegal value\n*          > 0:  the algorithm failed to find all of the eigenvalues in\n*                a total of 30*N iterations; if INFO = i, then i\n*                elements of E have not converged to zero.\n*\n\n*  =====================================================================\n*\n\n");
      rb_exit(0);
    }
    if (rb_hash_aref(rb_options, sUsage) == Qtrue) {
      printf("%s\n", "USAGE:\n  info, d, e = NumRu::Lapack.dsterf( d, e, [:usage => usage, :help => help])\n");
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
  if (NA_TYPE(rblapack_d) != NA_DFLOAT)
    rblapack_d = na_change_type(rblapack_d, NA_DFLOAT);
  d = NA_PTR_TYPE(rblapack_d, doublereal*);
  if (!NA_IsNArray(rblapack_e))
    rb_raise(rb_eArgError, "e (2th argument) must be NArray");
  if (NA_RANK(rblapack_e) != 1)
    rb_raise(rb_eArgError, "rank of e (2th argument) must be %d", 1);
  if (NA_SHAPE0(rblapack_e) != (n-1))
    rb_raise(rb_eRuntimeError, "shape 0 of e must be %d", n-1);
  if (NA_TYPE(rblapack_e) != NA_DFLOAT)
    rblapack_e = na_change_type(rblapack_e, NA_DFLOAT);
  e = NA_PTR_TYPE(rblapack_e, doublereal*);
  {
    int shape[1];
    shape[0] = n;
    rblapack_d_out__ = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  d_out__ = NA_PTR_TYPE(rblapack_d_out__, doublereal*);
  MEMCPY(d_out__, d, doublereal, NA_TOTAL(rblapack_d));
  rblapack_d = rblapack_d_out__;
  d = d_out__;
  {
    int shape[1];
    shape[0] = n-1;
    rblapack_e_out__ = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  e_out__ = NA_PTR_TYPE(rblapack_e_out__, doublereal*);
  MEMCPY(e_out__, e, doublereal, NA_TOTAL(rblapack_e));
  rblapack_e = rblapack_e_out__;
  e = e_out__;

  dsterf_(&n, d, e, &info);

  rblapack_info = INT2NUM(info);
  return rb_ary_new3(3, rblapack_info, rblapack_d, rblapack_e);
}

void
init_lapack_dsterf(VALUE mLapack){
  rb_define_module_function(mLapack, "dsterf", rblapack_dsterf, -1);
  sHelp = ID2SYM(rb_intern("help"));
  sUsage = ID2SYM(rb_intern("usage"));
}
