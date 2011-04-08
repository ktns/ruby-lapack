#include "rb_lapack.h"

extern VOID dptts2_(integer *n, integer *nrhs, doublereal *d, doublereal *e, doublereal *b, integer *ldb);

static VALUE sHelp, sUsage;

static VALUE
rblapack_dptts2(int argc, VALUE *argv, VALUE self){
  VALUE rblapack_d;
  doublereal *d; 
  VALUE rblapack_e;
  doublereal *e; 
  VALUE rblapack_b;
  doublereal *b; 
  VALUE rblapack_b_out__;
  doublereal *b_out__;

  integer n;
  integer ldb;
  integer nrhs;

  VALUE rb_options;
  if (argc > 0 && TYPE(argv[argc-1]) == T_HASH) {
    argc--;
    rb_options = argv[argc];
    if (rb_hash_aref(rb_options, sHelp) == Qtrue) {
      printf("%s\n", "USAGE:\n  b = NumRu::Lapack.dptts2( d, e, b, [:usage => usage, :help => help])\n\n\nFORTRAN MANUAL\n      SUBROUTINE DPTTS2( N, NRHS, D, E, B, LDB )\n\n*  Purpose\n*  =======\n*\n*  DPTTS2 solves a tridiagonal system of the form\n*     A * X = B\n*  using the L*D*L' factorization of A computed by DPTTRF.  D is a\n*  diagonal matrix specified in the vector D, L is a unit bidiagonal\n*  matrix whose subdiagonal is specified in the vector E, and X and B\n*  are N by NRHS matrices.\n*\n\n*  Arguments\n*  =========\n*\n*  N       (input) INTEGER\n*          The order of the tridiagonal matrix A.  N >= 0.\n*\n*  NRHS    (input) INTEGER\n*          The number of right hand sides, i.e., the number of columns\n*          of the matrix B.  NRHS >= 0.\n*\n*  D       (input) DOUBLE PRECISION array, dimension (N)\n*          The n diagonal elements of the diagonal matrix D from the\n*          L*D*L' factorization of A.\n*\n*  E       (input) DOUBLE PRECISION array, dimension (N-1)\n*          The (n-1) subdiagonal elements of the unit bidiagonal factor\n*          L from the L*D*L' factorization of A.  E can also be regarded\n*          as the superdiagonal of the unit bidiagonal factor U from the\n*          factorization A = U'*D*U.\n*\n*  B       (input/output) DOUBLE PRECISION array, dimension (LDB,NRHS)\n*          On entry, the right hand side vectors B for the system of\n*          linear equations.\n*          On exit, the solution vectors, X.\n*\n*  LDB     (input) INTEGER\n*          The leading dimension of the array B.  LDB >= max(1,N).\n*\n\n*  =====================================================================\n*\n*     .. Local Scalars ..\n      INTEGER            I, J\n*     ..\n*     .. External Subroutines ..\n      EXTERNAL           DSCAL\n*     ..\n\n");
      rb_exit(0);
    }
    if (rb_hash_aref(rb_options, sUsage) == Qtrue) {
      printf("%s\n", "USAGE:\n  b = NumRu::Lapack.dptts2( d, e, b, [:usage => usage, :help => help])\n");
      rb_exit(0);
    } 
  } else
    rb_options = Qnil;
  if (argc != 3)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 3)", argc);
  rblapack_d = argv[0];
  rblapack_e = argv[1];
  rblapack_b = argv[2];
  if (rb_options != Qnil) {
  }

  if (!NA_IsNArray(rblapack_b))
    rb_raise(rb_eArgError, "b (3th argument) must be NArray");
  if (NA_RANK(rblapack_b) != 2)
    rb_raise(rb_eArgError, "rank of b (3th argument) must be %d", 2);
  nrhs = NA_SHAPE1(rblapack_b);
  ldb = NA_SHAPE0(rblapack_b);
  if (NA_TYPE(rblapack_b) != NA_DFLOAT)
    rblapack_b = na_change_type(rblapack_b, NA_DFLOAT);
  b = NA_PTR_TYPE(rblapack_b, doublereal*);
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
    int shape[2];
    shape[0] = ldb;
    shape[1] = nrhs;
    rblapack_b_out__ = na_make_object(NA_DFLOAT, 2, shape, cNArray);
  }
  b_out__ = NA_PTR_TYPE(rblapack_b_out__, doublereal*);
  MEMCPY(b_out__, b, doublereal, NA_TOTAL(rblapack_b));
  rblapack_b = rblapack_b_out__;
  b = b_out__;

  dptts2_(&n, &nrhs, d, e, b, &ldb);

  return rblapack_b;
}

void
init_lapack_dptts2(VALUE mLapack){
  rb_define_module_function(mLapack, "dptts2", rblapack_dptts2, -1);
  sHelp = ID2SYM(rb_intern("help"));
  sUsage = ID2SYM(rb_intern("usage"));
}
