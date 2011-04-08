#include "rb_lapack.h"

extern VOID cptts2_(integer *iuplo, integer *n, integer *nrhs, real *d, complex *e, complex *b, integer *ldb);

static VALUE sHelp, sUsage;

static VALUE
rblapack_cptts2(int argc, VALUE *argv, VALUE self){
  VALUE rblapack_iuplo;
  integer iuplo; 
  VALUE rblapack_d;
  real *d; 
  VALUE rblapack_e;
  complex *e; 
  VALUE rblapack_b;
  complex *b; 
  VALUE rblapack_b_out__;
  complex *b_out__;

  integer n;
  integer ldb;
  integer nrhs;

  VALUE rb_options;
  if (argc > 0 && TYPE(argv[argc-1]) == T_HASH) {
    argc--;
    rb_options = argv[argc];
    if (rb_hash_aref(rb_options, sHelp) == Qtrue) {
      printf("%s\n", "USAGE:\n  b = NumRu::Lapack.cptts2( iuplo, d, e, b, [:usage => usage, :help => help])\n\n\nFORTRAN MANUAL\n      SUBROUTINE CPTTS2( IUPLO, N, NRHS, D, E, B, LDB )\n\n*  Purpose\n*  =======\n*\n*  CPTTS2 solves a tridiagonal system of the form\n*     A * X = B\n*  using the factorization A = U'*D*U or A = L*D*L' computed by CPTTRF.\n*  D is a diagonal matrix specified in the vector D, U (or L) is a unit\n*  bidiagonal matrix whose superdiagonal (subdiagonal) is specified in\n*  the vector E, and X and B are N by NRHS matrices.\n*\n\n*  Arguments\n*  =========\n*\n*  IUPLO   (input) INTEGER\n*          Specifies the form of the factorization and whether the\n*          vector E is the superdiagonal of the upper bidiagonal factor\n*          U or the subdiagonal of the lower bidiagonal factor L.\n*          = 1:  A = U'*D*U, E is the superdiagonal of U\n*          = 0:  A = L*D*L', E is the subdiagonal of L\n*\n*  N       (input) INTEGER\n*          The order of the tridiagonal matrix A.  N >= 0.\n*\n*  NRHS    (input) INTEGER\n*          The number of right hand sides, i.e., the number of columns\n*          of the matrix B.  NRHS >= 0.\n*\n*  D       (input) REAL array, dimension (N)\n*          The n diagonal elements of the diagonal matrix D from the\n*          factorization A = U'*D*U or A = L*D*L'.\n*\n*  E       (input) COMPLEX array, dimension (N-1)\n*          If IUPLO = 1, the (n-1) superdiagonal elements of the unit\n*          bidiagonal factor U from the factorization A = U'*D*U.\n*          If IUPLO = 0, the (n-1) subdiagonal elements of the unit\n*          bidiagonal factor L from the factorization A = L*D*L'.\n*\n*  B       (input/output) REAL array, dimension (LDB,NRHS)\n*          On entry, the right hand side vectors B for the system of\n*          linear equations.\n*          On exit, the solution vectors, X.\n*\n*  LDB     (input) INTEGER\n*          The leading dimension of the array B.  LDB >= max(1,N).\n*\n\n*  =====================================================================\n*\n*     .. Local Scalars ..\n      INTEGER            I, J\n*     ..\n*     .. External Subroutines ..\n      EXTERNAL           CSSCAL\n*     ..\n*     .. Intrinsic Functions ..\n      INTRINSIC          CONJG\n*     ..\n\n");
      rb_exit(0);
    }
    if (rb_hash_aref(rb_options, sUsage) == Qtrue) {
      printf("%s\n", "USAGE:\n  b = NumRu::Lapack.cptts2( iuplo, d, e, b, [:usage => usage, :help => help])\n");
      rb_exit(0);
    } 
  } else
    rb_options = Qnil;
  if (argc != 4)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 4)", argc);
  rblapack_iuplo = argv[0];
  rblapack_d = argv[1];
  rblapack_e = argv[2];
  rblapack_b = argv[3];
  if (rb_options != Qnil) {
  }

  if (!NA_IsNArray(rblapack_b))
    rb_raise(rb_eArgError, "b (4th argument) must be NArray");
  if (NA_RANK(rblapack_b) != 2)
    rb_raise(rb_eArgError, "rank of b (4th argument) must be %d", 2);
  nrhs = NA_SHAPE1(rblapack_b);
  ldb = NA_SHAPE0(rblapack_b);
  if (NA_TYPE(rblapack_b) != NA_SCOMPLEX)
    rblapack_b = na_change_type(rblapack_b, NA_SCOMPLEX);
  b = NA_PTR_TYPE(rblapack_b, complex*);
  if (!NA_IsNArray(rblapack_d))
    rb_raise(rb_eArgError, "d (2th argument) must be NArray");
  if (NA_RANK(rblapack_d) != 1)
    rb_raise(rb_eArgError, "rank of d (2th argument) must be %d", 1);
  n = NA_SHAPE0(rblapack_d);
  if (NA_TYPE(rblapack_d) != NA_SFLOAT)
    rblapack_d = na_change_type(rblapack_d, NA_SFLOAT);
  d = NA_PTR_TYPE(rblapack_d, real*);
  iuplo = NUM2INT(rblapack_iuplo);
  if (!NA_IsNArray(rblapack_e))
    rb_raise(rb_eArgError, "e (3th argument) must be NArray");
  if (NA_RANK(rblapack_e) != 1)
    rb_raise(rb_eArgError, "rank of e (3th argument) must be %d", 1);
  if (NA_SHAPE0(rblapack_e) != (n-1))
    rb_raise(rb_eRuntimeError, "shape 0 of e must be %d", n-1);
  if (NA_TYPE(rblapack_e) != NA_SCOMPLEX)
    rblapack_e = na_change_type(rblapack_e, NA_SCOMPLEX);
  e = NA_PTR_TYPE(rblapack_e, complex*);
  {
    int shape[2];
    shape[0] = ldb;
    shape[1] = nrhs;
    rblapack_b_out__ = na_make_object(NA_SCOMPLEX, 2, shape, cNArray);
  }
  b_out__ = NA_PTR_TYPE(rblapack_b_out__, complex*);
  MEMCPY(b_out__, b, complex, NA_TOTAL(rblapack_b));
  rblapack_b = rblapack_b_out__;
  b = b_out__;

  cptts2_(&iuplo, &n, &nrhs, d, e, b, &ldb);

  return rblapack_b;
}

void
init_lapack_cptts2(VALUE mLapack){
  rb_define_module_function(mLapack, "cptts2", rblapack_cptts2, -1);
  sHelp = ID2SYM(rb_intern("help"));
  sUsage = ID2SYM(rb_intern("usage"));
}
