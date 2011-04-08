#include "rb_lapack.h"

extern VOID cgtts2_(integer *itrans, integer *n, integer *nrhs, complex *dl, complex *d, complex *du, complex *du2, integer *ipiv, complex *b, integer *ldb);

static VALUE sHelp, sUsage;

static VALUE
rblapack_cgtts2(int argc, VALUE *argv, VALUE self){
  VALUE rblapack_itrans;
  integer itrans; 
  VALUE rblapack_dl;
  complex *dl; 
  VALUE rblapack_d;
  complex *d; 
  VALUE rblapack_du;
  complex *du; 
  VALUE rblapack_du2;
  complex *du2; 
  VALUE rblapack_ipiv;
  integer *ipiv; 
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
      printf("%s\n", "USAGE:\n  b = NumRu::Lapack.cgtts2( itrans, dl, d, du, du2, ipiv, b, [:usage => usage, :help => help])\n\n\nFORTRAN MANUAL\n      SUBROUTINE CGTTS2( ITRANS, N, NRHS, DL, D, DU, DU2, IPIV, B, LDB )\n\n*  Purpose\n*  =======\n*\n*  CGTTS2 solves one of the systems of equations\n*     A * X = B,  A**T * X = B,  or  A**H * X = B,\n*  with a tridiagonal matrix A using the LU factorization computed\n*  by CGTTRF.\n*\n\n*  Arguments\n*  =========\n*\n*  ITRANS  (input) INTEGER\n*          Specifies the form of the system of equations.\n*          = 0:  A * X = B     (No transpose)\n*          = 1:  A**T * X = B  (Transpose)\n*          = 2:  A**H * X = B  (Conjugate transpose)\n*\n*  N       (input) INTEGER\n*          The order of the matrix A.\n*\n*  NRHS    (input) INTEGER\n*          The number of right hand sides, i.e., the number of columns\n*          of the matrix B.  NRHS >= 0.\n*\n*  DL      (input) COMPLEX array, dimension (N-1)\n*          The (n-1) multipliers that define the matrix L from the\n*          LU factorization of A.\n*\n*  D       (input) COMPLEX array, dimension (N)\n*          The n diagonal elements of the upper triangular matrix U from\n*          the LU factorization of A.\n*\n*  DU      (input) COMPLEX array, dimension (N-1)\n*          The (n-1) elements of the first super-diagonal of U.\n*\n*  DU2     (input) COMPLEX array, dimension (N-2)\n*          The (n-2) elements of the second super-diagonal of U.\n*\n*  IPIV    (input) INTEGER array, dimension (N)\n*          The pivot indices; for 1 <= i <= n, row i of the matrix was\n*          interchanged with row IPIV(i).  IPIV(i) will always be either\n*          i or i+1; IPIV(i) = i indicates a row interchange was not\n*          required.\n*\n*  B       (input/output) COMPLEX array, dimension (LDB,NRHS)\n*          On entry, the matrix of right hand side vectors B.\n*          On exit, B is overwritten by the solution vectors X.\n*\n*  LDB     (input) INTEGER\n*          The leading dimension of the array B.  LDB >= max(1,N).\n*\n\n*  =====================================================================\n*\n*     .. Local Scalars ..\n      INTEGER            I, J\n      COMPLEX            TEMP\n*     ..\n*     .. Intrinsic Functions ..\n      INTRINSIC          CONJG\n*     ..\n\n");
      rb_exit(0);
    }
    if (rb_hash_aref(rb_options, sUsage) == Qtrue) {
      printf("%s\n", "USAGE:\n  b = NumRu::Lapack.cgtts2( itrans, dl, d, du, du2, ipiv, b, [:usage => usage, :help => help])\n");
      rb_exit(0);
    } 
  } else
    rb_options = Qnil;
  if (argc != 7)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 7)", argc);
  rblapack_itrans = argv[0];
  rblapack_dl = argv[1];
  rblapack_d = argv[2];
  rblapack_du = argv[3];
  rblapack_du2 = argv[4];
  rblapack_ipiv = argv[5];
  rblapack_b = argv[6];
  if (rb_options != Qnil) {
  }

  if (!NA_IsNArray(rblapack_ipiv))
    rb_raise(rb_eArgError, "ipiv (6th argument) must be NArray");
  if (NA_RANK(rblapack_ipiv) != 1)
    rb_raise(rb_eArgError, "rank of ipiv (6th argument) must be %d", 1);
  n = NA_SHAPE0(rblapack_ipiv);
  if (NA_TYPE(rblapack_ipiv) != NA_LINT)
    rblapack_ipiv = na_change_type(rblapack_ipiv, NA_LINT);
  ipiv = NA_PTR_TYPE(rblapack_ipiv, integer*);
  if (!NA_IsNArray(rblapack_d))
    rb_raise(rb_eArgError, "d (3th argument) must be NArray");
  if (NA_RANK(rblapack_d) != 1)
    rb_raise(rb_eArgError, "rank of d (3th argument) must be %d", 1);
  if (NA_SHAPE0(rblapack_d) != n)
    rb_raise(rb_eRuntimeError, "shape 0 of d must be the same as shape 0 of ipiv");
  if (NA_TYPE(rblapack_d) != NA_SCOMPLEX)
    rblapack_d = na_change_type(rblapack_d, NA_SCOMPLEX);
  d = NA_PTR_TYPE(rblapack_d, complex*);
  if (!NA_IsNArray(rblapack_b))
    rb_raise(rb_eArgError, "b (7th argument) must be NArray");
  if (NA_RANK(rblapack_b) != 2)
    rb_raise(rb_eArgError, "rank of b (7th argument) must be %d", 2);
  nrhs = NA_SHAPE1(rblapack_b);
  ldb = NA_SHAPE0(rblapack_b);
  if (NA_TYPE(rblapack_b) != NA_SCOMPLEX)
    rblapack_b = na_change_type(rblapack_b, NA_SCOMPLEX);
  b = NA_PTR_TYPE(rblapack_b, complex*);
  itrans = NUM2INT(rblapack_itrans);
  if (!NA_IsNArray(rblapack_du2))
    rb_raise(rb_eArgError, "du2 (5th argument) must be NArray");
  if (NA_RANK(rblapack_du2) != 1)
    rb_raise(rb_eArgError, "rank of du2 (5th argument) must be %d", 1);
  if (NA_SHAPE0(rblapack_du2) != (n-2))
    rb_raise(rb_eRuntimeError, "shape 0 of du2 must be %d", n-2);
  if (NA_TYPE(rblapack_du2) != NA_SCOMPLEX)
    rblapack_du2 = na_change_type(rblapack_du2, NA_SCOMPLEX);
  du2 = NA_PTR_TYPE(rblapack_du2, complex*);
  if (!NA_IsNArray(rblapack_du))
    rb_raise(rb_eArgError, "du (4th argument) must be NArray");
  if (NA_RANK(rblapack_du) != 1)
    rb_raise(rb_eArgError, "rank of du (4th argument) must be %d", 1);
  if (NA_SHAPE0(rblapack_du) != (n-1))
    rb_raise(rb_eRuntimeError, "shape 0 of du must be %d", n-1);
  if (NA_TYPE(rblapack_du) != NA_SCOMPLEX)
    rblapack_du = na_change_type(rblapack_du, NA_SCOMPLEX);
  du = NA_PTR_TYPE(rblapack_du, complex*);
  if (!NA_IsNArray(rblapack_dl))
    rb_raise(rb_eArgError, "dl (2th argument) must be NArray");
  if (NA_RANK(rblapack_dl) != 1)
    rb_raise(rb_eArgError, "rank of dl (2th argument) must be %d", 1);
  if (NA_SHAPE0(rblapack_dl) != (n-1))
    rb_raise(rb_eRuntimeError, "shape 0 of dl must be %d", n-1);
  if (NA_TYPE(rblapack_dl) != NA_SCOMPLEX)
    rblapack_dl = na_change_type(rblapack_dl, NA_SCOMPLEX);
  dl = NA_PTR_TYPE(rblapack_dl, complex*);
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

  cgtts2_(&itrans, &n, &nrhs, dl, d, du, du2, ipiv, b, &ldb);

  return rblapack_b;
}

void
init_lapack_cgtts2(VALUE mLapack){
  rb_define_module_function(mLapack, "cgtts2", rblapack_cgtts2, -1);
  sHelp = ID2SYM(rb_intern("help"));
  sUsage = ID2SYM(rb_intern("usage"));
}
