#include "rb_lapack.h"

extern VOID zptsv_(integer *n, integer *nrhs, doublereal *d, doublecomplex *e, doublecomplex *b, integer *ldb, integer *info);

static VALUE sHelp, sUsage;

static VALUE
rblapack_zptsv(int argc, VALUE *argv, VALUE self){
  VALUE rblapack_nrhs;
  integer nrhs; 
  VALUE rblapack_d;
  doublereal *d; 
  VALUE rblapack_e;
  doublecomplex *e; 
  VALUE rblapack_b;
  doublecomplex *b; 
  VALUE rblapack_info;
  integer info; 
  VALUE rblapack_d_out__;
  doublereal *d_out__;
  VALUE rblapack_e_out__;
  doublecomplex *e_out__;
  VALUE rblapack_b_out__;
  doublecomplex *b_out__;

  integer n;
  integer ldb;

  VALUE rb_options;
  if (argc > 0 && TYPE(argv[argc-1]) == T_HASH) {
    argc--;
    rb_options = argv[argc];
    if (rb_hash_aref(rb_options, sHelp) == Qtrue) {
      printf("%s\n", "USAGE:\n  info, d, e, b = NumRu::Lapack.zptsv( nrhs, d, e, b, [:usage => usage, :help => help])\n\n\nFORTRAN MANUAL\n      SUBROUTINE ZPTSV( N, NRHS, D, E, B, LDB, INFO )\n\n*  Purpose\n*  =======\n*\n*  ZPTSV computes the solution to a complex system of linear equations\n*  A*X = B, where A is an N-by-N Hermitian positive definite tridiagonal\n*  matrix, and X and B are N-by-NRHS matrices.\n*\n*  A is factored as A = L*D*L**H, and the factored form of A is then\n*  used to solve the system of equations.\n*\n\n*  Arguments\n*  =========\n*\n*  N       (input) INTEGER\n*          The order of the matrix A.  N >= 0.\n*\n*  NRHS    (input) INTEGER\n*          The number of right hand sides, i.e., the number of columns\n*          of the matrix B.  NRHS >= 0.\n*\n*  D       (input/output) DOUBLE PRECISION array, dimension (N)\n*          On entry, the n diagonal elements of the tridiagonal matrix\n*          A.  On exit, the n diagonal elements of the diagonal matrix\n*          D from the factorization A = L*D*L**H.\n*\n*  E       (input/output) COMPLEX*16 array, dimension (N-1)\n*          On entry, the (n-1) subdiagonal elements of the tridiagonal\n*          matrix A.  On exit, the (n-1) subdiagonal elements of the\n*          unit bidiagonal factor L from the L*D*L**H factorization of\n*          A.  E can also be regarded as the superdiagonal of the unit\n*          bidiagonal factor U from the U**H*D*U factorization of A.\n*\n*  B       (input/output) COMPLEX*16 array, dimension (LDB,N)\n*          On entry, the N-by-NRHS right hand side matrix B.\n*          On exit, if INFO = 0, the N-by-NRHS solution matrix X.\n*\n*  LDB     (input) INTEGER\n*          The leading dimension of the array B.  LDB >= max(1,N).\n*\n*  INFO    (output) INTEGER\n*          = 0:  successful exit\n*          < 0:  if INFO = -i, the i-th argument had an illegal value\n*          > 0:  if INFO = i, the leading minor of order i is not\n*                positive definite, and the solution has not been\n*                computed.  The factorization has not been completed\n*                unless i = N.\n*\n\n*  =====================================================================\n*\n*     .. External Subroutines ..\n      EXTERNAL           XERBLA, ZPTTRF, ZPTTRS\n*     ..\n*     .. Intrinsic Functions ..\n      INTRINSIC          MAX\n*     ..\n\n");
      rb_exit(0);
    }
    if (rb_hash_aref(rb_options, sUsage) == Qtrue) {
      printf("%s\n", "USAGE:\n  info, d, e, b = NumRu::Lapack.zptsv( nrhs, d, e, b, [:usage => usage, :help => help])\n");
      rb_exit(0);
    } 
  } else
    rb_options = Qnil;
  if (argc != 4)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 4)", argc);
  rblapack_nrhs = argv[0];
  rblapack_d = argv[1];
  rblapack_e = argv[2];
  rblapack_b = argv[3];
  if (rb_options != Qnil) {
  }

  if (!NA_IsNArray(rblapack_b))
    rb_raise(rb_eArgError, "b (4th argument) must be NArray");
  if (NA_RANK(rblapack_b) != 2)
    rb_raise(rb_eArgError, "rank of b (4th argument) must be %d", 2);
  n = NA_SHAPE1(rblapack_b);
  ldb = NA_SHAPE0(rblapack_b);
  if (NA_TYPE(rblapack_b) != NA_DCOMPLEX)
    rblapack_b = na_change_type(rblapack_b, NA_DCOMPLEX);
  b = NA_PTR_TYPE(rblapack_b, doublecomplex*);
  if (!NA_IsNArray(rblapack_d))
    rb_raise(rb_eArgError, "d (2th argument) must be NArray");
  if (NA_RANK(rblapack_d) != 1)
    rb_raise(rb_eArgError, "rank of d (2th argument) must be %d", 1);
  if (NA_SHAPE0(rblapack_d) != n)
    rb_raise(rb_eRuntimeError, "shape 0 of d must be the same as shape 1 of b");
  if (NA_TYPE(rblapack_d) != NA_DFLOAT)
    rblapack_d = na_change_type(rblapack_d, NA_DFLOAT);
  d = NA_PTR_TYPE(rblapack_d, doublereal*);
  nrhs = NUM2INT(rblapack_nrhs);
  if (!NA_IsNArray(rblapack_e))
    rb_raise(rb_eArgError, "e (3th argument) must be NArray");
  if (NA_RANK(rblapack_e) != 1)
    rb_raise(rb_eArgError, "rank of e (3th argument) must be %d", 1);
  if (NA_SHAPE0(rblapack_e) != (n-1))
    rb_raise(rb_eRuntimeError, "shape 0 of e must be %d", n-1);
  if (NA_TYPE(rblapack_e) != NA_DCOMPLEX)
    rblapack_e = na_change_type(rblapack_e, NA_DCOMPLEX);
  e = NA_PTR_TYPE(rblapack_e, doublecomplex*);
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
    rblapack_e_out__ = na_make_object(NA_DCOMPLEX, 1, shape, cNArray);
  }
  e_out__ = NA_PTR_TYPE(rblapack_e_out__, doublecomplex*);
  MEMCPY(e_out__, e, doublecomplex, NA_TOTAL(rblapack_e));
  rblapack_e = rblapack_e_out__;
  e = e_out__;
  {
    int shape[2];
    shape[0] = ldb;
    shape[1] = n;
    rblapack_b_out__ = na_make_object(NA_DCOMPLEX, 2, shape, cNArray);
  }
  b_out__ = NA_PTR_TYPE(rblapack_b_out__, doublecomplex*);
  MEMCPY(b_out__, b, doublecomplex, NA_TOTAL(rblapack_b));
  rblapack_b = rblapack_b_out__;
  b = b_out__;

  zptsv_(&n, &nrhs, d, e, b, &ldb, &info);

  rblapack_info = INT2NUM(info);
  return rb_ary_new3(4, rblapack_info, rblapack_d, rblapack_e, rblapack_b);
}

void
init_lapack_zptsv(VALUE mLapack){
  rb_define_module_function(mLapack, "zptsv", rblapack_zptsv, -1);
  sHelp = ID2SYM(rb_intern("help"));
  sUsage = ID2SYM(rb_intern("usage"));
}
