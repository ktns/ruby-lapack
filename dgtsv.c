#include "rb_lapack.h"

extern VOID dgtsv_(integer *n, integer *nrhs, doublereal *dl, doublereal *d, doublereal *du, doublereal *b, integer *ldb, integer *info);

static VALUE sHelp, sUsage;

static VALUE
rblapack_dgtsv(int argc, VALUE *argv, VALUE self){
  VALUE rblapack_dl;
  doublereal *dl; 
  VALUE rblapack_d;
  doublereal *d; 
  VALUE rblapack_du;
  doublereal *du; 
  VALUE rblapack_b;
  doublereal *b; 
  VALUE rblapack_info;
  integer info; 
  VALUE rblapack_dl_out__;
  doublereal *dl_out__;
  VALUE rblapack_d_out__;
  doublereal *d_out__;
  VALUE rblapack_du_out__;
  doublereal *du_out__;
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
      printf("%s\n", "USAGE:\n  info, dl, d, du, b = NumRu::Lapack.dgtsv( dl, d, du, b, [:usage => usage, :help => help])\n\n\nFORTRAN MANUAL\n      SUBROUTINE DGTSV( N, NRHS, DL, D, DU, B, LDB, INFO )\n\n*  Purpose\n*  =======\n*\n*  DGTSV  solves the equation\n*\n*     A*X = B,\n*\n*  where A is an n by n tridiagonal matrix, by Gaussian elimination with\n*  partial pivoting.\n*\n*  Note that the equation  A'*X = B  may be solved by interchanging the\n*  order of the arguments DU and DL.\n*\n\n*  Arguments\n*  =========\n*\n*  N       (input) INTEGER\n*          The order of the matrix A.  N >= 0.\n*\n*  NRHS    (input) INTEGER\n*          The number of right hand sides, i.e., the number of columns\n*          of the matrix B.  NRHS >= 0.\n*\n*  DL      (input/output) DOUBLE PRECISION array, dimension (N-1)\n*          On entry, DL must contain the (n-1) sub-diagonal elements of\n*          A.\n*\n*          On exit, DL is overwritten by the (n-2) elements of the\n*          second super-diagonal of the upper triangular matrix U from\n*          the LU factorization of A, in DL(1), ..., DL(n-2).\n*\n*  D       (input/output) DOUBLE PRECISION array, dimension (N)\n*          On entry, D must contain the diagonal elements of A.\n*\n*          On exit, D is overwritten by the n diagonal elements of U.\n*\n*  DU      (input/output) DOUBLE PRECISION array, dimension (N-1)\n*          On entry, DU must contain the (n-1) super-diagonal elements\n*          of A.\n*\n*          On exit, DU is overwritten by the (n-1) elements of the first\n*          super-diagonal of U.\n*\n*  B       (input/output) DOUBLE PRECISION array, dimension (LDB,NRHS)\n*          On entry, the N by NRHS matrix of right hand side matrix B.\n*          On exit, if INFO = 0, the N by NRHS solution matrix X.\n*\n*  LDB     (input) INTEGER\n*          The leading dimension of the array B.  LDB >= max(1,N).\n*\n*  INFO    (output) INTEGER\n*          = 0: successful exit\n*          < 0: if INFO = -i, the i-th argument had an illegal value\n*          > 0: if INFO = i, U(i,i) is exactly zero, and the solution\n*               has not been computed.  The factorization has not been\n*               completed unless i = N.\n*\n\n*  =====================================================================\n*\n\n");
      rb_exit(0);
    }
    if (rb_hash_aref(rb_options, sUsage) == Qtrue) {
      printf("%s\n", "USAGE:\n  info, dl, d, du, b = NumRu::Lapack.dgtsv( dl, d, du, b, [:usage => usage, :help => help])\n");
      rb_exit(0);
    } 
  } else
    rb_options = Qnil;
  if (argc != 4)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 4)", argc);
  rblapack_dl = argv[0];
  rblapack_d = argv[1];
  rblapack_du = argv[2];
  rblapack_b = argv[3];
  if (rb_options != Qnil) {
  }

  if (!NA_IsNArray(rblapack_b))
    rb_raise(rb_eArgError, "b (4th argument) must be NArray");
  if (NA_RANK(rblapack_b) != 2)
    rb_raise(rb_eArgError, "rank of b (4th argument) must be %d", 2);
  nrhs = NA_SHAPE1(rblapack_b);
  ldb = NA_SHAPE0(rblapack_b);
  if (NA_TYPE(rblapack_b) != NA_DFLOAT)
    rblapack_b = na_change_type(rblapack_b, NA_DFLOAT);
  b = NA_PTR_TYPE(rblapack_b, doublereal*);
  if (!NA_IsNArray(rblapack_d))
    rb_raise(rb_eArgError, "d (2th argument) must be NArray");
  if (NA_RANK(rblapack_d) != 1)
    rb_raise(rb_eArgError, "rank of d (2th argument) must be %d", 1);
  n = NA_SHAPE0(rblapack_d);
  if (NA_TYPE(rblapack_d) != NA_DFLOAT)
    rblapack_d = na_change_type(rblapack_d, NA_DFLOAT);
  d = NA_PTR_TYPE(rblapack_d, doublereal*);
  if (!NA_IsNArray(rblapack_du))
    rb_raise(rb_eArgError, "du (3th argument) must be NArray");
  if (NA_RANK(rblapack_du) != 1)
    rb_raise(rb_eArgError, "rank of du (3th argument) must be %d", 1);
  if (NA_SHAPE0(rblapack_du) != (n-1))
    rb_raise(rb_eRuntimeError, "shape 0 of du must be %d", n-1);
  if (NA_TYPE(rblapack_du) != NA_DFLOAT)
    rblapack_du = na_change_type(rblapack_du, NA_DFLOAT);
  du = NA_PTR_TYPE(rblapack_du, doublereal*);
  if (!NA_IsNArray(rblapack_dl))
    rb_raise(rb_eArgError, "dl (1th argument) must be NArray");
  if (NA_RANK(rblapack_dl) != 1)
    rb_raise(rb_eArgError, "rank of dl (1th argument) must be %d", 1);
  if (NA_SHAPE0(rblapack_dl) != (n-1))
    rb_raise(rb_eRuntimeError, "shape 0 of dl must be %d", n-1);
  if (NA_TYPE(rblapack_dl) != NA_DFLOAT)
    rblapack_dl = na_change_type(rblapack_dl, NA_DFLOAT);
  dl = NA_PTR_TYPE(rblapack_dl, doublereal*);
  {
    int shape[1];
    shape[0] = n-1;
    rblapack_dl_out__ = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  dl_out__ = NA_PTR_TYPE(rblapack_dl_out__, doublereal*);
  MEMCPY(dl_out__, dl, doublereal, NA_TOTAL(rblapack_dl));
  rblapack_dl = rblapack_dl_out__;
  dl = dl_out__;
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
    rblapack_du_out__ = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  du_out__ = NA_PTR_TYPE(rblapack_du_out__, doublereal*);
  MEMCPY(du_out__, du, doublereal, NA_TOTAL(rblapack_du));
  rblapack_du = rblapack_du_out__;
  du = du_out__;
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

  dgtsv_(&n, &nrhs, dl, d, du, b, &ldb, &info);

  rblapack_info = INT2NUM(info);
  return rb_ary_new3(5, rblapack_info, rblapack_dl, rblapack_d, rblapack_du, rblapack_b);
}

void
init_lapack_dgtsv(VALUE mLapack){
  rb_define_module_function(mLapack, "dgtsv", rblapack_dgtsv, -1);
  sHelp = ID2SYM(rb_intern("help"));
  sUsage = ID2SYM(rb_intern("usage"));
}
