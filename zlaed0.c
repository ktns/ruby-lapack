#include "rb_lapack.h"

extern VOID zlaed0_(integer *qsiz, integer *n, doublereal *d, doublereal *e, doublecomplex *q, integer *ldq, doublecomplex *qstore, integer *ldqs, doublereal *rwork, integer *iwork, integer *info);

static VALUE sHelp, sUsage;

static VALUE
rblapack_zlaed0(int argc, VALUE *argv, VALUE self){
  VALUE rblapack_qsiz;
  integer qsiz; 
  VALUE rblapack_d;
  doublereal *d; 
  VALUE rblapack_e;
  doublereal *e; 
  VALUE rblapack_q;
  doublecomplex *q; 
  VALUE rblapack_info;
  integer info; 
  VALUE rblapack_d_out__;
  doublereal *d_out__;
  VALUE rblapack_e_out__;
  doublereal *e_out__;
  VALUE rblapack_q_out__;
  doublecomplex *q_out__;
  doublecomplex *qstore;
  doublereal *rwork;
  integer *iwork;

  integer n;
  integer ldq;
  integer ldqs;

  VALUE rb_options;
  if (argc > 0 && TYPE(argv[argc-1]) == T_HASH) {
    argc--;
    rb_options = argv[argc];
    if (rb_hash_aref(rb_options, sHelp) == Qtrue) {
      printf("%s\n", "USAGE:\n  info, d, e, q = NumRu::Lapack.zlaed0( qsiz, d, e, q, [:usage => usage, :help => help])\n\n\nFORTRAN MANUAL\n      SUBROUTINE ZLAED0( QSIZ, N, D, E, Q, LDQ, QSTORE, LDQS, RWORK, IWORK, INFO )\n\n*  Purpose\n*  =======\n*\n*  Using the divide and conquer method, ZLAED0 computes all eigenvalues\n*  of a symmetric tridiagonal matrix which is one diagonal block of\n*  those from reducing a dense or band Hermitian matrix and\n*  corresponding eigenvectors of the dense or band matrix.\n*\n\n*  Arguments\n*  =========\n*\n*  QSIZ   (input) INTEGER\n*         The dimension of the unitary matrix used to reduce\n*         the full matrix to tridiagonal form.  QSIZ >= N if ICOMPQ = 1.\n*\n*  N      (input) INTEGER\n*         The dimension of the symmetric tridiagonal matrix.  N >= 0.\n*\n*  D      (input/output) DOUBLE PRECISION array, dimension (N)\n*         On entry, the diagonal elements of the tridiagonal matrix.\n*         On exit, the eigenvalues in ascending order.\n*\n*  E      (input/output) DOUBLE PRECISION array, dimension (N-1)\n*         On entry, the off-diagonal elements of the tridiagonal matrix.\n*         On exit, E has been destroyed.\n*\n*  Q      (input/output) COMPLEX*16 array, dimension (LDQ,N)\n*         On entry, Q must contain an QSIZ x N matrix whose columns\n*         unitarily orthonormal. It is a part of the unitary matrix\n*         that reduces the full dense Hermitian matrix to a\n*         (reducible) symmetric tridiagonal matrix.\n*\n*  LDQ    (input) INTEGER\n*         The leading dimension of the array Q.  LDQ >= max(1,N).\n*\n*  IWORK  (workspace) INTEGER array,\n*         the dimension of IWORK must be at least\n*                      6 + 6*N + 5*N*lg N\n*                      ( lg( N ) = smallest integer k\n*                                  such that 2^k >= N )\n*\n*  RWORK  (workspace) DOUBLE PRECISION array,\n*                               dimension (1 + 3*N + 2*N*lg N + 3*N**2)\n*                        ( lg( N ) = smallest integer k\n*                                    such that 2^k >= N )\n*\n*  QSTORE (workspace) COMPLEX*16 array, dimension (LDQS, N)\n*         Used to store parts of\n*         the eigenvector matrix when the updating matrix multiplies\n*         take place.\n*\n*  LDQS   (input) INTEGER\n*         The leading dimension of the array QSTORE.\n*         LDQS >= max(1,N).\n*\n*  INFO   (output) INTEGER\n*          = 0:  successful exit.\n*          < 0:  if INFO = -i, the i-th argument had an illegal value.\n*          > 0:  The algorithm failed to compute an eigenvalue while\n*                working on the submatrix lying in rows and columns\n*                INFO/(N+1) through mod(INFO,N+1).\n*\n\n*  =====================================================================\n*\n*  Warning:      N could be as big as QSIZ!\n*\n\n");
      rb_exit(0);
    }
    if (rb_hash_aref(rb_options, sUsage) == Qtrue) {
      printf("%s\n", "USAGE:\n  info, d, e, q = NumRu::Lapack.zlaed0( qsiz, d, e, q, [:usage => usage, :help => help])\n");
      rb_exit(0);
    } 
  } else
    rb_options = Qnil;
  if (argc != 4)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 4)", argc);
  rblapack_qsiz = argv[0];
  rblapack_d = argv[1];
  rblapack_e = argv[2];
  rblapack_q = argv[3];
  if (rb_options != Qnil) {
  }

  qsiz = NUM2INT(rblapack_qsiz);
  if (!NA_IsNArray(rblapack_d))
    rb_raise(rb_eArgError, "d (2th argument) must be NArray");
  if (NA_RANK(rblapack_d) != 1)
    rb_raise(rb_eArgError, "rank of d (2th argument) must be %d", 1);
  n = NA_SHAPE0(rblapack_d);
  if (NA_TYPE(rblapack_d) != NA_DFLOAT)
    rblapack_d = na_change_type(rblapack_d, NA_DFLOAT);
  d = NA_PTR_TYPE(rblapack_d, doublereal*);
  if (!NA_IsNArray(rblapack_q))
    rb_raise(rb_eArgError, "q (4th argument) must be NArray");
  if (NA_RANK(rblapack_q) != 2)
    rb_raise(rb_eArgError, "rank of q (4th argument) must be %d", 2);
  if (NA_SHAPE1(rblapack_q) != n)
    rb_raise(rb_eRuntimeError, "shape 1 of q must be the same as shape 0 of d");
  ldq = NA_SHAPE0(rblapack_q);
  if (NA_TYPE(rblapack_q) != NA_DCOMPLEX)
    rblapack_q = na_change_type(rblapack_q, NA_DCOMPLEX);
  q = NA_PTR_TYPE(rblapack_q, doublecomplex*);
  if (!NA_IsNArray(rblapack_e))
    rb_raise(rb_eArgError, "e (3th argument) must be NArray");
  if (NA_RANK(rblapack_e) != 1)
    rb_raise(rb_eArgError, "rank of e (3th argument) must be %d", 1);
  if (NA_SHAPE0(rblapack_e) != (n-1))
    rb_raise(rb_eRuntimeError, "shape 0 of e must be %d", n-1);
  if (NA_TYPE(rblapack_e) != NA_DFLOAT)
    rblapack_e = na_change_type(rblapack_e, NA_DFLOAT);
  e = NA_PTR_TYPE(rblapack_e, doublereal*);
  ldqs = MAX(1,n);
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
  {
    int shape[2];
    shape[0] = ldq;
    shape[1] = n;
    rblapack_q_out__ = na_make_object(NA_DCOMPLEX, 2, shape, cNArray);
  }
  q_out__ = NA_PTR_TYPE(rblapack_q_out__, doublecomplex*);
  MEMCPY(q_out__, q, doublecomplex, NA_TOTAL(rblapack_q));
  rblapack_q = rblapack_q_out__;
  q = q_out__;
  qstore = ALLOC_N(doublecomplex, (ldqs)*(n));
  rwork = ALLOC_N(doublereal, (1 + 3*n + 2*n*LG(n) + 3*pow(n,2)));
  iwork = ALLOC_N(integer, (6 + 6*n + 5*n*LG(n)));

  zlaed0_(&qsiz, &n, d, e, q, &ldq, qstore, &ldqs, rwork, iwork, &info);

  free(qstore);
  free(rwork);
  free(iwork);
  rblapack_info = INT2NUM(info);
  return rb_ary_new3(4, rblapack_info, rblapack_d, rblapack_e, rblapack_q);
}

void
init_lapack_zlaed0(VALUE mLapack){
  rb_define_module_function(mLapack, "zlaed0", rblapack_zlaed0, -1);
  sHelp = ID2SYM(rb_intern("help"));
  sUsage = ID2SYM(rb_intern("usage"));
}
