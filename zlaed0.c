#include "rb_lapack.h"

static VALUE
rb_zlaed0(int argc, VALUE *argv, VALUE self){
  VALUE rb_qsiz;
  integer qsiz; 
  VALUE rb_d;
  doublereal *d; 
  VALUE rb_e;
  doublereal *e; 
  VALUE rb_q;
  doublecomplex *q; 
  VALUE rb_info;
  integer info; 
  VALUE rb_d_out__;
  doublereal *d_out__;
  VALUE rb_e_out__;
  doublereal *e_out__;
  VALUE rb_q_out__;
  doublecomplex *q_out__;
  doublecomplex *qstore;
  doublereal *rwork;
  integer *iwork;

  integer n;
  integer ldq;
  integer ldqs;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  info, d, e, q = NumRu::Lapack.zlaed0( qsiz, d, e, q)\n    or\n  NumRu::Lapack.zlaed0  # print help\n\n\nFORTRAN MANUAL\n      SUBROUTINE ZLAED0( QSIZ, N, D, E, Q, LDQ, QSTORE, LDQS, RWORK, IWORK, INFO )\n\n*  Purpose\n*  =======\n*\n*  Using the divide and conquer method, ZLAED0 computes all eigenvalues\n*  of a symmetric tridiagonal matrix which is one diagonal block of\n*  those from reducing a dense or band Hermitian matrix and\n*  corresponding eigenvectors of the dense or band matrix.\n*\n\n*  Arguments\n*  =========\n*\n*  QSIZ   (input) INTEGER\n*         The dimension of the unitary matrix used to reduce\n*         the full matrix to tridiagonal form.  QSIZ >= N if ICOMPQ = 1.\n*\n*  N      (input) INTEGER\n*         The dimension of the symmetric tridiagonal matrix.  N >= 0.\n*\n*  D      (input/output) DOUBLE PRECISION array, dimension (N)\n*         On entry, the diagonal elements of the tridiagonal matrix.\n*         On exit, the eigenvalues in ascending order.\n*\n*  E      (input/output) DOUBLE PRECISION array, dimension (N-1)\n*         On entry, the off-diagonal elements of the tridiagonal matrix.\n*         On exit, E has been destroyed.\n*\n*  Q      (input/output) COMPLEX*16 array, dimension (LDQ,N)\n*         On entry, Q must contain an QSIZ x N matrix whose columns\n*         unitarily orthonormal. It is a part of the unitary matrix\n*         that reduces the full dense Hermitian matrix to a\n*         (reducible) symmetric tridiagonal matrix.\n*\n*  LDQ    (input) INTEGER\n*         The leading dimension of the array Q.  LDQ >= max(1,N).\n*\n*  IWORK  (workspace) INTEGER array,\n*         the dimension of IWORK must be at least\n*                      6 + 6*N + 5*N*lg N\n*                      ( lg( N ) = smallest integer k\n*                                  such that 2^k >= N )\n*\n*  RWORK  (workspace) DOUBLE PRECISION array,\n*                               dimension (1 + 3*N + 2*N*lg N + 3*N**2)\n*                        ( lg( N ) = smallest integer k\n*                                    such that 2^k >= N )\n*\n*  QSTORE (workspace) COMPLEX*16 array, dimension (LDQS, N)\n*         Used to store parts of\n*         the eigenvector matrix when the updating matrix multiplies\n*         take place.\n*\n*  LDQS   (input) INTEGER\n*         The leading dimension of the array QSTORE.\n*         LDQS >= max(1,N).\n*\n*  INFO   (output) INTEGER\n*          = 0:  successful exit.\n*          < 0:  if INFO = -i, the i-th argument had an illegal value.\n*          > 0:  The algorithm failed to compute an eigenvalue while\n*                working on the submatrix lying in rows and columns\n*                INFO/(N+1) through mod(INFO,N+1).\n*\n\n*  =====================================================================\n*\n*  Warning:      N could be as big as QSIZ!\n*\n\n");
    return Qnil;
  }
  if (argc != 4)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 4)", argc);
  rb_qsiz = argv[0];
  rb_d = argv[1];
  rb_e = argv[2];
  rb_q = argv[3];

  qsiz = NUM2INT(rb_qsiz);
  if (!NA_IsNArray(rb_d))
    rb_raise(rb_eArgError, "d (2th argument) must be NArray");
  if (NA_RANK(rb_d) != 1)
    rb_raise(rb_eArgError, "rank of d (2th argument) must be %d", 1);
  n = NA_SHAPE0(rb_d);
  if (NA_TYPE(rb_d) != NA_DFLOAT)
    rb_d = na_change_type(rb_d, NA_DFLOAT);
  d = NA_PTR_TYPE(rb_d, doublereal*);
  if (!NA_IsNArray(rb_e))
    rb_raise(rb_eArgError, "e (3th argument) must be NArray");
  if (NA_RANK(rb_e) != 1)
    rb_raise(rb_eArgError, "rank of e (3th argument) must be %d", 1);
  if (NA_SHAPE0(rb_e) != (n-1))
    rb_raise(rb_eRuntimeError, "shape 0 of e must be %d", n-1);
  if (NA_TYPE(rb_e) != NA_DFLOAT)
    rb_e = na_change_type(rb_e, NA_DFLOAT);
  e = NA_PTR_TYPE(rb_e, doublereal*);
  if (!NA_IsNArray(rb_q))
    rb_raise(rb_eArgError, "q (4th argument) must be NArray");
  if (NA_RANK(rb_q) != 2)
    rb_raise(rb_eArgError, "rank of q (4th argument) must be %d", 2);
  ldq = NA_SHAPE0(rb_q);
  if (NA_SHAPE1(rb_q) != n)
    rb_raise(rb_eRuntimeError, "shape 1 of q must be the same as shape 0 of d");
  if (NA_TYPE(rb_q) != NA_DCOMPLEX)
    rb_q = na_change_type(rb_q, NA_DCOMPLEX);
  q = NA_PTR_TYPE(rb_q, doublecomplex*);
  {
    int shape[1];
    shape[0] = n;
    rb_d_out__ = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  d_out__ = NA_PTR_TYPE(rb_d_out__, doublereal*);
  MEMCPY(d_out__, d, doublereal, NA_TOTAL(rb_d));
  rb_d = rb_d_out__;
  d = d_out__;
  {
    int shape[1];
    shape[0] = n-1;
    rb_e_out__ = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  e_out__ = NA_PTR_TYPE(rb_e_out__, doublereal*);
  MEMCPY(e_out__, e, doublereal, NA_TOTAL(rb_e));
  rb_e = rb_e_out__;
  e = e_out__;
  {
    int shape[2];
    shape[0] = ldq;
    shape[1] = n;
    rb_q_out__ = na_make_object(NA_DCOMPLEX, 2, shape, cNArray);
  }
  q_out__ = NA_PTR_TYPE(rb_q_out__, doublecomplex*);
  MEMCPY(q_out__, q, doublecomplex, NA_TOTAL(rb_q));
  rb_q = rb_q_out__;
  q = q_out__;
  ldqs = MAX(1,n);
  qstore = ALLOC_N(doublecomplex, (ldqs)*(n));
  rwork = ALLOC_N(doublereal, (1 + 3*n + 2*n*LG(n) + 3*pow(n,2)));
  iwork = ALLOC_N(integer, (6 + 6*n + 5*n*LG(n)));

  zlaed0_(&qsiz, &n, d, e, q, &ldq, qstore, &ldqs, rwork, iwork, &info);

  free(qstore);
  free(rwork);
  free(iwork);
  rb_info = INT2NUM(info);
  return rb_ary_new3(4, rb_info, rb_d, rb_e, rb_q);
}

void
init_lapack_zlaed0(VALUE mLapack){
  rb_define_module_function(mLapack, "zlaed0", rb_zlaed0, -1);
}
