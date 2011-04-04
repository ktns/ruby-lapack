#include "rb_lapack.h"

extern VOID dlaed0_(integer *icompq, integer *qsiz, integer *n, doublereal *d, doublereal *e, doublereal *q, integer *ldq, doublereal *qstore, integer *ldqs, doublereal *work, integer *iwork, integer *info);

static VALUE
rb_dlaed0(int argc, VALUE *argv, VALUE self){
  VALUE rb_icompq;
  integer icompq; 
  VALUE rb_qsiz;
  integer qsiz; 
  VALUE rb_d;
  doublereal *d; 
  VALUE rb_e;
  doublereal *e; 
  VALUE rb_q;
  doublereal *q; 
  VALUE rb_info;
  integer info; 
  VALUE rb_d_out__;
  doublereal *d_out__;
  VALUE rb_q_out__;
  doublereal *q_out__;
  doublereal *qstore;
  doublereal *work;
  integer *iwork;

  integer n;
  integer ldq;
  integer ldqs;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  info, d, q = NumRu::Lapack.dlaed0( icompq, qsiz, d, e, q)\n    or\n  NumRu::Lapack.dlaed0  # print help\n\n\nFORTRAN MANUAL\n      SUBROUTINE DLAED0( ICOMPQ, QSIZ, N, D, E, Q, LDQ, QSTORE, LDQS, WORK, IWORK, INFO )\n\n*  Purpose\n*  =======\n*\n*  DLAED0 computes all eigenvalues and corresponding eigenvectors of a\n*  symmetric tridiagonal matrix using the divide and conquer method.\n*\n\n*  Arguments\n*  =========\n*\n*  ICOMPQ  (input) INTEGER\n*          = 0:  Compute eigenvalues only.\n*          = 1:  Compute eigenvectors of original dense symmetric matrix\n*                also.  On entry, Q contains the orthogonal matrix used\n*                to reduce the original matrix to tridiagonal form.\n*          = 2:  Compute eigenvalues and eigenvectors of tridiagonal\n*                matrix.\n*\n*  QSIZ   (input) INTEGER\n*         The dimension of the orthogonal matrix used to reduce\n*         the full matrix to tridiagonal form.  QSIZ >= N if ICOMPQ = 1.\n*\n*  N      (input) INTEGER\n*         The dimension of the symmetric tridiagonal matrix.  N >= 0.\n*\n*  D      (input/output) DOUBLE PRECISION array, dimension (N)\n*         On entry, the main diagonal of the tridiagonal matrix.\n*         On exit, its eigenvalues.\n*\n*  E      (input) DOUBLE PRECISION array, dimension (N-1)\n*         The off-diagonal elements of the tridiagonal matrix.\n*         On exit, E has been destroyed.\n*\n*  Q      (input/output) DOUBLE PRECISION array, dimension (LDQ, N)\n*         On entry, Q must contain an N-by-N orthogonal matrix.\n*         If ICOMPQ = 0    Q is not referenced.\n*         If ICOMPQ = 1    On entry, Q is a subset of the columns of the\n*                          orthogonal matrix used to reduce the full\n*                          matrix to tridiagonal form corresponding to\n*                          the subset of the full matrix which is being\n*                          decomposed at this time.\n*         If ICOMPQ = 2    On entry, Q will be the identity matrix.\n*                          On exit, Q contains the eigenvectors of the\n*                          tridiagonal matrix.\n*\n*  LDQ    (input) INTEGER\n*         The leading dimension of the array Q.  If eigenvectors are\n*         desired, then  LDQ >= max(1,N).  In any case,  LDQ >= 1.\n*\n*  QSTORE (workspace) DOUBLE PRECISION array, dimension (LDQS, N)\n*         Referenced only when ICOMPQ = 1.  Used to store parts of\n*         the eigenvector matrix when the updating matrix multiplies\n*         take place.\n*\n*  LDQS   (input) INTEGER\n*         The leading dimension of the array QSTORE.  If ICOMPQ = 1,\n*         then  LDQS >= max(1,N).  In any case,  LDQS >= 1.\n*\n*  WORK   (workspace) DOUBLE PRECISION array,\n*         If ICOMPQ = 0 or 1, the dimension of WORK must be at least\n*                     1 + 3*N + 2*N*lg N + 2*N**2\n*                     ( lg( N ) = smallest integer k\n*                                 such that 2^k >= N )\n*         If ICOMPQ = 2, the dimension of WORK must be at least\n*                     4*N + N**2.\n*\n*  IWORK  (workspace) INTEGER array,\n*         If ICOMPQ = 0 or 1, the dimension of IWORK must be at least\n*                        6 + 6*N + 5*N*lg N.\n*                        ( lg( N ) = smallest integer k\n*                                    such that 2^k >= N )\n*         If ICOMPQ = 2, the dimension of IWORK must be at least\n*                        3 + 5*N.\n*\n*  INFO   (output) INTEGER\n*          = 0:  successful exit.\n*          < 0:  if INFO = -i, the i-th argument had an illegal value.\n*          > 0:  The algorithm failed to compute an eigenvalue while\n*                working on the submatrix lying in rows and columns\n*                INFO/(N+1) through mod(INFO,N+1).\n*\n\n*  Further Details\n*  ===============\n*\n*  Based on contributions by\n*     Jeff Rutter, Computer Science Division, University of California\n*     at Berkeley, USA\n*\n*  =====================================================================\n*\n\n");
    return Qnil;
  }
  if (argc != 5)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 5)", argc);
  rb_icompq = argv[0];
  rb_qsiz = argv[1];
  rb_d = argv[2];
  rb_e = argv[3];
  rb_q = argv[4];

  qsiz = NUM2INT(rb_qsiz);
  if (!NA_IsNArray(rb_d))
    rb_raise(rb_eArgError, "d (3th argument) must be NArray");
  if (NA_RANK(rb_d) != 1)
    rb_raise(rb_eArgError, "rank of d (3th argument) must be %d", 1);
  n = NA_SHAPE0(rb_d);
  if (NA_TYPE(rb_d) != NA_DFLOAT)
    rb_d = na_change_type(rb_d, NA_DFLOAT);
  d = NA_PTR_TYPE(rb_d, doublereal*);
  if (!NA_IsNArray(rb_q))
    rb_raise(rb_eArgError, "q (5th argument) must be NArray");
  if (NA_RANK(rb_q) != 2)
    rb_raise(rb_eArgError, "rank of q (5th argument) must be %d", 2);
  if (NA_SHAPE1(rb_q) != n)
    rb_raise(rb_eRuntimeError, "shape 1 of q must be the same as shape 0 of d");
  ldq = NA_SHAPE0(rb_q);
  if (NA_TYPE(rb_q) != NA_DFLOAT)
    rb_q = na_change_type(rb_q, NA_DFLOAT);
  q = NA_PTR_TYPE(rb_q, doublereal*);
  icompq = NUM2INT(rb_icompq);
  if (!NA_IsNArray(rb_e))
    rb_raise(rb_eArgError, "e (4th argument) must be NArray");
  if (NA_RANK(rb_e) != 1)
    rb_raise(rb_eArgError, "rank of e (4th argument) must be %d", 1);
  if (NA_SHAPE0(rb_e) != (n-1))
    rb_raise(rb_eRuntimeError, "shape 0 of e must be %d", n-1);
  if (NA_TYPE(rb_e) != NA_DFLOAT)
    rb_e = na_change_type(rb_e, NA_DFLOAT);
  e = NA_PTR_TYPE(rb_e, doublereal*);
  ldqs = icompq == 1 ? MAX(1,n) : 1;
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
    int shape[2];
    shape[0] = ldq;
    shape[1] = n;
    rb_q_out__ = na_make_object(NA_DFLOAT, 2, shape, cNArray);
  }
  q_out__ = NA_PTR_TYPE(rb_q_out__, doublereal*);
  MEMCPY(q_out__, q, doublereal, NA_TOTAL(rb_q));
  rb_q = rb_q_out__;
  q = q_out__;
  qstore = ALLOC_N(doublereal, (ldqs)*(n));
  work = ALLOC_N(doublereal, (((icompq == 0) || (icompq == 1)) ? 1 + 3*n + 2*n*LG(n) + 2*pow(n,2) : icompq == 2 ? 4*n + pow(n,2) : 0));
  iwork = ALLOC_N(integer, (((icompq == 0) || (icompq == 1)) ? 6 + 6*n + 5*n*LG(n) : icompq == 2 ? 3 + 5*n : 0));

  dlaed0_(&icompq, &qsiz, &n, d, e, q, &ldq, qstore, &ldqs, work, iwork, &info);

  free(qstore);
  free(work);
  free(iwork);
  rb_info = INT2NUM(info);
  return rb_ary_new3(3, rb_info, rb_d, rb_q);
}

void
init_lapack_dlaed0(VALUE mLapack){
  rb_define_module_function(mLapack, "dlaed0", rb_dlaed0, -1);
}
