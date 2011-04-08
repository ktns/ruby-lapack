#include "rb_lapack.h"

extern VOID slaed0_(integer *icompq, integer *qsiz, integer *n, real *d, real *e, real *q, integer *ldq, real *qstore, integer *ldqs, real *work, integer *iwork, integer *info);

static VALUE sHelp, sUsage;

static VALUE
rblapack_slaed0(int argc, VALUE *argv, VALUE self){
  VALUE rblapack_icompq;
  integer icompq; 
  VALUE rblapack_qsiz;
  integer qsiz; 
  VALUE rblapack_d;
  real *d; 
  VALUE rblapack_e;
  real *e; 
  VALUE rblapack_q;
  real *q; 
  VALUE rblapack_info;
  integer info; 
  VALUE rblapack_d_out__;
  real *d_out__;
  VALUE rblapack_q_out__;
  real *q_out__;
  real *qstore;
  real *work;
  integer *iwork;

  integer n;
  integer ldq;
  integer ldqs;

  VALUE rb_options;
  if (argc > 0 && TYPE(argv[argc-1]) == T_HASH) {
    argc--;
    rb_options = argv[argc];
    if (rb_hash_aref(rb_options, sHelp) == Qtrue) {
      printf("%s\n", "USAGE:\n  info, d, q = NumRu::Lapack.slaed0( icompq, qsiz, d, e, q, [:usage => usage, :help => help])\n\n\nFORTRAN MANUAL\n      SUBROUTINE SLAED0( ICOMPQ, QSIZ, N, D, E, Q, LDQ, QSTORE, LDQS, WORK, IWORK, INFO )\n\n*  Purpose\n*  =======\n*\n*  SLAED0 computes all eigenvalues and corresponding eigenvectors of a\n*  symmetric tridiagonal matrix using the divide and conquer method.\n*\n\n*  Arguments\n*  =========\n*\n*  ICOMPQ  (input) INTEGER\n*          = 0:  Compute eigenvalues only.\n*          = 1:  Compute eigenvectors of original dense symmetric matrix\n*                also.  On entry, Q contains the orthogonal matrix used\n*                to reduce the original matrix to tridiagonal form.\n*          = 2:  Compute eigenvalues and eigenvectors of tridiagonal\n*                matrix.\n*\n*  QSIZ   (input) INTEGER\n*         The dimension of the orthogonal matrix used to reduce\n*         the full matrix to tridiagonal form.  QSIZ >= N if ICOMPQ = 1.\n*\n*  N      (input) INTEGER\n*         The dimension of the symmetric tridiagonal matrix.  N >= 0.\n*\n*  D      (input/output) REAL array, dimension (N)\n*         On entry, the main diagonal of the tridiagonal matrix.\n*         On exit, its eigenvalues.\n*\n*  E      (input) REAL array, dimension (N-1)\n*         The off-diagonal elements of the tridiagonal matrix.\n*         On exit, E has been destroyed.\n*\n*  Q      (input/output) REAL array, dimension (LDQ, N)\n*         On entry, Q must contain an N-by-N orthogonal matrix.\n*         If ICOMPQ = 0    Q is not referenced.\n*         If ICOMPQ = 1    On entry, Q is a subset of the columns of the\n*                          orthogonal matrix used to reduce the full\n*                          matrix to tridiagonal form corresponding to\n*                          the subset of the full matrix which is being\n*                          decomposed at this time.\n*         If ICOMPQ = 2    On entry, Q will be the identity matrix.\n*                          On exit, Q contains the eigenvectors of the\n*                          tridiagonal matrix.\n*\n*  LDQ    (input) INTEGER\n*         The leading dimension of the array Q.  If eigenvectors are\n*         desired, then  LDQ >= max(1,N).  In any case,  LDQ >= 1.\n*\n*  QSTORE (workspace) REAL array, dimension (LDQS, N)\n*         Referenced only when ICOMPQ = 1.  Used to store parts of\n*         the eigenvector matrix when the updating matrix multiplies\n*         take place.\n*\n*  LDQS   (input) INTEGER\n*         The leading dimension of the array QSTORE.  If ICOMPQ = 1,\n*         then  LDQS >= max(1,N).  In any case,  LDQS >= 1.\n*\n*  WORK   (workspace) REAL array,\n*         If ICOMPQ = 0 or 1, the dimension of WORK must be at least\n*                     1 + 3*N + 2*N*lg N + 2*N**2\n*                     ( lg( N ) = smallest integer k\n*                                 such that 2^k >= N )\n*         If ICOMPQ = 2, the dimension of WORK must be at least\n*                     4*N + N**2.\n*\n*  IWORK  (workspace) INTEGER array,\n*         If ICOMPQ = 0 or 1, the dimension of IWORK must be at least\n*                        6 + 6*N + 5*N*lg N.\n*                        ( lg( N ) = smallest integer k\n*                                    such that 2^k >= N )\n*         If ICOMPQ = 2, the dimension of IWORK must be at least\n*                        3 + 5*N.\n*\n*  INFO   (output) INTEGER\n*          = 0:  successful exit.\n*          < 0:  if INFO = -i, the i-th argument had an illegal value.\n*          > 0:  The algorithm failed to compute an eigenvalue while\n*                working on the submatrix lying in rows and columns\n*                INFO/(N+1) through mod(INFO,N+1).\n*\n\n*  Further Details\n*  ===============\n*\n*  Based on contributions by\n*     Jeff Rutter, Computer Science Division, University of California\n*     at Berkeley, USA\n*\n*  =====================================================================\n*\n\n");
      rb_exit(0);
    }
    if (rb_hash_aref(rb_options, sUsage) == Qtrue) {
      printf("%s\n", "USAGE:\n  info, d, q = NumRu::Lapack.slaed0( icompq, qsiz, d, e, q, [:usage => usage, :help => help])\n");
      rb_exit(0);
    } 
  } else
    rb_options = Qnil;
  if (argc != 5)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 5)", argc);
  rblapack_icompq = argv[0];
  rblapack_qsiz = argv[1];
  rblapack_d = argv[2];
  rblapack_e = argv[3];
  rblapack_q = argv[4];
  if (rb_options != Qnil) {
  }

  qsiz = NUM2INT(rblapack_qsiz);
  if (!NA_IsNArray(rblapack_d))
    rb_raise(rb_eArgError, "d (3th argument) must be NArray");
  if (NA_RANK(rblapack_d) != 1)
    rb_raise(rb_eArgError, "rank of d (3th argument) must be %d", 1);
  n = NA_SHAPE0(rblapack_d);
  if (NA_TYPE(rblapack_d) != NA_SFLOAT)
    rblapack_d = na_change_type(rblapack_d, NA_SFLOAT);
  d = NA_PTR_TYPE(rblapack_d, real*);
  if (!NA_IsNArray(rblapack_q))
    rb_raise(rb_eArgError, "q (5th argument) must be NArray");
  if (NA_RANK(rblapack_q) != 2)
    rb_raise(rb_eArgError, "rank of q (5th argument) must be %d", 2);
  if (NA_SHAPE1(rblapack_q) != n)
    rb_raise(rb_eRuntimeError, "shape 1 of q must be the same as shape 0 of d");
  ldq = NA_SHAPE0(rblapack_q);
  if (NA_TYPE(rblapack_q) != NA_SFLOAT)
    rblapack_q = na_change_type(rblapack_q, NA_SFLOAT);
  q = NA_PTR_TYPE(rblapack_q, real*);
  icompq = NUM2INT(rblapack_icompq);
  if (!NA_IsNArray(rblapack_e))
    rb_raise(rb_eArgError, "e (4th argument) must be NArray");
  if (NA_RANK(rblapack_e) != 1)
    rb_raise(rb_eArgError, "rank of e (4th argument) must be %d", 1);
  if (NA_SHAPE0(rblapack_e) != (n-1))
    rb_raise(rb_eRuntimeError, "shape 0 of e must be %d", n-1);
  if (NA_TYPE(rblapack_e) != NA_SFLOAT)
    rblapack_e = na_change_type(rblapack_e, NA_SFLOAT);
  e = NA_PTR_TYPE(rblapack_e, real*);
  ldqs = icompq == 1 ? MAX(1,n) : 1;
  {
    int shape[1];
    shape[0] = n;
    rblapack_d_out__ = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  d_out__ = NA_PTR_TYPE(rblapack_d_out__, real*);
  MEMCPY(d_out__, d, real, NA_TOTAL(rblapack_d));
  rblapack_d = rblapack_d_out__;
  d = d_out__;
  {
    int shape[2];
    shape[0] = ldq;
    shape[1] = n;
    rblapack_q_out__ = na_make_object(NA_SFLOAT, 2, shape, cNArray);
  }
  q_out__ = NA_PTR_TYPE(rblapack_q_out__, real*);
  MEMCPY(q_out__, q, real, NA_TOTAL(rblapack_q));
  rblapack_q = rblapack_q_out__;
  q = q_out__;
  qstore = ALLOC_N(real, (ldqs)*(n));
  work = ALLOC_N(real, (((icompq == 0) || (icompq == 1)) ? 1 + 3*n + 2*n*LG(n) + 2*pow(n,2) : icompq == 2 ? 4*n + pow(n,2) : 0));
  iwork = ALLOC_N(integer, (((icompq == 0) || (icompq == 1)) ? 6 + 6*n + 5*n*LG(n) : icompq == 2 ? 3 + 5*n : 0));

  slaed0_(&icompq, &qsiz, &n, d, e, q, &ldq, qstore, &ldqs, work, iwork, &info);

  free(qstore);
  free(work);
  free(iwork);
  rblapack_info = INT2NUM(info);
  return rb_ary_new3(3, rblapack_info, rblapack_d, rblapack_q);
}

void
init_lapack_slaed0(VALUE mLapack){
  rb_define_module_function(mLapack, "slaed0", rblapack_slaed0, -1);
  sHelp = ID2SYM(rb_intern("help"));
  sUsage = ID2SYM(rb_intern("usage"));
}
