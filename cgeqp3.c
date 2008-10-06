#include "rb_lapack.h"

static VALUE
rb_cgeqp3(int argc, VALUE *argv, VALUE self){
  VALUE rb_m;
  integer m; 
  VALUE rb_a;
  complex *a; 
  VALUE rb_jpvt;
  integer *jpvt; 
  VALUE rb_lwork;
  integer lwork; 
  VALUE rb_tau;
  complex *tau; 
  VALUE rb_work;
  complex *work; 
  VALUE rb_info;
  integer info; 
  VALUE rb_a_out__;
  complex *a_out__;
  VALUE rb_jpvt_out__;
  integer *jpvt_out__;
  real *rwork;

  integer lda;
  integer n;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  tau, work, info, a, jpvt = NumRu::Lapack.cgeqp3( m, a, jpvt, lwork)\n    or\n  NumRu::Lapack.cgeqp3  # print help\n\n\nFORTRAN MANUAL\n      SUBROUTINE CGEQP3( M, N, A, LDA, JPVT, TAU, WORK, LWORK, RWORK, INFO )\n\n*  Purpose\n*  =======\n*\n*  CGEQP3 computes a QR factorization with column pivoting of a\n*  matrix A:  A*P = Q*R  using Level 3 BLAS.\n*\n\n*  Arguments\n*  =========\n*\n*  M       (input) INTEGER\n*          The number of rows of the matrix A. M >= 0.\n*\n*  N       (input) INTEGER\n*          The number of columns of the matrix A.  N >= 0.\n*\n*  A       (input/output) COMPLEX array, dimension (LDA,N)\n*          On entry, the M-by-N matrix A.\n*          On exit, the upper triangle of the array contains the\n*          min(M,N)-by-N upper trapezoidal matrix R; the elements below\n*          the diagonal, together with the array TAU, represent the\n*          unitary matrix Q as a product of min(M,N) elementary\n*          reflectors.\n*\n*  LDA     (input) INTEGER\n*          The leading dimension of the array A. LDA >= max(1,M).\n*\n*  JPVT    (input/output) INTEGER array, dimension (N)\n*          On entry, if JPVT(J).ne.0, the J-th column of A is permuted\n*          to the front of A*P (a leading column); if JPVT(J)=0,\n*          the J-th column of A is a free column.\n*          On exit, if JPVT(J)=K, then the J-th column of A*P was the\n*          the K-th column of A.\n*\n*  TAU     (output) COMPLEX array, dimension (min(M,N))\n*          The scalar factors of the elementary reflectors.\n*\n*  WORK    (workspace/output) COMPLEX array, dimension (MAX(1,LWORK))\n*          On exit, if INFO=0, WORK(1) returns the optimal LWORK.\n*\n*  LWORK   (input) INTEGER\n*          The dimension of the array WORK. LWORK >= N+1.\n*          For optimal performance LWORK >= ( N+1 )*NB, where NB\n*          is the optimal blocksize.\n*\n*          If LWORK = -1, then a workspace query is assumed; the routine\n*          only calculates the optimal size of the WORK array, returns\n*          this value as the first entry of the WORK array, and no error\n*          message related to LWORK is issued by XERBLA.\n*\n*  RWORK   (workspace) REAL array, dimension (2*N)\n*\n*  INFO    (output) INTEGER\n*          = 0: successful exit.\n*          < 0: if INFO = -i, the i-th argument had an illegal value.\n*\n\n*  Further Details\n*  ===============\n*\n*  The matrix Q is represented as a product of elementary reflectors\n*\n*     Q = H(1) H(2) . . . H(k), where k = min(m,n).\n*\n*  Each H(i) has the form\n*\n*     H(i) = I - tau * v * v'\n*\n*  where tau is a real/complex scalar, and v is a real/complex vector\n*  with v(1:i-1) = 0 and v(i) = 1; v(i+1:m) is stored on exit in\n*  A(i+1:m,i), and tau in TAU(i).\n*\n*  Based on contributions by\n*    G. Quintana-Orti, Depto. de Informatica, Universidad Jaime I, Spain\n*    X. Sun, Computer Science Dept., Duke University, USA\n*\n*  =====================================================================\n*\n\n");
    return Qnil;
  }
  if (argc != 4)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 4)", argc);
  rb_m = argv[0];
  rb_a = argv[1];
  rb_jpvt = argv[2];
  rb_lwork = argv[3];

  m = NUM2INT(rb_m);
  lwork = NUM2INT(rb_lwork);
  if (!NA_IsNArray(rb_a))
    rb_raise(rb_eArgError, "a (2th argument) must be NArray");
  if (NA_RANK(rb_a) != 2)
    rb_raise(rb_eArgError, "rank of a (2th argument) must be %d", 2);
  lda = NA_SHAPE0(rb_a);
  n = NA_SHAPE1(rb_a);
  if (NA_TYPE(rb_a) != NA_SCOMPLEX)
    rb_a = na_change_type(rb_a, NA_SCOMPLEX);
  a = NA_PTR_TYPE(rb_a, complex*);
  if (!NA_IsNArray(rb_jpvt))
    rb_raise(rb_eArgError, "jpvt (3th argument) must be NArray");
  if (NA_RANK(rb_jpvt) != 1)
    rb_raise(rb_eArgError, "rank of jpvt (3th argument) must be %d", 1);
  if (NA_SHAPE0(rb_jpvt) != n)
    rb_raise(rb_eRuntimeError, "shape 0 of jpvt must be the same as shape 1 of a");
  if (NA_TYPE(rb_jpvt) != NA_LINT)
    rb_jpvt = na_change_type(rb_jpvt, NA_LINT);
  jpvt = NA_PTR_TYPE(rb_jpvt, integer*);
  {
    int shape[1];
    shape[0] = MIN(m,n);
    rb_tau = na_make_object(NA_SCOMPLEX, 1, shape, cNArray);
  }
  tau = NA_PTR_TYPE(rb_tau, complex*);
  {
    int shape[1];
    shape[0] = MAX(1,lwork);
    rb_work = na_make_object(NA_SCOMPLEX, 1, shape, cNArray);
  }
  work = NA_PTR_TYPE(rb_work, complex*);
  {
    int shape[2];
    shape[0] = lda;
    shape[1] = n;
    rb_a_out__ = na_make_object(NA_SCOMPLEX, 2, shape, cNArray);
  }
  a_out__ = NA_PTR_TYPE(rb_a_out__, complex*);
  MEMCPY(a_out__, a, complex, NA_TOTAL(rb_a));
  rb_a = rb_a_out__;
  a = a_out__;
  {
    int shape[1];
    shape[0] = n;
    rb_jpvt_out__ = na_make_object(NA_LINT, 1, shape, cNArray);
  }
  jpvt_out__ = NA_PTR_TYPE(rb_jpvt_out__, integer*);
  MEMCPY(jpvt_out__, jpvt, integer, NA_TOTAL(rb_jpvt));
  rb_jpvt = rb_jpvt_out__;
  jpvt = jpvt_out__;
  rwork = ALLOC_N(real, (2*n));

  cgeqp3_(&m, &n, a, &lda, jpvt, tau, work, &lwork, rwork, &info);

  free(rwork);
  rb_info = INT2NUM(info);
  return rb_ary_new3(5, rb_tau, rb_work, rb_info, rb_a, rb_jpvt);
}

void
init_lapack_cgeqp3(VALUE mLapack){
  rb_define_module_function(mLapack, "cgeqp3", rb_cgeqp3, -1);
}
