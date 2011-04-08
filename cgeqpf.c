#include "rb_lapack.h"

extern VOID cgeqpf_(integer *m, integer *n, complex *a, integer *lda, integer *jpvt, complex *tau, complex *work, real *rwork, integer *info);

static VALUE sHelp, sUsage;

static VALUE
rblapack_cgeqpf(int argc, VALUE *argv, VALUE self){
  VALUE rblapack_m;
  integer m; 
  VALUE rblapack_a;
  complex *a; 
  VALUE rblapack_jpvt;
  integer *jpvt; 
  VALUE rblapack_tau;
  complex *tau; 
  VALUE rblapack_info;
  integer info; 
  VALUE rblapack_a_out__;
  complex *a_out__;
  VALUE rblapack_jpvt_out__;
  integer *jpvt_out__;
  complex *work;
  real *rwork;

  integer lda;
  integer n;

  VALUE rb_options;
  if (argc > 0 && TYPE(argv[argc-1]) == T_HASH) {
    argc--;
    rb_options = argv[argc];
    if (rb_hash_aref(rb_options, sHelp) == Qtrue) {
      printf("%s\n", "USAGE:\n  tau, info, a, jpvt = NumRu::Lapack.cgeqpf( m, a, jpvt, [:usage => usage, :help => help])\n\n\nFORTRAN MANUAL\n      SUBROUTINE CGEQPF( M, N, A, LDA, JPVT, TAU, WORK, RWORK, INFO )\n\n*  Purpose\n*  =======\n*\n*  This routine is deprecated and has been replaced by routine CGEQP3.\n*\n*  CGEQPF computes a QR factorization with column pivoting of a\n*  complex M-by-N matrix A: A*P = Q*R.\n*\n\n*  Arguments\n*  =========\n*\n*  M       (input) INTEGER\n*          The number of rows of the matrix A. M >= 0.\n*\n*  N       (input) INTEGER\n*          The number of columns of the matrix A. N >= 0\n*\n*  A       (input/output) COMPLEX array, dimension (LDA,N)\n*          On entry, the M-by-N matrix A.\n*          On exit, the upper triangle of the array contains the\n*          min(M,N)-by-N upper triangular matrix R; the elements\n*          below the diagonal, together with the array TAU,\n*          represent the unitary matrix Q as a product of\n*          min(m,n) elementary reflectors.\n*\n*  LDA     (input) INTEGER\n*          The leading dimension of the array A. LDA >= max(1,M).\n*\n*  JPVT    (input/output) INTEGER array, dimension (N)\n*          On entry, if JPVT(i) .ne. 0, the i-th column of A is permuted\n*          to the front of A*P (a leading column); if JPVT(i) = 0,\n*          the i-th column of A is a free column.\n*          On exit, if JPVT(i) = k, then the i-th column of A*P\n*          was the k-th column of A.\n*\n*  TAU     (output) COMPLEX array, dimension (min(M,N))\n*          The scalar factors of the elementary reflectors.\n*\n*  WORK    (workspace) COMPLEX array, dimension (N)\n*\n*  RWORK   (workspace) REAL array, dimension (2*N)\n*\n*  INFO    (output) INTEGER\n*          = 0:  successful exit\n*          < 0:  if INFO = -i, the i-th argument had an illegal value\n*\n\n*  Further Details\n*  ===============\n*\n*  The matrix Q is represented as a product of elementary reflectors\n*\n*     Q = H(1) H(2) . . . H(n)\n*\n*  Each H(i) has the form\n*\n*     H = I - tau * v * v'\n*\n*  where tau is a complex scalar, and v is a complex vector with\n*  v(1:i-1) = 0 and v(i) = 1; v(i+1:m) is stored on exit in A(i+1:m,i).\n*\n*  The matrix P is represented in jpvt as follows: If\n*     jpvt(j) = i\n*  then the jth column of P is the ith canonical unit vector.\n*\n*  Partial column norm updating strategy modified by\n*    Z. Drmac and Z. Bujanovic, Dept. of Mathematics,\n*    University of Zagreb, Croatia.\n*     June 2010\n*  For more details see LAPACK Working Note 176.\n*\n*  =====================================================================\n*\n\n");
      rb_exit(0);
    }
    if (rb_hash_aref(rb_options, sUsage) == Qtrue) {
      printf("%s\n", "USAGE:\n  tau, info, a, jpvt = NumRu::Lapack.cgeqpf( m, a, jpvt, [:usage => usage, :help => help])\n");
      rb_exit(0);
    } 
  } else
    rb_options = Qnil;
  if (argc != 3)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 3)", argc);
  rblapack_m = argv[0];
  rblapack_a = argv[1];
  rblapack_jpvt = argv[2];
  if (rb_options != Qnil) {
  }

  if (!NA_IsNArray(rblapack_a))
    rb_raise(rb_eArgError, "a (2th argument) must be NArray");
  if (NA_RANK(rblapack_a) != 2)
    rb_raise(rb_eArgError, "rank of a (2th argument) must be %d", 2);
  n = NA_SHAPE1(rblapack_a);
  lda = NA_SHAPE0(rblapack_a);
  if (NA_TYPE(rblapack_a) != NA_SCOMPLEX)
    rblapack_a = na_change_type(rblapack_a, NA_SCOMPLEX);
  a = NA_PTR_TYPE(rblapack_a, complex*);
  m = NUM2INT(rblapack_m);
  if (!NA_IsNArray(rblapack_jpvt))
    rb_raise(rb_eArgError, "jpvt (3th argument) must be NArray");
  if (NA_RANK(rblapack_jpvt) != 1)
    rb_raise(rb_eArgError, "rank of jpvt (3th argument) must be %d", 1);
  if (NA_SHAPE0(rblapack_jpvt) != n)
    rb_raise(rb_eRuntimeError, "shape 0 of jpvt must be the same as shape 1 of a");
  if (NA_TYPE(rblapack_jpvt) != NA_LINT)
    rblapack_jpvt = na_change_type(rblapack_jpvt, NA_LINT);
  jpvt = NA_PTR_TYPE(rblapack_jpvt, integer*);
  {
    int shape[1];
    shape[0] = MIN(m,n);
    rblapack_tau = na_make_object(NA_SCOMPLEX, 1, shape, cNArray);
  }
  tau = NA_PTR_TYPE(rblapack_tau, complex*);
  {
    int shape[2];
    shape[0] = lda;
    shape[1] = n;
    rblapack_a_out__ = na_make_object(NA_SCOMPLEX, 2, shape, cNArray);
  }
  a_out__ = NA_PTR_TYPE(rblapack_a_out__, complex*);
  MEMCPY(a_out__, a, complex, NA_TOTAL(rblapack_a));
  rblapack_a = rblapack_a_out__;
  a = a_out__;
  {
    int shape[1];
    shape[0] = n;
    rblapack_jpvt_out__ = na_make_object(NA_LINT, 1, shape, cNArray);
  }
  jpvt_out__ = NA_PTR_TYPE(rblapack_jpvt_out__, integer*);
  MEMCPY(jpvt_out__, jpvt, integer, NA_TOTAL(rblapack_jpvt));
  rblapack_jpvt = rblapack_jpvt_out__;
  jpvt = jpvt_out__;
  work = ALLOC_N(complex, (n));
  rwork = ALLOC_N(real, (2*n));

  cgeqpf_(&m, &n, a, &lda, jpvt, tau, work, rwork, &info);

  free(work);
  free(rwork);
  rblapack_info = INT2NUM(info);
  return rb_ary_new3(4, rblapack_tau, rblapack_info, rblapack_a, rblapack_jpvt);
}

void
init_lapack_cgeqpf(VALUE mLapack){
  rb_define_module_function(mLapack, "cgeqpf", rblapack_cgeqpf, -1);
  sHelp = ID2SYM(rb_intern("help"));
  sUsage = ID2SYM(rb_intern("usage"));
}
