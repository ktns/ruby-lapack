#include "rb_lapack.h"

extern VOID slaqp2_(integer *m, integer *n, integer *offset, real *a, integer *lda, integer *jpvt, real *tau, real *vn1, real *vn2, real *work);

static VALUE sHelp, sUsage;

static VALUE
rblapack_slaqp2(int argc, VALUE *argv, VALUE self){
  VALUE rblapack_m;
  integer m; 
  VALUE rblapack_offset;
  integer offset; 
  VALUE rblapack_a;
  real *a; 
  VALUE rblapack_jpvt;
  integer *jpvt; 
  VALUE rblapack_vn1;
  real *vn1; 
  VALUE rblapack_vn2;
  real *vn2; 
  VALUE rblapack_tau;
  real *tau; 
  VALUE rblapack_a_out__;
  real *a_out__;
  VALUE rblapack_jpvt_out__;
  integer *jpvt_out__;
  VALUE rblapack_vn1_out__;
  real *vn1_out__;
  VALUE rblapack_vn2_out__;
  real *vn2_out__;
  real *work;

  integer lda;
  integer n;

  VALUE rb_options;
  if (argc > 0 && TYPE(argv[argc-1]) == T_HASH) {
    argc--;
    rb_options = argv[argc];
    if (rb_hash_aref(rb_options, sHelp) == Qtrue) {
      printf("%s\n", "USAGE:\n  tau, a, jpvt, vn1, vn2 = NumRu::Lapack.slaqp2( m, offset, a, jpvt, vn1, vn2, [:usage => usage, :help => help])\n\n\nFORTRAN MANUAL\n      SUBROUTINE SLAQP2( M, N, OFFSET, A, LDA, JPVT, TAU, VN1, VN2, WORK )\n\n*  Purpose\n*  =======\n*\n*  SLAQP2 computes a QR factorization with column pivoting of\n*  the block A(OFFSET+1:M,1:N).\n*  The block A(1:OFFSET,1:N) is accordingly pivoted, but not factorized.\n*\n\n*  Arguments\n*  =========\n*\n*  M       (input) INTEGER\n*          The number of rows of the matrix A. M >= 0.\n*\n*  N       (input) INTEGER\n*          The number of columns of the matrix A. N >= 0.\n*\n*  OFFSET  (input) INTEGER\n*          The number of rows of the matrix A that must be pivoted\n*          but no factorized. OFFSET >= 0.\n*\n*  A       (input/output) REAL array, dimension (LDA,N)\n*          On entry, the M-by-N matrix A.\n*          On exit, the upper triangle of block A(OFFSET+1:M,1:N) is \n*          the triangular factor obtained; the elements in block \n*          A(OFFSET+1:M,1:N) below the diagonal, together with the \n*          array TAU, represent the orthogonal matrix Q as a product of\n*          elementary reflectors. Block A(1:OFFSET,1:N) has been\n*          accordingly pivoted, but no factorized.\n*\n*  LDA     (input) INTEGER\n*          The leading dimension of the array A. LDA >= max(1,M).\n*\n*  JPVT    (input/output) INTEGER array, dimension (N)\n*          On entry, if JPVT(i) .ne. 0, the i-th column of A is permuted\n*          to the front of A*P (a leading column); if JPVT(i) = 0,\n*          the i-th column of A is a free column.\n*          On exit, if JPVT(i) = k, then the i-th column of A*P\n*          was the k-th column of A.\n*\n*  TAU     (output) REAL array, dimension (min(M,N))\n*          The scalar factors of the elementary reflectors.\n*\n*  VN1     (input/output) REAL array, dimension (N)\n*          The vector with the partial column norms.\n*\n*  VN2     (input/output) REAL array, dimension (N)\n*          The vector with the exact column norms.\n*\n*  WORK    (workspace) REAL array, dimension (N)\n*\n\n*  Further Details\n*  ===============\n*\n*  Based on contributions by\n*    G. Quintana-Orti, Depto. de Informatica, Universidad Jaime I, Spain\n*    X. Sun, Computer Science Dept., Duke University, USA\n*\n*  Partial column norm updating strategy modified by\n*    Z. Drmac and Z. Bujanovic, Dept. of Mathematics,\n*    University of Zagreb, Croatia.\n*     June 2010\n*  For more details see LAPACK Working Note 176.\n*  =====================================================================\n*\n\n");
      rb_exit(0);
    }
    if (rb_hash_aref(rb_options, sUsage) == Qtrue) {
      printf("%s\n", "USAGE:\n  tau, a, jpvt, vn1, vn2 = NumRu::Lapack.slaqp2( m, offset, a, jpvt, vn1, vn2, [:usage => usage, :help => help])\n");
      rb_exit(0);
    } 
  } else
    rb_options = Qnil;
  if (argc != 6)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 6)", argc);
  rblapack_m = argv[0];
  rblapack_offset = argv[1];
  rblapack_a = argv[2];
  rblapack_jpvt = argv[3];
  rblapack_vn1 = argv[4];
  rblapack_vn2 = argv[5];
  if (rb_options != Qnil) {
  }

  if (!NA_IsNArray(rblapack_a))
    rb_raise(rb_eArgError, "a (3th argument) must be NArray");
  if (NA_RANK(rblapack_a) != 2)
    rb_raise(rb_eArgError, "rank of a (3th argument) must be %d", 2);
  n = NA_SHAPE1(rblapack_a);
  lda = NA_SHAPE0(rblapack_a);
  if (NA_TYPE(rblapack_a) != NA_SFLOAT)
    rblapack_a = na_change_type(rblapack_a, NA_SFLOAT);
  a = NA_PTR_TYPE(rblapack_a, real*);
  m = NUM2INT(rblapack_m);
  if (!NA_IsNArray(rblapack_vn1))
    rb_raise(rb_eArgError, "vn1 (5th argument) must be NArray");
  if (NA_RANK(rblapack_vn1) != 1)
    rb_raise(rb_eArgError, "rank of vn1 (5th argument) must be %d", 1);
  if (NA_SHAPE0(rblapack_vn1) != n)
    rb_raise(rb_eRuntimeError, "shape 0 of vn1 must be the same as shape 1 of a");
  if (NA_TYPE(rblapack_vn1) != NA_SFLOAT)
    rblapack_vn1 = na_change_type(rblapack_vn1, NA_SFLOAT);
  vn1 = NA_PTR_TYPE(rblapack_vn1, real*);
  if (!NA_IsNArray(rblapack_vn2))
    rb_raise(rb_eArgError, "vn2 (6th argument) must be NArray");
  if (NA_RANK(rblapack_vn2) != 1)
    rb_raise(rb_eArgError, "rank of vn2 (6th argument) must be %d", 1);
  if (NA_SHAPE0(rblapack_vn2) != n)
    rb_raise(rb_eRuntimeError, "shape 0 of vn2 must be the same as shape 1 of a");
  if (NA_TYPE(rblapack_vn2) != NA_SFLOAT)
    rblapack_vn2 = na_change_type(rblapack_vn2, NA_SFLOAT);
  vn2 = NA_PTR_TYPE(rblapack_vn2, real*);
  if (!NA_IsNArray(rblapack_jpvt))
    rb_raise(rb_eArgError, "jpvt (4th argument) must be NArray");
  if (NA_RANK(rblapack_jpvt) != 1)
    rb_raise(rb_eArgError, "rank of jpvt (4th argument) must be %d", 1);
  if (NA_SHAPE0(rblapack_jpvt) != n)
    rb_raise(rb_eRuntimeError, "shape 0 of jpvt must be the same as shape 1 of a");
  if (NA_TYPE(rblapack_jpvt) != NA_LINT)
    rblapack_jpvt = na_change_type(rblapack_jpvt, NA_LINT);
  jpvt = NA_PTR_TYPE(rblapack_jpvt, integer*);
  offset = NUM2INT(rblapack_offset);
  {
    int shape[1];
    shape[0] = MIN(m,n);
    rblapack_tau = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  tau = NA_PTR_TYPE(rblapack_tau, real*);
  {
    int shape[2];
    shape[0] = lda;
    shape[1] = n;
    rblapack_a_out__ = na_make_object(NA_SFLOAT, 2, shape, cNArray);
  }
  a_out__ = NA_PTR_TYPE(rblapack_a_out__, real*);
  MEMCPY(a_out__, a, real, NA_TOTAL(rblapack_a));
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
  {
    int shape[1];
    shape[0] = n;
    rblapack_vn1_out__ = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  vn1_out__ = NA_PTR_TYPE(rblapack_vn1_out__, real*);
  MEMCPY(vn1_out__, vn1, real, NA_TOTAL(rblapack_vn1));
  rblapack_vn1 = rblapack_vn1_out__;
  vn1 = vn1_out__;
  {
    int shape[1];
    shape[0] = n;
    rblapack_vn2_out__ = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  vn2_out__ = NA_PTR_TYPE(rblapack_vn2_out__, real*);
  MEMCPY(vn2_out__, vn2, real, NA_TOTAL(rblapack_vn2));
  rblapack_vn2 = rblapack_vn2_out__;
  vn2 = vn2_out__;
  work = ALLOC_N(real, (n));

  slaqp2_(&m, &n, &offset, a, &lda, jpvt, tau, vn1, vn2, work);

  free(work);
  return rb_ary_new3(5, rblapack_tau, rblapack_a, rblapack_jpvt, rblapack_vn1, rblapack_vn2);
}

void
init_lapack_slaqp2(VALUE mLapack){
  rb_define_module_function(mLapack, "slaqp2", rblapack_slaqp2, -1);
  sHelp = ID2SYM(rb_intern("help"));
  sUsage = ID2SYM(rb_intern("usage"));
}
