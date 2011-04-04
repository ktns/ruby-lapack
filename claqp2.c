#include "rb_lapack.h"

extern VOID claqp2_(integer *m, integer *n, integer *offset, complex *a, integer *lda, integer *jpvt, complex *tau, real *vn1, real *vn2, complex *work);

static VALUE
rb_claqp2(int argc, VALUE *argv, VALUE self){
  VALUE rb_m;
  integer m; 
  VALUE rb_offset;
  integer offset; 
  VALUE rb_a;
  complex *a; 
  VALUE rb_jpvt;
  integer *jpvt; 
  VALUE rb_vn1;
  real *vn1; 
  VALUE rb_vn2;
  real *vn2; 
  VALUE rb_tau;
  complex *tau; 
  VALUE rb_a_out__;
  complex *a_out__;
  VALUE rb_jpvt_out__;
  integer *jpvt_out__;
  VALUE rb_vn1_out__;
  real *vn1_out__;
  VALUE rb_vn2_out__;
  real *vn2_out__;
  complex *work;

  integer lda;
  integer n;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  tau, a, jpvt, vn1, vn2 = NumRu::Lapack.claqp2( m, offset, a, jpvt, vn1, vn2)\n    or\n  NumRu::Lapack.claqp2  # print help\n\n\nFORTRAN MANUAL\n      SUBROUTINE CLAQP2( M, N, OFFSET, A, LDA, JPVT, TAU, VN1, VN2, WORK )\n\n*  Purpose\n*  =======\n*\n*  CLAQP2 computes a QR factorization with column pivoting of\n*  the block A(OFFSET+1:M,1:N).\n*  The block A(1:OFFSET,1:N) is accordingly pivoted, but not factorized.\n*\n\n*  Arguments\n*  =========\n*\n*  M       (input) INTEGER\n*          The number of rows of the matrix A. M >= 0.\n*\n*  N       (input) INTEGER\n*          The number of columns of the matrix A. N >= 0.\n*\n*  OFFSET  (input) INTEGER\n*          The number of rows of the matrix A that must be pivoted\n*          but no factorized. OFFSET >= 0.\n*\n*  A       (input/output) COMPLEX array, dimension (LDA,N)\n*          On entry, the M-by-N matrix A.\n*          On exit, the upper triangle of block A(OFFSET+1:M,1:N) is \n*          the triangular factor obtained; the elements in block\n*          A(OFFSET+1:M,1:N) below the diagonal, together with the\n*          array TAU, represent the orthogonal matrix Q as a product of\n*          elementary reflectors. Block A(1:OFFSET,1:N) has been\n*          accordingly pivoted, but no factorized.\n*\n*  LDA     (input) INTEGER\n*          The leading dimension of the array A. LDA >= max(1,M).\n*\n*  JPVT    (input/output) INTEGER array, dimension (N)\n*          On entry, if JPVT(i) .ne. 0, the i-th column of A is permuted\n*          to the front of A*P (a leading column); if JPVT(i) = 0,\n*          the i-th column of A is a free column.\n*          On exit, if JPVT(i) = k, then the i-th column of A*P\n*          was the k-th column of A.\n*\n*  TAU     (output) COMPLEX array, dimension (min(M,N))\n*          The scalar factors of the elementary reflectors.\n*\n*  VN1     (input/output) REAL array, dimension (N)\n*          The vector with the partial column norms.\n*\n*  VN2     (input/output) REAL array, dimension (N)\n*          The vector with the exact column norms.\n*\n*  WORK    (workspace) COMPLEX array, dimension (N)\n*\n\n*  Further Details\n*  ===============\n*\n*  Based on contributions by\n*    G. Quintana-Orti, Depto. de Informatica, Universidad Jaime I, Spain\n*    X. Sun, Computer Science Dept., Duke University, USA\n*\n*  Partial column norm updating strategy modified by\n*    Z. Drmac and Z. Bujanovic, Dept. of Mathematics,\n*    University of Zagreb, Croatia.\n*     June 2010\n*  For more details see LAPACK Working Note 176.\n*  =====================================================================\n*\n\n");
    return Qnil;
  }
  if (argc != 6)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 6)", argc);
  rb_m = argv[0];
  rb_offset = argv[1];
  rb_a = argv[2];
  rb_jpvt = argv[3];
  rb_vn1 = argv[4];
  rb_vn2 = argv[5];

  if (!NA_IsNArray(rb_a))
    rb_raise(rb_eArgError, "a (3th argument) must be NArray");
  if (NA_RANK(rb_a) != 2)
    rb_raise(rb_eArgError, "rank of a (3th argument) must be %d", 2);
  n = NA_SHAPE1(rb_a);
  lda = NA_SHAPE0(rb_a);
  if (NA_TYPE(rb_a) != NA_SCOMPLEX)
    rb_a = na_change_type(rb_a, NA_SCOMPLEX);
  a = NA_PTR_TYPE(rb_a, complex*);
  m = NUM2INT(rb_m);
  if (!NA_IsNArray(rb_vn1))
    rb_raise(rb_eArgError, "vn1 (5th argument) must be NArray");
  if (NA_RANK(rb_vn1) != 1)
    rb_raise(rb_eArgError, "rank of vn1 (5th argument) must be %d", 1);
  if (NA_SHAPE0(rb_vn1) != n)
    rb_raise(rb_eRuntimeError, "shape 0 of vn1 must be the same as shape 1 of a");
  if (NA_TYPE(rb_vn1) != NA_SFLOAT)
    rb_vn1 = na_change_type(rb_vn1, NA_SFLOAT);
  vn1 = NA_PTR_TYPE(rb_vn1, real*);
  if (!NA_IsNArray(rb_vn2))
    rb_raise(rb_eArgError, "vn2 (6th argument) must be NArray");
  if (NA_RANK(rb_vn2) != 1)
    rb_raise(rb_eArgError, "rank of vn2 (6th argument) must be %d", 1);
  if (NA_SHAPE0(rb_vn2) != n)
    rb_raise(rb_eRuntimeError, "shape 0 of vn2 must be the same as shape 1 of a");
  if (NA_TYPE(rb_vn2) != NA_SFLOAT)
    rb_vn2 = na_change_type(rb_vn2, NA_SFLOAT);
  vn2 = NA_PTR_TYPE(rb_vn2, real*);
  if (!NA_IsNArray(rb_jpvt))
    rb_raise(rb_eArgError, "jpvt (4th argument) must be NArray");
  if (NA_RANK(rb_jpvt) != 1)
    rb_raise(rb_eArgError, "rank of jpvt (4th argument) must be %d", 1);
  if (NA_SHAPE0(rb_jpvt) != n)
    rb_raise(rb_eRuntimeError, "shape 0 of jpvt must be the same as shape 1 of a");
  if (NA_TYPE(rb_jpvt) != NA_LINT)
    rb_jpvt = na_change_type(rb_jpvt, NA_LINT);
  jpvt = NA_PTR_TYPE(rb_jpvt, integer*);
  offset = NUM2INT(rb_offset);
  {
    int shape[1];
    shape[0] = MIN(m,n);
    rb_tau = na_make_object(NA_SCOMPLEX, 1, shape, cNArray);
  }
  tau = NA_PTR_TYPE(rb_tau, complex*);
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
  {
    int shape[1];
    shape[0] = n;
    rb_vn1_out__ = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  vn1_out__ = NA_PTR_TYPE(rb_vn1_out__, real*);
  MEMCPY(vn1_out__, vn1, real, NA_TOTAL(rb_vn1));
  rb_vn1 = rb_vn1_out__;
  vn1 = vn1_out__;
  {
    int shape[1];
    shape[0] = n;
    rb_vn2_out__ = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  vn2_out__ = NA_PTR_TYPE(rb_vn2_out__, real*);
  MEMCPY(vn2_out__, vn2, real, NA_TOTAL(rb_vn2));
  rb_vn2 = rb_vn2_out__;
  vn2 = vn2_out__;
  work = ALLOC_N(complex, (n));

  claqp2_(&m, &n, &offset, a, &lda, jpvt, tau, vn1, vn2, work);

  free(work);
  return rb_ary_new3(5, rb_tau, rb_a, rb_jpvt, rb_vn1, rb_vn2);
}

void
init_lapack_claqp2(VALUE mLapack){
  rb_define_module_function(mLapack, "claqp2", rb_claqp2, -1);
}
