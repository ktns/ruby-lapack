#include "rb_lapack.h"

extern VOID dgetri_(integer *n, doublereal *a, integer *lda, integer *ipiv, doublereal *work, integer *lwork, integer *info);

static VALUE
rb_dgetri(int argc, VALUE *argv, VALUE self){
  VALUE rb_a;
  doublereal *a; 
  VALUE rb_ipiv;
  integer *ipiv; 
  VALUE rb_lwork;
  integer lwork; 
  VALUE rb_work;
  doublereal *work; 
  VALUE rb_info;
  integer info; 
  VALUE rb_a_out__;
  doublereal *a_out__;

  integer lda;
  integer n;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  work, info, a = NumRu::Lapack.dgetri( a, ipiv, lwork)\n    or\n  NumRu::Lapack.dgetri  # print help\n\n\nFORTRAN MANUAL\n      SUBROUTINE DGETRI( N, A, LDA, IPIV, WORK, LWORK, INFO )\n\n*  Purpose\n*  =======\n*\n*  DGETRI computes the inverse of a matrix using the LU factorization\n*  computed by DGETRF.\n*\n*  This method inverts U and then computes inv(A) by solving the system\n*  inv(A)*L = inv(U) for inv(A).\n*\n\n*  Arguments\n*  =========\n*\n*  N       (input) INTEGER\n*          The order of the matrix A.  N >= 0.\n*\n*  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)\n*          On entry, the factors L and U from the factorization\n*          A = P*L*U as computed by DGETRF.\n*          On exit, if INFO = 0, the inverse of the original matrix A.\n*\n*  LDA     (input) INTEGER\n*          The leading dimension of the array A.  LDA >= max(1,N).\n*\n*  IPIV    (input) INTEGER array, dimension (N)\n*          The pivot indices from DGETRF; for 1<=i<=N, row i of the\n*          matrix was interchanged with row IPIV(i).\n*\n*  WORK    (workspace/output) DOUBLE PRECISION array, dimension (MAX(1,LWORK))\n*          On exit, if INFO=0, then WORK(1) returns the optimal LWORK.\n*\n*  LWORK   (input) INTEGER\n*          The dimension of the array WORK.  LWORK >= max(1,N).\n*          For optimal performance LWORK >= N*NB, where NB is\n*          the optimal blocksize returned by ILAENV.\n*\n*          If LWORK = -1, then a workspace query is assumed; the routine\n*          only calculates the optimal size of the WORK array, returns\n*          this value as the first entry of the WORK array, and no error\n*          message related to LWORK is issued by XERBLA.\n*\n*  INFO    (output) INTEGER\n*          = 0:  successful exit\n*          < 0:  if INFO = -i, the i-th argument had an illegal value\n*          > 0:  if INFO = i, U(i,i) is exactly zero; the matrix is\n*                singular and its inverse could not be computed.\n*\n\n*  =====================================================================\n*\n\n");
    return Qnil;
  }
  if (argc != 3)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 3)", argc);
  rb_a = argv[0];
  rb_ipiv = argv[1];
  rb_lwork = argv[2];

  if (!NA_IsNArray(rb_ipiv))
    rb_raise(rb_eArgError, "ipiv (2th argument) must be NArray");
  if (NA_RANK(rb_ipiv) != 1)
    rb_raise(rb_eArgError, "rank of ipiv (2th argument) must be %d", 1);
  n = NA_SHAPE0(rb_ipiv);
  if (NA_TYPE(rb_ipiv) != NA_LINT)
    rb_ipiv = na_change_type(rb_ipiv, NA_LINT);
  ipiv = NA_PTR_TYPE(rb_ipiv, integer*);
  if (!NA_IsNArray(rb_a))
    rb_raise(rb_eArgError, "a (1th argument) must be NArray");
  if (NA_RANK(rb_a) != 2)
    rb_raise(rb_eArgError, "rank of a (1th argument) must be %d", 2);
  if (NA_SHAPE1(rb_a) != n)
    rb_raise(rb_eRuntimeError, "shape 1 of a must be the same as shape 0 of ipiv");
  lda = NA_SHAPE0(rb_a);
  if (NA_TYPE(rb_a) != NA_DFLOAT)
    rb_a = na_change_type(rb_a, NA_DFLOAT);
  a = NA_PTR_TYPE(rb_a, doublereal*);
  lwork = NUM2INT(rb_lwork);
  {
    int shape[1];
    shape[0] = MAX(1,lwork);
    rb_work = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  work = NA_PTR_TYPE(rb_work, doublereal*);
  {
    int shape[2];
    shape[0] = lda;
    shape[1] = n;
    rb_a_out__ = na_make_object(NA_DFLOAT, 2, shape, cNArray);
  }
  a_out__ = NA_PTR_TYPE(rb_a_out__, doublereal*);
  MEMCPY(a_out__, a, doublereal, NA_TOTAL(rb_a));
  rb_a = rb_a_out__;
  a = a_out__;

  dgetri_(&n, a, &lda, ipiv, work, &lwork, &info);

  rb_info = INT2NUM(info);
  return rb_ary_new3(3, rb_work, rb_info, rb_a);
}

void
init_lapack_dgetri(VALUE mLapack){
  rb_define_module_function(mLapack, "dgetri", rb_dgetri, -1);
}
