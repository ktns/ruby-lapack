#include "rb_lapack.h"

static VALUE
rb_dgeev(int argc, VALUE *argv, VALUE self){
  VALUE rb_jobvl;
  char jobvl; 
  VALUE rb_jobvr;
  char jobvr; 
  VALUE rb_a;
  doublereal *a; 
  VALUE rb_lwork;
  integer lwork; 
  VALUE rb_wr;
  doublereal *wr; 
  VALUE rb_wi;
  doublereal *wi; 
  VALUE rb_vl;
  doublereal *vl; 
  VALUE rb_vr;
  doublereal *vr; 
  VALUE rb_work;
  doublereal *work; 
  VALUE rb_info;
  integer info; 
  VALUE rb_a_out__;
  doublereal *a_out__;

  integer lda;
  integer n;
  integer ldvl;
  integer ldvr;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  wr, wi, vl, vr, work, info, a = NumRu::Lapack.dgeev( jobvl, jobvr, a, lwork)\n    or\n  NumRu::Lapack.dgeev  # print help\n\n\nFORTRAN MANUAL\n      SUBROUTINE DGEEV( JOBVL, JOBVR, N, A, LDA, WR, WI, VL, LDVL, VR, LDVR, WORK, LWORK, INFO )\n\n*  Purpose\n*  =======\n*\n*  DGEEV computes for an N-by-N real nonsymmetric matrix A, the\n*  eigenvalues and, optionally, the left and/or right eigenvectors.\n*\n*  The right eigenvector v(j) of A satisfies\n*                   A * v(j) = lambda(j) * v(j)\n*  where lambda(j) is its eigenvalue.\n*  The left eigenvector u(j) of A satisfies\n*                u(j)**H * A = lambda(j) * u(j)**H\n*  where u(j)**H denotes the conjugate transpose of u(j).\n*\n*  The computed eigenvectors are normalized to have Euclidean norm\n*  equal to 1 and largest component real.\n*\n\n*  Arguments\n*  =========\n*\n*  JOBVL   (input) CHARACTER*1\n*          = 'N': left eigenvectors of A are not computed;\n*          = 'V': left eigenvectors of A are computed.\n*\n*  JOBVR   (input) CHARACTER*1\n*          = 'N': right eigenvectors of A are not computed;\n*          = 'V': right eigenvectors of A are computed.\n*\n*  N       (input) INTEGER\n*          The order of the matrix A. N >= 0.\n*\n*  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)\n*          On entry, the N-by-N matrix A.\n*          On exit, A has been overwritten.\n*\n*  LDA     (input) INTEGER\n*          The leading dimension of the array A.  LDA >= max(1,N).\n*\n*  WR      (output) DOUBLE PRECISION array, dimension (N)\n*  WI      (output) DOUBLE PRECISION array, dimension (N)\n*          WR and WI contain the real and imaginary parts,\n*          respectively, of the computed eigenvalues.  Complex\n*          conjugate pairs of eigenvalues appear consecutively\n*          with the eigenvalue having the positive imaginary part\n*          first.\n*\n*  VL      (output) DOUBLE PRECISION array, dimension (LDVL,N)\n*          If JOBVL = 'V', the left eigenvectors u(j) are stored one\n*          after another in the columns of VL, in the same order\n*          as their eigenvalues.\n*          If JOBVL = 'N', VL is not referenced.\n*          If the j-th eigenvalue is real, then u(j) = VL(:,j),\n*          the j-th column of VL.\n*          If the j-th and (j+1)-st eigenvalues form a complex\n*          conjugate pair, then u(j) = VL(:,j) + i*VL(:,j+1) and\n*          u(j+1) = VL(:,j) - i*VL(:,j+1).\n*\n*  LDVL    (input) INTEGER\n*          The leading dimension of the array VL.  LDVL >= 1; if\n*          JOBVL = 'V', LDVL >= N.\n*\n*  VR      (output) DOUBLE PRECISION array, dimension (LDVR,N)\n*          If JOBVR = 'V', the right eigenvectors v(j) are stored one\n*          after another in the columns of VR, in the same order\n*          as their eigenvalues.\n*          If JOBVR = 'N', VR is not referenced.\n*          If the j-th eigenvalue is real, then v(j) = VR(:,j),\n*          the j-th column of VR.\n*          If the j-th and (j+1)-st eigenvalues form a complex\n*          conjugate pair, then v(j) = VR(:,j) + i*VR(:,j+1) and\n*          v(j+1) = VR(:,j) - i*VR(:,j+1).\n*\n*  LDVR    (input) INTEGER\n*          The leading dimension of the array VR.  LDVR >= 1; if\n*          JOBVR = 'V', LDVR >= N.\n*\n*  WORK    (workspace/output) DOUBLE PRECISION array, dimension (MAX(1,LWORK))\n*          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.\n*\n*  LWORK   (input) INTEGER\n*          The dimension of the array WORK.  LWORK >= max(1,3*N), and\n*          if JOBVL = 'V' or JOBVR = 'V', LWORK >= 4*N.  For good\n*          performance, LWORK must generally be larger.\n*\n*          If LWORK = -1, then a workspace query is assumed; the routine\n*          only calculates the optimal size of the WORK array, returns\n*          this value as the first entry of the WORK array, and no error\n*          message related to LWORK is issued by XERBLA.\n*\n*  INFO    (output) INTEGER\n*          = 0:  successful exit\n*          < 0:  if INFO = -i, the i-th argument had an illegal value.\n*          > 0:  if INFO = i, the QR algorithm failed to compute all the\n*                eigenvalues, and no eigenvectors have been computed;\n*                elements i+1:N of WR and WI contain eigenvalues which\n*                have converged.\n*\n\n*  =====================================================================\n*\n\n");
    return Qnil;
  }
  if (argc != 4)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 4)", argc);
  rb_jobvl = argv[0];
  rb_jobvr = argv[1];
  rb_a = argv[2];
  rb_lwork = argv[3];

  jobvl = StringValueCStr(rb_jobvl)[0];
  jobvr = StringValueCStr(rb_jobvr)[0];
  lwork = NUM2INT(rb_lwork);
  if (!NA_IsNArray(rb_a))
    rb_raise(rb_eArgError, "a (3th argument) must be NArray");
  if (NA_RANK(rb_a) != 2)
    rb_raise(rb_eArgError, "rank of a (3th argument) must be %d", 2);
  lda = NA_SHAPE0(rb_a);
  n = NA_SHAPE1(rb_a);
  if (NA_TYPE(rb_a) != NA_DFLOAT)
    rb_a = na_change_type(rb_a, NA_DFLOAT);
  a = NA_PTR_TYPE(rb_a, doublereal*);
  {
    int shape[1];
    shape[0] = n;
    rb_wr = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  wr = NA_PTR_TYPE(rb_wr, doublereal*);
  {
    int shape[1];
    shape[0] = n;
    rb_wi = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  wi = NA_PTR_TYPE(rb_wi, doublereal*);
  ldvl = lsame_(&jobvl,"V") ? n : 1;
  {
    int shape[2];
    shape[0] = ldvl;
    shape[1] = n;
    rb_vl = na_make_object(NA_DFLOAT, 2, shape, cNArray);
  }
  vl = NA_PTR_TYPE(rb_vl, doublereal*);
  ldvr = lsame_(&jobvr,"V") ? n : 1;
  {
    int shape[2];
    shape[0] = ldvr;
    shape[1] = n;
    rb_vr = na_make_object(NA_DFLOAT, 2, shape, cNArray);
  }
  vr = NA_PTR_TYPE(rb_vr, doublereal*);
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

  dgeev_(&jobvl, &jobvr, &n, a, &lda, wr, wi, vl, &ldvl, vr, &ldvr, work, &lwork, &info);

  rb_info = INT2NUM(info);
  return rb_ary_new3(7, rb_wr, rb_wi, rb_vl, rb_vr, rb_work, rb_info, rb_a);
}

void
init_lapack_dgeev(VALUE mLapack){
  rb_define_module_function(mLapack, "dgeev", rb_dgeev, -1);
}
