#include "rb_lapack.h"

extern VOID zgeev_(char *jobvl, char *jobvr, integer *n, doublecomplex *a, integer *lda, doublecomplex *w, doublecomplex *vl, integer *ldvl, doublecomplex *vr, integer *ldvr, doublecomplex *work, integer *lwork, doublereal *rwork, integer *info);

static VALUE sHelp, sUsage;

static VALUE
rblapack_zgeev(int argc, VALUE *argv, VALUE self){
  VALUE rblapack_jobvl;
  char jobvl; 
  VALUE rblapack_jobvr;
  char jobvr; 
  VALUE rblapack_a;
  doublecomplex *a; 
  VALUE rblapack_lwork;
  integer lwork; 
  VALUE rblapack_w;
  doublecomplex *w; 
  VALUE rblapack_vl;
  doublecomplex *vl; 
  VALUE rblapack_vr;
  doublecomplex *vr; 
  VALUE rblapack_work;
  doublecomplex *work; 
  VALUE rblapack_info;
  integer info; 
  VALUE rblapack_a_out__;
  doublecomplex *a_out__;
  doublereal *rwork;

  integer lda;
  integer n;
  integer ldvl;
  integer ldvr;

  VALUE rb_options;
  if (argc > 0 && TYPE(argv[argc-1]) == T_HASH) {
    argc--;
    rb_options = argv[argc];
    if (rb_hash_aref(rb_options, sHelp) == Qtrue) {
      printf("%s\n", "USAGE:\n  w, vl, vr, work, info, a = NumRu::Lapack.zgeev( jobvl, jobvr, a, lwork, [:usage => usage, :help => help])\n\n\nFORTRAN MANUAL\n      SUBROUTINE ZGEEV( JOBVL, JOBVR, N, A, LDA, W, VL, LDVL, VR, LDVR, WORK, LWORK, RWORK, INFO )\n\n*  Purpose\n*  =======\n*\n*  ZGEEV computes for an N-by-N complex nonsymmetric matrix A, the\n*  eigenvalues and, optionally, the left and/or right eigenvectors.\n*\n*  The right eigenvector v(j) of A satisfies\n*                   A * v(j) = lambda(j) * v(j)\n*  where lambda(j) is its eigenvalue.\n*  The left eigenvector u(j) of A satisfies\n*                u(j)**H * A = lambda(j) * u(j)**H\n*  where u(j)**H denotes the conjugate transpose of u(j).\n*\n*  The computed eigenvectors are normalized to have Euclidean norm\n*  equal to 1 and largest component real.\n*\n\n*  Arguments\n*  =========\n*\n*  JOBVL   (input) CHARACTER*1\n*          = 'N': left eigenvectors of A are not computed;\n*          = 'V': left eigenvectors of are computed.\n*\n*  JOBVR   (input) CHARACTER*1\n*          = 'N': right eigenvectors of A are not computed;\n*          = 'V': right eigenvectors of A are computed.\n*\n*  N       (input) INTEGER\n*          The order of the matrix A. N >= 0.\n*\n*  A       (input/output) COMPLEX*16 array, dimension (LDA,N)\n*          On entry, the N-by-N matrix A.\n*          On exit, A has been overwritten.\n*\n*  LDA     (input) INTEGER\n*          The leading dimension of the array A.  LDA >= max(1,N).\n*\n*  W       (output) COMPLEX*16 array, dimension (N)\n*          W contains the computed eigenvalues.\n*\n*  VL      (output) COMPLEX*16 array, dimension (LDVL,N)\n*          If JOBVL = 'V', the left eigenvectors u(j) are stored one\n*          after another in the columns of VL, in the same order\n*          as their eigenvalues.\n*          If JOBVL = 'N', VL is not referenced.\n*          u(j) = VL(:,j), the j-th column of VL.\n*\n*  LDVL    (input) INTEGER\n*          The leading dimension of the array VL.  LDVL >= 1; if\n*          JOBVL = 'V', LDVL >= N.\n*\n*  VR      (output) COMPLEX*16 array, dimension (LDVR,N)\n*          If JOBVR = 'V', the right eigenvectors v(j) are stored one\n*          after another in the columns of VR, in the same order\n*          as their eigenvalues.\n*          If JOBVR = 'N', VR is not referenced.\n*          v(j) = VR(:,j), the j-th column of VR.\n*\n*  LDVR    (input) INTEGER\n*          The leading dimension of the array VR.  LDVR >= 1; if\n*          JOBVR = 'V', LDVR >= N.\n*\n*  WORK    (workspace/output) COMPLEX*16 array, dimension (MAX(1,LWORK))\n*          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.\n*\n*  LWORK   (input) INTEGER\n*          The dimension of the array WORK.  LWORK >= max(1,2*N).\n*          For good performance, LWORK must generally be larger.\n*\n*          If LWORK = -1, then a workspace query is assumed; the routine\n*          only calculates the optimal size of the WORK array, returns\n*          this value as the first entry of the WORK array, and no error\n*          message related to LWORK is issued by XERBLA.\n*\n*  RWORK   (workspace) DOUBLE PRECISION array, dimension (2*N)\n*\n*  INFO    (output) INTEGER\n*          = 0:  successful exit\n*          < 0:  if INFO = -i, the i-th argument had an illegal value.\n*          > 0:  if INFO = i, the QR algorithm failed to compute all the\n*                eigenvalues, and no eigenvectors have been computed;\n*                elements and i+1:N of W contain eigenvalues which have\n*                converged.\n*\n\n*  =====================================================================\n*\n\n");
      rb_exit(0);
    }
    if (rb_hash_aref(rb_options, sUsage) == Qtrue) {
      printf("%s\n", "USAGE:\n  w, vl, vr, work, info, a = NumRu::Lapack.zgeev( jobvl, jobvr, a, lwork, [:usage => usage, :help => help])\n");
      rb_exit(0);
    } 
  } else
    rb_options = Qnil;
  if (argc != 4)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 4)", argc);
  rblapack_jobvl = argv[0];
  rblapack_jobvr = argv[1];
  rblapack_a = argv[2];
  rblapack_lwork = argv[3];
  if (rb_options != Qnil) {
  }

  if (!NA_IsNArray(rblapack_a))
    rb_raise(rb_eArgError, "a (3th argument) must be NArray");
  if (NA_RANK(rblapack_a) != 2)
    rb_raise(rb_eArgError, "rank of a (3th argument) must be %d", 2);
  n = NA_SHAPE1(rblapack_a);
  lda = NA_SHAPE0(rblapack_a);
  if (NA_TYPE(rblapack_a) != NA_DCOMPLEX)
    rblapack_a = na_change_type(rblapack_a, NA_DCOMPLEX);
  a = NA_PTR_TYPE(rblapack_a, doublecomplex*);
  jobvr = StringValueCStr(rblapack_jobvr)[0];
  jobvl = StringValueCStr(rblapack_jobvl)[0];
  lwork = NUM2INT(rblapack_lwork);
  ldvr = lsame_(&jobvr,"V") ? n : 1;
  ldvl = lsame_(&jobvl,"V") ? n : 1;
  {
    int shape[1];
    shape[0] = n;
    rblapack_w = na_make_object(NA_DCOMPLEX, 1, shape, cNArray);
  }
  w = NA_PTR_TYPE(rblapack_w, doublecomplex*);
  {
    int shape[2];
    shape[0] = ldvl;
    shape[1] = n;
    rblapack_vl = na_make_object(NA_DCOMPLEX, 2, shape, cNArray);
  }
  vl = NA_PTR_TYPE(rblapack_vl, doublecomplex*);
  {
    int shape[2];
    shape[0] = ldvr;
    shape[1] = n;
    rblapack_vr = na_make_object(NA_DCOMPLEX, 2, shape, cNArray);
  }
  vr = NA_PTR_TYPE(rblapack_vr, doublecomplex*);
  {
    int shape[1];
    shape[0] = MAX(1,lwork);
    rblapack_work = na_make_object(NA_DCOMPLEX, 1, shape, cNArray);
  }
  work = NA_PTR_TYPE(rblapack_work, doublecomplex*);
  {
    int shape[2];
    shape[0] = lda;
    shape[1] = n;
    rblapack_a_out__ = na_make_object(NA_DCOMPLEX, 2, shape, cNArray);
  }
  a_out__ = NA_PTR_TYPE(rblapack_a_out__, doublecomplex*);
  MEMCPY(a_out__, a, doublecomplex, NA_TOTAL(rblapack_a));
  rblapack_a = rblapack_a_out__;
  a = a_out__;
  rwork = ALLOC_N(doublereal, (2*n));

  zgeev_(&jobvl, &jobvr, &n, a, &lda, w, vl, &ldvl, vr, &ldvr, work, &lwork, rwork, &info);

  free(rwork);
  rblapack_info = INT2NUM(info);
  return rb_ary_new3(6, rblapack_w, rblapack_vl, rblapack_vr, rblapack_work, rblapack_info, rblapack_a);
}

void
init_lapack_zgeev(VALUE mLapack){
  rb_define_module_function(mLapack, "zgeev", rblapack_zgeev, -1);
  sHelp = ID2SYM(rb_intern("help"));
  sUsage = ID2SYM(rb_intern("usage"));
}
