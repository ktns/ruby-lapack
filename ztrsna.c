#include "rb_lapack.h"

extern VOID ztrsna_(char *job, char *howmny, logical *select, integer *n, doublecomplex *t, integer *ldt, doublecomplex *vl, integer *ldvl, doublecomplex *vr, integer *ldvr, doublereal *s, doublereal *sep, integer *mm, integer *m, doublecomplex *work, integer *ldwork, doublereal *rwork, integer *info);

static VALUE
rb_ztrsna(int argc, VALUE *argv, VALUE self){
  VALUE rb_job;
  char job; 
  VALUE rb_howmny;
  char howmny; 
  VALUE rb_select;
  logical *select; 
  VALUE rb_t;
  doublecomplex *t; 
  VALUE rb_vl;
  doublecomplex *vl; 
  VALUE rb_vr;
  doublecomplex *vr; 
  VALUE rb_ldwork;
  integer ldwork; 
  VALUE rb_s;
  doublereal *s; 
  VALUE rb_sep;
  doublereal *sep; 
  VALUE rb_m;
  integer m; 
  VALUE rb_info;
  integer info; 
  doublecomplex *work;
  doublereal *rwork;

  integer n;
  integer ldt;
  integer ldvl;
  integer ldvr;
  integer mm;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  s, sep, m, info = NumRu::Lapack.ztrsna( job, howmny, select, t, vl, vr, ldwork)\n    or\n  NumRu::Lapack.ztrsna  # print help\n\n\nFORTRAN MANUAL\n      SUBROUTINE ZTRSNA( JOB, HOWMNY, SELECT, N, T, LDT, VL, LDVL, VR, LDVR, S, SEP, MM, M, WORK, LDWORK, RWORK, INFO )\n\n*  Purpose\n*  =======\n*\n*  ZTRSNA estimates reciprocal condition numbers for specified\n*  eigenvalues and/or right eigenvectors of a complex upper triangular\n*  matrix T (or of any matrix Q*T*Q**H with Q unitary).\n*\n\n*  Arguments\n*  =========\n*\n*  JOB     (input) CHARACTER*1\n*          Specifies whether condition numbers are required for\n*          eigenvalues (S) or eigenvectors (SEP):\n*          = 'E': for eigenvalues only (S);\n*          = 'V': for eigenvectors only (SEP);\n*          = 'B': for both eigenvalues and eigenvectors (S and SEP).\n*\n*  HOWMNY  (input) CHARACTER*1\n*          = 'A': compute condition numbers for all eigenpairs;\n*          = 'S': compute condition numbers for selected eigenpairs\n*                 specified by the array SELECT.\n*\n*  SELECT  (input) LOGICAL array, dimension (N)\n*          If HOWMNY = 'S', SELECT specifies the eigenpairs for which\n*          condition numbers are required. To select condition numbers\n*          for the j-th eigenpair, SELECT(j) must be set to .TRUE..\n*          If HOWMNY = 'A', SELECT is not referenced.\n*\n*  N       (input) INTEGER\n*          The order of the matrix T. N >= 0.\n*\n*  T       (input) COMPLEX*16 array, dimension (LDT,N)\n*          The upper triangular matrix T.\n*\n*  LDT     (input) INTEGER\n*          The leading dimension of the array T. LDT >= max(1,N).\n*\n*  VL      (input) COMPLEX*16 array, dimension (LDVL,M)\n*          If JOB = 'E' or 'B', VL must contain left eigenvectors of T\n*          (or of any Q*T*Q**H with Q unitary), corresponding to the\n*          eigenpairs specified by HOWMNY and SELECT. The eigenvectors\n*          must be stored in consecutive columns of VL, as returned by\n*          ZHSEIN or ZTREVC.\n*          If JOB = 'V', VL is not referenced.\n*\n*  LDVL    (input) INTEGER\n*          The leading dimension of the array VL.\n*          LDVL >= 1; and if JOB = 'E' or 'B', LDVL >= N.\n*\n*  VR      (input) COMPLEX*16 array, dimension (LDVR,M)\n*          If JOB = 'E' or 'B', VR must contain right eigenvectors of T\n*          (or of any Q*T*Q**H with Q unitary), corresponding to the\n*          eigenpairs specified by HOWMNY and SELECT. The eigenvectors\n*          must be stored in consecutive columns of VR, as returned by\n*          ZHSEIN or ZTREVC.\n*          If JOB = 'V', VR is not referenced.\n*\n*  LDVR    (input) INTEGER\n*          The leading dimension of the array VR.\n*          LDVR >= 1; and if JOB = 'E' or 'B', LDVR >= N.\n*\n*  S       (output) DOUBLE PRECISION array, dimension (MM)\n*          If JOB = 'E' or 'B', the reciprocal condition numbers of the\n*          selected eigenvalues, stored in consecutive elements of the\n*          array. Thus S(j), SEP(j), and the j-th columns of VL and VR\n*          all correspond to the same eigenpair (but not in general the\n*          j-th eigenpair, unless all eigenpairs are selected).\n*          If JOB = 'V', S is not referenced.\n*\n*  SEP     (output) DOUBLE PRECISION array, dimension (MM)\n*          If JOB = 'V' or 'B', the estimated reciprocal condition\n*          numbers of the selected eigenvectors, stored in consecutive\n*          elements of the array.\n*          If JOB = 'E', SEP is not referenced.\n*\n*  MM      (input) INTEGER\n*          The number of elements in the arrays S (if JOB = 'E' or 'B')\n*           and/or SEP (if JOB = 'V' or 'B'). MM >= M.\n*\n*  M       (output) INTEGER\n*          The number of elements of the arrays S and/or SEP actually\n*          used to store the estimated condition numbers.\n*          If HOWMNY = 'A', M is set to N.\n*\n*  WORK    (workspace) COMPLEX*16 array, dimension (LDWORK,N+6)\n*          If JOB = 'E', WORK is not referenced.\n*\n*  LDWORK  (input) INTEGER\n*          The leading dimension of the array WORK.\n*          LDWORK >= 1; and if JOB = 'V' or 'B', LDWORK >= N.\n*\n*  RWORK   (workspace) DOUBLE PRECISION array, dimension (N)\n*          If JOB = 'E', RWORK is not referenced.\n*\n*  INFO    (output) INTEGER\n*          = 0: successful exit\n*          < 0: if INFO = -i, the i-th argument had an illegal value\n*\n\n*  Further Details\n*  ===============\n*\n*  The reciprocal of the condition number of an eigenvalue lambda is\n*  defined as\n*\n*          S(lambda) = |v'*u| / (norm(u)*norm(v))\n*\n*  where u and v are the right and left eigenvectors of T corresponding\n*  to lambda; v' denotes the conjugate transpose of v, and norm(u)\n*  denotes the Euclidean norm. These reciprocal condition numbers always\n*  lie between zero (very badly conditioned) and one (very well\n*  conditioned). If n = 1, S(lambda) is defined to be 1.\n*\n*  An approximate error bound for a computed eigenvalue W(i) is given by\n*\n*                      EPS * norm(T) / S(i)\n*\n*  where EPS is the machine precision.\n*\n*  The reciprocal of the condition number of the right eigenvector u\n*  corresponding to lambda is defined as follows. Suppose\n*\n*              T = ( lambda  c  )\n*                  (   0    T22 )\n*\n*  Then the reciprocal condition number is\n*\n*          SEP( lambda, T22 ) = sigma-min( T22 - lambda*I )\n*\n*  where sigma-min denotes the smallest singular value. We approximate\n*  the smallest singular value by the reciprocal of an estimate of the\n*  one-norm of the inverse of T22 - lambda*I. If n = 1, SEP(1) is\n*  defined to be abs(T(1,1)).\n*\n*  An approximate error bound for a computed right eigenvector VR(i)\n*  is given by\n*\n*                      EPS * norm(T) / SEP(i)\n*\n*  =====================================================================\n*\n\n");
    return Qnil;
  }
  if (argc != 7)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 7)", argc);
  rb_job = argv[0];
  rb_howmny = argv[1];
  rb_select = argv[2];
  rb_t = argv[3];
  rb_vl = argv[4];
  rb_vr = argv[5];
  rb_ldwork = argv[6];

  howmny = StringValueCStr(rb_howmny)[0];
  if (!NA_IsNArray(rb_t))
    rb_raise(rb_eArgError, "t (4th argument) must be NArray");
  if (NA_RANK(rb_t) != 2)
    rb_raise(rb_eArgError, "rank of t (4th argument) must be %d", 2);
  n = NA_SHAPE1(rb_t);
  ldt = NA_SHAPE0(rb_t);
  if (NA_TYPE(rb_t) != NA_DCOMPLEX)
    rb_t = na_change_type(rb_t, NA_DCOMPLEX);
  t = NA_PTR_TYPE(rb_t, doublecomplex*);
  if (!NA_IsNArray(rb_vl))
    rb_raise(rb_eArgError, "vl (5th argument) must be NArray");
  if (NA_RANK(rb_vl) != 2)
    rb_raise(rb_eArgError, "rank of vl (5th argument) must be %d", 2);
  m = NA_SHAPE1(rb_vl);
  ldvl = NA_SHAPE0(rb_vl);
  if (NA_TYPE(rb_vl) != NA_DCOMPLEX)
    rb_vl = na_change_type(rb_vl, NA_DCOMPLEX);
  vl = NA_PTR_TYPE(rb_vl, doublecomplex*);
  if (!NA_IsNArray(rb_vr))
    rb_raise(rb_eArgError, "vr (6th argument) must be NArray");
  if (NA_RANK(rb_vr) != 2)
    rb_raise(rb_eArgError, "rank of vr (6th argument) must be %d", 2);
  if (NA_SHAPE1(rb_vr) != m)
    rb_raise(rb_eRuntimeError, "shape 1 of vr must be the same as shape 1 of vl");
  ldvr = NA_SHAPE0(rb_vr);
  if (NA_TYPE(rb_vr) != NA_DCOMPLEX)
    rb_vr = na_change_type(rb_vr, NA_DCOMPLEX);
  vr = NA_PTR_TYPE(rb_vr, doublecomplex*);
  if (!NA_IsNArray(rb_select))
    rb_raise(rb_eArgError, "select (3th argument) must be NArray");
  if (NA_RANK(rb_select) != 1)
    rb_raise(rb_eArgError, "rank of select (3th argument) must be %d", 1);
  if (NA_SHAPE0(rb_select) != n)
    rb_raise(rb_eRuntimeError, "shape 0 of select must be the same as shape 1 of t");
  if (NA_TYPE(rb_select) != NA_LINT)
    rb_select = na_change_type(rb_select, NA_LINT);
  select = NA_PTR_TYPE(rb_select, logical*);
  job = StringValueCStr(rb_job)[0];
  ldwork = ((lsame_(&job,"V")) || (lsame_(&job,"B"))) ? n : 1;
  mm = m;
  {
    int shape[1];
    shape[0] = mm;
    rb_s = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  s = NA_PTR_TYPE(rb_s, doublereal*);
  {
    int shape[1];
    shape[0] = mm;
    rb_sep = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  sep = NA_PTR_TYPE(rb_sep, doublereal*);
  work = ALLOC_N(doublecomplex, (lsame_(&job,"E") ? 0 : ldwork)*(lsame_(&job,"E") ? 0 : n+6));
  rwork = ALLOC_N(doublereal, (lsame_(&job,"E") ? 0 : n));

  ztrsna_(&job, &howmny, select, &n, t, &ldt, vl, &ldvl, vr, &ldvr, s, sep, &mm, &m, work, &ldwork, rwork, &info);

  free(work);
  free(rwork);
  rb_m = INT2NUM(m);
  rb_info = INT2NUM(info);
  return rb_ary_new3(4, rb_s, rb_sep, rb_m, rb_info);
}

void
init_lapack_ztrsna(VALUE mLapack){
  rb_define_module_function(mLapack, "ztrsna", rb_ztrsna, -1);
}
