#include "rb_lapack.h"

extern VOID sstevx_(char *jobz, char *range, integer *n, real *d, real *e, real *vl, real *vu, integer *il, integer *iu, real *abstol, integer *m, real *w, real *z, integer *ldz, real *work, integer *iwork, integer *ifail, integer *info);

static VALUE
rb_sstevx(int argc, VALUE *argv, VALUE self){
  VALUE rb_jobz;
  char jobz; 
  VALUE rb_range;
  char range; 
  VALUE rb_d;
  real *d; 
  VALUE rb_e;
  real *e; 
  VALUE rb_vl;
  real vl; 
  VALUE rb_vu;
  real vu; 
  VALUE rb_il;
  integer il; 
  VALUE rb_iu;
  integer iu; 
  VALUE rb_abstol;
  real abstol; 
  VALUE rb_m;
  integer m; 
  VALUE rb_w;
  real *w; 
  VALUE rb_z;
  real *z; 
  VALUE rb_ifail;
  integer *ifail; 
  VALUE rb_info;
  integer info; 
  VALUE rb_d_out__;
  real *d_out__;
  VALUE rb_e_out__;
  real *e_out__;
  real *work;
  integer *iwork;

  integer n;
  integer ldz;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  m, w, z, ifail, info, d, e = NumRu::Lapack.sstevx( jobz, range, d, e, vl, vu, il, iu, abstol)\n    or\n  NumRu::Lapack.sstevx  # print help\n\n\nFORTRAN MANUAL\n      SUBROUTINE SSTEVX( JOBZ, RANGE, N, D, E, VL, VU, IL, IU, ABSTOL, M, W, Z, LDZ, WORK, IWORK, IFAIL, INFO )\n\n*  Purpose\n*  =======\n*\n*  SSTEVX computes selected eigenvalues and, optionally, eigenvectors\n*  of a real symmetric tridiagonal matrix A.  Eigenvalues and\n*  eigenvectors can be selected by specifying either a range of values\n*  or a range of indices for the desired eigenvalues.\n*\n\n*  Arguments\n*  =========\n*\n*  JOBZ    (input) CHARACTER*1\n*          = 'N':  Compute eigenvalues only;\n*          = 'V':  Compute eigenvalues and eigenvectors.\n*\n*  RANGE   (input) CHARACTER*1\n*          = 'A': all eigenvalues will be found.\n*          = 'V': all eigenvalues in the half-open interval (VL,VU]\n*                 will be found.\n*          = 'I': the IL-th through IU-th eigenvalues will be found.\n*\n*  N       (input) INTEGER\n*          The order of the matrix.  N >= 0.\n*\n*  D       (input/output) REAL array, dimension (N)\n*          On entry, the n diagonal elements of the tridiagonal matrix\n*          A.\n*          On exit, D may be multiplied by a constant factor chosen\n*          to avoid over/underflow in computing the eigenvalues.\n*\n*  E       (input/output) REAL array, dimension (max(1,N-1))\n*          On entry, the (n-1) subdiagonal elements of the tridiagonal\n*          matrix A in elements 1 to N-1 of E.\n*          On exit, E may be multiplied by a constant factor chosen\n*          to avoid over/underflow in computing the eigenvalues.\n*\n*  VL      (input) REAL\n*  VU      (input) REAL\n*          If RANGE='V', the lower and upper bounds of the interval to\n*          be searched for eigenvalues. VL < VU.\n*          Not referenced if RANGE = 'A' or 'I'.\n*\n*  IL      (input) INTEGER\n*  IU      (input) INTEGER\n*          If RANGE='I', the indices (in ascending order) of the\n*          smallest and largest eigenvalues to be returned.\n*          1 <= IL <= IU <= N, if N > 0; IL = 1 and IU = 0 if N = 0.\n*          Not referenced if RANGE = 'A' or 'V'.\n*\n*  ABSTOL  (input) REAL\n*          The absolute error tolerance for the eigenvalues.\n*          An approximate eigenvalue is accepted as converged\n*          when it is determined to lie in an interval [a,b]\n*          of width less than or equal to\n*\n*                  ABSTOL + EPS *   max( |a|,|b| ) ,\n*\n*          where EPS is the machine precision.  If ABSTOL is less\n*          than or equal to zero, then  EPS*|T|  will be used in\n*          its place, where |T| is the 1-norm of the tridiagonal\n*          matrix.\n*\n*          Eigenvalues will be computed most accurately when ABSTOL is\n*          set to twice the underflow threshold 2*SLAMCH('S'), not zero.\n*          If this routine returns with INFO>0, indicating that some\n*          eigenvectors did not converge, try setting ABSTOL to\n*          2*SLAMCH('S').\n*\n*          See \"Computing Small Singular Values of Bidiagonal Matrices\n*          with Guaranteed High Relative Accuracy,\" by Demmel and\n*          Kahan, LAPACK Working Note #3.\n*\n*  M       (output) INTEGER\n*          The total number of eigenvalues found.  0 <= M <= N.\n*          If RANGE = 'A', M = N, and if RANGE = 'I', M = IU-IL+1.\n*\n*  W       (output) REAL array, dimension (N)\n*          The first M elements contain the selected eigenvalues in\n*          ascending order.\n*\n*  Z       (output) REAL array, dimension (LDZ, max(1,M) )\n*          If JOBZ = 'V', then if INFO = 0, the first M columns of Z\n*          contain the orthonormal eigenvectors of the matrix A\n*          corresponding to the selected eigenvalues, with the i-th\n*          column of Z holding the eigenvector associated with W(i).\n*          If an eigenvector fails to converge (INFO > 0), then that\n*          column of Z contains the latest approximation to the\n*          eigenvector, and the index of the eigenvector is returned\n*          in IFAIL.  If JOBZ = 'N', then Z is not referenced.\n*          Note: the user must ensure that at least max(1,M) columns are\n*          supplied in the array Z; if RANGE = 'V', the exact value of M\n*          is not known in advance and an upper bound must be used.\n*\n*  LDZ     (input) INTEGER\n*          The leading dimension of the array Z.  LDZ >= 1, and if\n*          JOBZ = 'V', LDZ >= max(1,N).\n*\n*  WORK    (workspace) REAL array, dimension (5*N)\n*\n*  IWORK   (workspace) INTEGER array, dimension (5*N)\n*\n*  IFAIL   (output) INTEGER array, dimension (N)\n*          If JOBZ = 'V', then if INFO = 0, the first M elements of\n*          IFAIL are zero.  If INFO > 0, then IFAIL contains the\n*          indices of the eigenvectors that failed to converge.\n*          If JOBZ = 'N', then IFAIL is not referenced.\n*\n*  INFO    (output) INTEGER\n*          = 0:  successful exit\n*          < 0:  if INFO = -i, the i-th argument had an illegal value\n*          > 0:  if INFO = i, then i eigenvectors failed to converge.\n*                Their indices are stored in array IFAIL.\n*\n\n*  =====================================================================\n*\n\n");
    return Qnil;
  }
  if (argc != 9)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 9)", argc);
  rb_jobz = argv[0];
  rb_range = argv[1];
  rb_d = argv[2];
  rb_e = argv[3];
  rb_vl = argv[4];
  rb_vu = argv[5];
  rb_il = argv[6];
  rb_iu = argv[7];
  rb_abstol = argv[8];

  abstol = (real)NUM2DBL(rb_abstol);
  vl = (real)NUM2DBL(rb_vl);
  iu = NUM2INT(rb_iu);
  jobz = StringValueCStr(rb_jobz)[0];
  vu = (real)NUM2DBL(rb_vu);
  range = StringValueCStr(rb_range)[0];
  if (!NA_IsNArray(rb_d))
    rb_raise(rb_eArgError, "d (3th argument) must be NArray");
  if (NA_RANK(rb_d) != 1)
    rb_raise(rb_eArgError, "rank of d (3th argument) must be %d", 1);
  n = NA_SHAPE0(rb_d);
  if (NA_TYPE(rb_d) != NA_SFLOAT)
    rb_d = na_change_type(rb_d, NA_SFLOAT);
  d = NA_PTR_TYPE(rb_d, real*);
  il = NUM2INT(rb_il);
  if (!NA_IsNArray(rb_e))
    rb_raise(rb_eArgError, "e (4th argument) must be NArray");
  if (NA_RANK(rb_e) != 1)
    rb_raise(rb_eArgError, "rank of e (4th argument) must be %d", 1);
  if (NA_SHAPE0(rb_e) != (MAX(1,n-1)))
    rb_raise(rb_eRuntimeError, "shape 0 of e must be %d", MAX(1,n-1));
  if (NA_TYPE(rb_e) != NA_SFLOAT)
    rb_e = na_change_type(rb_e, NA_SFLOAT);
  e = NA_PTR_TYPE(rb_e, real*);
  m = n;
  ldz = lsame_(&jobz,"V") ? MAX(1,n) : 1;
  {
    int shape[1];
    shape[0] = n;
    rb_w = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  w = NA_PTR_TYPE(rb_w, real*);
  {
    int shape[2];
    shape[0] = ldz;
    shape[1] = MAX(1,m);
    rb_z = na_make_object(NA_SFLOAT, 2, shape, cNArray);
  }
  z = NA_PTR_TYPE(rb_z, real*);
  {
    int shape[1];
    shape[0] = n;
    rb_ifail = na_make_object(NA_LINT, 1, shape, cNArray);
  }
  ifail = NA_PTR_TYPE(rb_ifail, integer*);
  {
    int shape[1];
    shape[0] = n;
    rb_d_out__ = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  d_out__ = NA_PTR_TYPE(rb_d_out__, real*);
  MEMCPY(d_out__, d, real, NA_TOTAL(rb_d));
  rb_d = rb_d_out__;
  d = d_out__;
  {
    int shape[1];
    shape[0] = MAX(1,n-1);
    rb_e_out__ = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  e_out__ = NA_PTR_TYPE(rb_e_out__, real*);
  MEMCPY(e_out__, e, real, NA_TOTAL(rb_e));
  rb_e = rb_e_out__;
  e = e_out__;
  work = ALLOC_N(real, (5*n));
  iwork = ALLOC_N(integer, (5*n));

  sstevx_(&jobz, &range, &n, d, e, &vl, &vu, &il, &iu, &abstol, &m, w, z, &ldz, work, iwork, ifail, &info);

  free(work);
  free(iwork);
  rb_m = INT2NUM(m);
  rb_info = INT2NUM(info);
  return rb_ary_new3(7, rb_m, rb_w, rb_z, rb_ifail, rb_info, rb_d, rb_e);
}

void
init_lapack_sstevx(VALUE mLapack){
  rb_define_module_function(mLapack, "sstevx", rb_sstevx, -1);
}
