#include "rb_lapack.h"

extern VOID dstevr_(char *jobz, char *range, integer *n, doublereal *d, doublereal *e, doublereal *vl, doublereal *vu, integer *il, integer *iu, doublereal *abstol, integer *m, doublereal *w, doublereal *z, integer *ldz, integer *isuppz, doublereal *work, integer *lwork, integer *iwork, integer *liwork, integer *info);

static VALUE sHelp, sUsage;

static VALUE
rblapack_dstevr(int argc, VALUE *argv, VALUE self){
  VALUE rblapack_jobz;
  char jobz; 
  VALUE rblapack_range;
  char range; 
  VALUE rblapack_d;
  doublereal *d; 
  VALUE rblapack_e;
  doublereal *e; 
  VALUE rblapack_vl;
  doublereal vl; 
  VALUE rblapack_vu;
  doublereal vu; 
  VALUE rblapack_il;
  integer il; 
  VALUE rblapack_iu;
  integer iu; 
  VALUE rblapack_abstol;
  doublereal abstol; 
  VALUE rblapack_lwork;
  integer lwork; 
  VALUE rblapack_liwork;
  integer liwork; 
  VALUE rblapack_m;
  integer m; 
  VALUE rblapack_w;
  doublereal *w; 
  VALUE rblapack_z;
  doublereal *z; 
  VALUE rblapack_isuppz;
  integer *isuppz; 
  VALUE rblapack_work;
  doublereal *work; 
  VALUE rblapack_iwork;
  integer *iwork; 
  VALUE rblapack_info;
  integer info; 
  VALUE rblapack_d_out__;
  doublereal *d_out__;
  VALUE rblapack_e_out__;
  doublereal *e_out__;

  integer n;
  integer ldz;

  VALUE rb_options;
  if (argc > 0 && TYPE(argv[argc-1]) == T_HASH) {
    argc--;
    rb_options = argv[argc];
    if (rb_hash_aref(rb_options, sHelp) == Qtrue) {
      printf("%s\n", "USAGE:\n  m, w, z, isuppz, work, iwork, info, d, e = NumRu::Lapack.dstevr( jobz, range, d, e, vl, vu, il, iu, abstol, lwork, liwork, [:usage => usage, :help => help])\n\n\nFORTRAN MANUAL\n      SUBROUTINE DSTEVR( JOBZ, RANGE, N, D, E, VL, VU, IL, IU, ABSTOL, M, W, Z, LDZ, ISUPPZ, WORK, LWORK, IWORK, LIWORK, INFO )\n\n*  Purpose\n*  =======\n*\n*  DSTEVR computes selected eigenvalues and, optionally, eigenvectors\n*  of a real symmetric tridiagonal matrix T.  Eigenvalues and\n*  eigenvectors can be selected by specifying either a range of values\n*  or a range of indices for the desired eigenvalues.\n*\n*  Whenever possible, DSTEVR calls DSTEMR to compute the\n*  eigenspectrum using Relatively Robust Representations.  DSTEMR\n*  computes eigenvalues by the dqds algorithm, while orthogonal\n*  eigenvectors are computed from various \"good\" L D L^T representations\n*  (also known as Relatively Robust Representations). Gram-Schmidt\n*  orthogonalization is avoided as far as possible. More specifically,\n*  the various steps of the algorithm are as follows. For the i-th\n*  unreduced block of T,\n*     (a) Compute T - sigma_i = L_i D_i L_i^T, such that L_i D_i L_i^T\n*          is a relatively robust representation,\n*     (b) Compute the eigenvalues, lambda_j, of L_i D_i L_i^T to high\n*         relative accuracy by the dqds algorithm,\n*     (c) If there is a cluster of close eigenvalues, \"choose\" sigma_i\n*         close to the cluster, and go to step (a),\n*     (d) Given the approximate eigenvalue lambda_j of L_i D_i L_i^T,\n*         compute the corresponding eigenvector by forming a\n*         rank-revealing twisted factorization.\n*  The desired accuracy of the output can be specified by the input\n*  parameter ABSTOL.\n*\n*  For more details, see \"A new O(n^2) algorithm for the symmetric\n*  tridiagonal eigenvalue/eigenvector problem\", by Inderjit Dhillon,\n*  Computer Science Division Technical Report No. UCB//CSD-97-971,\n*  UC Berkeley, May 1997.\n*\n*\n*  Note 1 : DSTEVR calls DSTEMR when the full spectrum is requested\n*  on machines which conform to the ieee-754 floating point standard.\n*  DSTEVR calls DSTEBZ and DSTEIN on non-ieee machines and\n*  when partial spectrum requests are made.\n*\n*  Normal execution of DSTEMR may create NaNs and infinities and\n*  hence may abort due to a floating point exception in environments\n*  which do not handle NaNs and infinities in the ieee standard default\n*  manner.\n*\n\n*  Arguments\n*  =========\n*\n*  JOBZ    (input) CHARACTER*1\n*          = 'N':  Compute eigenvalues only;\n*          = 'V':  Compute eigenvalues and eigenvectors.\n*\n*  RANGE   (input) CHARACTER*1\n*          = 'A': all eigenvalues will be found.\n*          = 'V': all eigenvalues in the half-open interval (VL,VU]\n*                 will be found.\n*          = 'I': the IL-th through IU-th eigenvalues will be found.\n********** For RANGE = 'V' or 'I' and IU - IL < N - 1, DSTEBZ and\n********** DSTEIN are called\n*\n*  N       (input) INTEGER\n*          The order of the matrix.  N >= 0.\n*\n*  D       (input/output) DOUBLE PRECISION array, dimension (N)\n*          On entry, the n diagonal elements of the tridiagonal matrix\n*          A.\n*          On exit, D may be multiplied by a constant factor chosen\n*          to avoid over/underflow in computing the eigenvalues.\n*\n*  E       (input/output) DOUBLE PRECISION array, dimension (max(1,N-1))\n*          On entry, the (n-1) subdiagonal elements of the tridiagonal\n*          matrix A in elements 1 to N-1 of E.\n*          On exit, E may be multiplied by a constant factor chosen\n*          to avoid over/underflow in computing the eigenvalues.\n*\n*  VL      (input) DOUBLE PRECISION\n*  VU      (input) DOUBLE PRECISION\n*          If RANGE='V', the lower and upper bounds of the interval to\n*          be searched for eigenvalues. VL < VU.\n*          Not referenced if RANGE = 'A' or 'I'.\n*\n*  IL      (input) INTEGER\n*  IU      (input) INTEGER\n*          If RANGE='I', the indices (in ascending order) of the\n*          smallest and largest eigenvalues to be returned.\n*          1 <= IL <= IU <= N, if N > 0; IL = 1 and IU = 0 if N = 0.\n*          Not referenced if RANGE = 'A' or 'V'.\n*\n*  ABSTOL  (input) DOUBLE PRECISION\n*          The absolute error tolerance for the eigenvalues.\n*          An approximate eigenvalue is accepted as converged\n*          when it is determined to lie in an interval [a,b]\n*          of width less than or equal to\n*\n*                  ABSTOL + EPS *   max( |a|,|b| ) ,\n*\n*          where EPS is the machine precision.  If ABSTOL is less than\n*          or equal to zero, then  EPS*|T|  will be used in its place,\n*          where |T| is the 1-norm of the tridiagonal matrix obtained\n*          by reducing A to tridiagonal form.\n*\n*          See \"Computing Small Singular Values of Bidiagonal Matrices\n*          with Guaranteed High Relative Accuracy,\" by Demmel and\n*          Kahan, LAPACK Working Note #3.\n*\n*          If high relative accuracy is important, set ABSTOL to\n*          DLAMCH( 'Safe minimum' ).  Doing so will guarantee that\n*          eigenvalues are computed to high relative accuracy when\n*          possible in future releases.  The current code does not\n*          make any guarantees about high relative accuracy, but\n*          future releases will. See J. Barlow and J. Demmel,\n*          \"Computing Accurate Eigensystems of Scaled Diagonally\n*          Dominant Matrices\", LAPACK Working Note #7, for a discussion\n*          of which matrices define their eigenvalues to high relative\n*          accuracy.\n*\n*  M       (output) INTEGER\n*          The total number of eigenvalues found.  0 <= M <= N.\n*          If RANGE = 'A', M = N, and if RANGE = 'I', M = IU-IL+1.\n*\n*  W       (output) DOUBLE PRECISION array, dimension (N)\n*          The first M elements contain the selected eigenvalues in\n*          ascending order.\n*\n*  Z       (output) DOUBLE PRECISION array, dimension (LDZ, max(1,M) )\n*          If JOBZ = 'V', then if INFO = 0, the first M columns of Z\n*          contain the orthonormal eigenvectors of the matrix A\n*          corresponding to the selected eigenvalues, with the i-th\n*          column of Z holding the eigenvector associated with W(i).\n*          Note: the user must ensure that at least max(1,M) columns are\n*          supplied in the array Z; if RANGE = 'V', the exact value of M\n*          is not known in advance and an upper bound must be used.\n*\n*  LDZ     (input) INTEGER\n*          The leading dimension of the array Z.  LDZ >= 1, and if\n*          JOBZ = 'V', LDZ >= max(1,N).\n*\n*  ISUPPZ  (output) INTEGER array, dimension ( 2*max(1,M) )\n*          The support of the eigenvectors in Z, i.e., the indices\n*          indicating the nonzero elements in Z. The i-th eigenvector\n*          is nonzero only in elements ISUPPZ( 2*i-1 ) through\n*          ISUPPZ( 2*i ).\n********** Implemented only for RANGE = 'A' or 'I' and IU - IL = N - 1\n*\n*  WORK    (workspace/output) DOUBLE PRECISION array, dimension (MAX(1,LWORK))\n*          On exit, if INFO = 0, WORK(1) returns the optimal (and\n*          minimal) LWORK.\n*\n*  LWORK   (input) INTEGER\n*          The dimension of the array WORK.  LWORK >= max(1,20*N).\n*\n*          If LWORK = -1, then a workspace query is assumed; the routine\n*          only calculates the optimal sizes of the WORK and IWORK\n*          arrays, returns these values as the first entries of the WORK\n*          and IWORK arrays, and no error message related to LWORK or\n*          LIWORK is issued by XERBLA.\n*\n*  IWORK   (workspace/output) INTEGER array, dimension (MAX(1,LIWORK))\n*          On exit, if INFO = 0, IWORK(1) returns the optimal (and\n*          minimal) LIWORK.\n*\n*  LIWORK  (input) INTEGER\n*          The dimension of the array IWORK.  LIWORK >= max(1,10*N).\n*\n*          If LIWORK = -1, then a workspace query is assumed; the\n*          routine only calculates the optimal sizes of the WORK and\n*          IWORK arrays, returns these values as the first entries of\n*          the WORK and IWORK arrays, and no error message related to\n*          LWORK or LIWORK is issued by XERBLA.\n*\n*  INFO    (output) INTEGER\n*          = 0:  successful exit\n*          < 0:  if INFO = -i, the i-th argument had an illegal value\n*          > 0:  Internal error\n*\n\n*  Further Details\n*  ===============\n*\n*  Based on contributions by\n*     Inderjit Dhillon, IBM Almaden, USA\n*     Osni Marques, LBNL/NERSC, USA\n*     Ken Stanley, Computer Science Division, University of\n*       California at Berkeley, USA\n*\n*  =====================================================================\n*\n\n");
      rb_exit(0);
    }
    if (rb_hash_aref(rb_options, sUsage) == Qtrue) {
      printf("%s\n", "USAGE:\n  m, w, z, isuppz, work, iwork, info, d, e = NumRu::Lapack.dstevr( jobz, range, d, e, vl, vu, il, iu, abstol, lwork, liwork, [:usage => usage, :help => help])\n");
      rb_exit(0);
    } 
  } else
    rb_options = Qnil;
  if (argc != 11)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 11)", argc);
  rblapack_jobz = argv[0];
  rblapack_range = argv[1];
  rblapack_d = argv[2];
  rblapack_e = argv[3];
  rblapack_vl = argv[4];
  rblapack_vu = argv[5];
  rblapack_il = argv[6];
  rblapack_iu = argv[7];
  rblapack_abstol = argv[8];
  rblapack_lwork = argv[9];
  rblapack_liwork = argv[10];
  if (rb_options != Qnil) {
  }

  abstol = NUM2DBL(rblapack_abstol);
  vl = NUM2DBL(rblapack_vl);
  iu = NUM2INT(rblapack_iu);
  jobz = StringValueCStr(rblapack_jobz)[0];
  vu = NUM2DBL(rblapack_vu);
  liwork = NUM2INT(rblapack_liwork);
  range = StringValueCStr(rblapack_range)[0];
  if (!NA_IsNArray(rblapack_d))
    rb_raise(rb_eArgError, "d (3th argument) must be NArray");
  if (NA_RANK(rblapack_d) != 1)
    rb_raise(rb_eArgError, "rank of d (3th argument) must be %d", 1);
  n = NA_SHAPE0(rblapack_d);
  if (NA_TYPE(rblapack_d) != NA_DFLOAT)
    rblapack_d = na_change_type(rblapack_d, NA_DFLOAT);
  d = NA_PTR_TYPE(rblapack_d, doublereal*);
  lwork = NUM2INT(rblapack_lwork);
  il = NUM2INT(rblapack_il);
  if (!NA_IsNArray(rblapack_e))
    rb_raise(rb_eArgError, "e (4th argument) must be NArray");
  if (NA_RANK(rblapack_e) != 1)
    rb_raise(rb_eArgError, "rank of e (4th argument) must be %d", 1);
  if (NA_SHAPE0(rblapack_e) != (MAX(1,n-1)))
    rb_raise(rb_eRuntimeError, "shape 0 of e must be %d", MAX(1,n-1));
  if (NA_TYPE(rblapack_e) != NA_DFLOAT)
    rblapack_e = na_change_type(rblapack_e, NA_DFLOAT);
  e = NA_PTR_TYPE(rblapack_e, doublereal*);
  m = lsame_(&range,"I") ? iu-il+1 : n;
  ldz = lsame_(&jobz,"V") ? MAX(1,n) : 1;
  {
    int shape[1];
    shape[0] = n;
    rblapack_w = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  w = NA_PTR_TYPE(rblapack_w, doublereal*);
  {
    int shape[2];
    shape[0] = ldz;
    shape[1] = MAX(1,m);
    rblapack_z = na_make_object(NA_DFLOAT, 2, shape, cNArray);
  }
  z = NA_PTR_TYPE(rblapack_z, doublereal*);
  {
    int shape[1];
    shape[0] = 2*MAX(1,m);
    rblapack_isuppz = na_make_object(NA_LINT, 1, shape, cNArray);
  }
  isuppz = NA_PTR_TYPE(rblapack_isuppz, integer*);
  {
    int shape[1];
    shape[0] = MAX(1,lwork);
    rblapack_work = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  work = NA_PTR_TYPE(rblapack_work, doublereal*);
  {
    int shape[1];
    shape[0] = MAX(1,liwork);
    rblapack_iwork = na_make_object(NA_LINT, 1, shape, cNArray);
  }
  iwork = NA_PTR_TYPE(rblapack_iwork, integer*);
  {
    int shape[1];
    shape[0] = n;
    rblapack_d_out__ = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  d_out__ = NA_PTR_TYPE(rblapack_d_out__, doublereal*);
  MEMCPY(d_out__, d, doublereal, NA_TOTAL(rblapack_d));
  rblapack_d = rblapack_d_out__;
  d = d_out__;
  {
    int shape[1];
    shape[0] = MAX(1,n-1);
    rblapack_e_out__ = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  e_out__ = NA_PTR_TYPE(rblapack_e_out__, doublereal*);
  MEMCPY(e_out__, e, doublereal, NA_TOTAL(rblapack_e));
  rblapack_e = rblapack_e_out__;
  e = e_out__;

  dstevr_(&jobz, &range, &n, d, e, &vl, &vu, &il, &iu, &abstol, &m, w, z, &ldz, isuppz, work, &lwork, iwork, &liwork, &info);

  rblapack_m = INT2NUM(m);
  rblapack_info = INT2NUM(info);
  return rb_ary_new3(9, rblapack_m, rblapack_w, rblapack_z, rblapack_isuppz, rblapack_work, rblapack_iwork, rblapack_info, rblapack_d, rblapack_e);
}

void
init_lapack_dstevr(VALUE mLapack){
  rb_define_module_function(mLapack, "dstevr", rblapack_dstevr, -1);
  sHelp = ID2SYM(rb_intern("help"));
  sUsage = ID2SYM(rb_intern("usage"));
}
