#include "rb_lapack.h"

extern VOID sstein_(integer *n, real *d, real *e, integer *m, real *w, integer *iblock, integer *isplit, real *z, integer *ldz, real *work, integer *iwork, integer *ifail, integer *info);

static VALUE sHelp, sUsage;

static VALUE
rblapack_sstein(int argc, VALUE *argv, VALUE self){
  VALUE rblapack_d;
  real *d; 
  VALUE rblapack_e;
  real *e; 
  VALUE rblapack_w;
  real *w; 
  VALUE rblapack_iblock;
  integer *iblock; 
  VALUE rblapack_isplit;
  integer *isplit; 
  VALUE rblapack_z;
  real *z; 
  VALUE rblapack_ifail;
  integer *ifail; 
  VALUE rblapack_info;
  integer info; 
  real *work;
  integer *iwork;

  integer n;
  integer ldz;
  integer m;

  VALUE rb_options;
  if (argc > 0 && TYPE(argv[argc-1]) == T_HASH) {
    argc--;
    rb_options = argv[argc];
    if (rb_hash_aref(rb_options, sHelp) == Qtrue) {
      printf("%s\n", "USAGE:\n  z, ifail, info = NumRu::Lapack.sstein( d, e, w, iblock, isplit, [:usage => usage, :help => help])\n\n\nFORTRAN MANUAL\n      SUBROUTINE SSTEIN( N, D, E, M, W, IBLOCK, ISPLIT, Z, LDZ, WORK, IWORK, IFAIL, INFO )\n\n*  Purpose\n*  =======\n*\n*  SSTEIN computes the eigenvectors of a real symmetric tridiagonal\n*  matrix T corresponding to specified eigenvalues, using inverse\n*  iteration.\n*\n*  The maximum number of iterations allowed for each eigenvector is\n*  specified by an internal parameter MAXITS (currently set to 5).\n*\n\n*  Arguments\n*  =========\n*\n*  N       (input) INTEGER\n*          The order of the matrix.  N >= 0.\n*\n*  D       (input) REAL array, dimension (N)\n*          The n diagonal elements of the tridiagonal matrix T.\n*\n*  E       (input) REAL array, dimension (N-1)\n*          The (n-1) subdiagonal elements of the tridiagonal matrix\n*          T, in elements 1 to N-1.\n*\n*  M       (input) INTEGER\n*          The number of eigenvectors to be found.  0 <= M <= N.\n*\n*  W       (input) REAL array, dimension (N)\n*          The first M elements of W contain the eigenvalues for\n*          which eigenvectors are to be computed.  The eigenvalues\n*          should be grouped by split-off block and ordered from\n*          smallest to largest within the block.  ( The output array\n*          W from SSTEBZ with ORDER = 'B' is expected here. )\n*\n*  IBLOCK  (input) INTEGER array, dimension (N)\n*          The submatrix indices associated with the corresponding\n*          eigenvalues in W; IBLOCK(i)=1 if eigenvalue W(i) belongs to\n*          the first submatrix from the top, =2 if W(i) belongs to\n*          the second submatrix, etc.  ( The output array IBLOCK\n*          from SSTEBZ is expected here. )\n*\n*  ISPLIT  (input) INTEGER array, dimension (N)\n*          The splitting points, at which T breaks up into submatrices.\n*          The first submatrix consists of rows/columns 1 to\n*          ISPLIT( 1 ), the second of rows/columns ISPLIT( 1 )+1\n*          through ISPLIT( 2 ), etc.\n*          ( The output array ISPLIT from SSTEBZ is expected here. )\n*\n*  Z       (output) REAL array, dimension (LDZ, M)\n*          The computed eigenvectors.  The eigenvector associated\n*          with the eigenvalue W(i) is stored in the i-th column of\n*          Z.  Any vector which fails to converge is set to its current\n*          iterate after MAXITS iterations.\n*\n*  LDZ     (input) INTEGER\n*          The leading dimension of the array Z.  LDZ >= max(1,N).\n*\n*  WORK    (workspace) REAL array, dimension (5*N)\n*\n*  IWORK   (workspace) INTEGER array, dimension (N)\n*\n*  IFAIL   (output) INTEGER array, dimension (M)\n*          On normal exit, all elements of IFAIL are zero.\n*          If one or more eigenvectors fail to converge after\n*          MAXITS iterations, then their indices are stored in\n*          array IFAIL.\n*\n*  INFO    (output) INTEGER\n*          = 0: successful exit.\n*          < 0: if INFO = -i, the i-th argument had an illegal value\n*          > 0: if INFO = i, then i eigenvectors failed to converge\n*               in MAXITS iterations.  Their indices are stored in\n*               array IFAIL.\n*\n*  Internal Parameters\n*  ===================\n*\n*  MAXITS  INTEGER, default = 5\n*          The maximum number of iterations performed.\n*\n*  EXTRA   INTEGER, default = 2\n*          The number of iterations performed after norm growth\n*          criterion is satisfied, should be at least 1.\n*\n\n*  =====================================================================\n*\n\n");
      rb_exit(0);
    }
    if (rb_hash_aref(rb_options, sUsage) == Qtrue) {
      printf("%s\n", "USAGE:\n  z, ifail, info = NumRu::Lapack.sstein( d, e, w, iblock, isplit, [:usage => usage, :help => help])\n");
      rb_exit(0);
    } 
  } else
    rb_options = Qnil;
  if (argc != 5)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 5)", argc);
  rblapack_d = argv[0];
  rblapack_e = argv[1];
  rblapack_w = argv[2];
  rblapack_iblock = argv[3];
  rblapack_isplit = argv[4];
  if (rb_options != Qnil) {
  }

  if (!NA_IsNArray(rblapack_w))
    rb_raise(rb_eArgError, "w (3th argument) must be NArray");
  if (NA_RANK(rblapack_w) != 1)
    rb_raise(rb_eArgError, "rank of w (3th argument) must be %d", 1);
  n = NA_SHAPE0(rblapack_w);
  if (NA_TYPE(rblapack_w) != NA_SFLOAT)
    rblapack_w = na_change_type(rblapack_w, NA_SFLOAT);
  w = NA_PTR_TYPE(rblapack_w, real*);
  if (!NA_IsNArray(rblapack_iblock))
    rb_raise(rb_eArgError, "iblock (4th argument) must be NArray");
  if (NA_RANK(rblapack_iblock) != 1)
    rb_raise(rb_eArgError, "rank of iblock (4th argument) must be %d", 1);
  if (NA_SHAPE0(rblapack_iblock) != n)
    rb_raise(rb_eRuntimeError, "shape 0 of iblock must be the same as shape 0 of w");
  if (NA_TYPE(rblapack_iblock) != NA_LINT)
    rblapack_iblock = na_change_type(rblapack_iblock, NA_LINT);
  iblock = NA_PTR_TYPE(rblapack_iblock, integer*);
  if (!NA_IsNArray(rblapack_isplit))
    rb_raise(rb_eArgError, "isplit (5th argument) must be NArray");
  if (NA_RANK(rblapack_isplit) != 1)
    rb_raise(rb_eArgError, "rank of isplit (5th argument) must be %d", 1);
  if (NA_SHAPE0(rblapack_isplit) != n)
    rb_raise(rb_eRuntimeError, "shape 0 of isplit must be the same as shape 0 of w");
  if (NA_TYPE(rblapack_isplit) != NA_LINT)
    rblapack_isplit = na_change_type(rblapack_isplit, NA_LINT);
  isplit = NA_PTR_TYPE(rblapack_isplit, integer*);
  if (!NA_IsNArray(rblapack_d))
    rb_raise(rb_eArgError, "d (1th argument) must be NArray");
  if (NA_RANK(rblapack_d) != 1)
    rb_raise(rb_eArgError, "rank of d (1th argument) must be %d", 1);
  if (NA_SHAPE0(rblapack_d) != n)
    rb_raise(rb_eRuntimeError, "shape 0 of d must be the same as shape 0 of w");
  if (NA_TYPE(rblapack_d) != NA_SFLOAT)
    rblapack_d = na_change_type(rblapack_d, NA_SFLOAT);
  d = NA_PTR_TYPE(rblapack_d, real*);
  if (!NA_IsNArray(rblapack_e))
    rb_raise(rb_eArgError, "e (2th argument) must be NArray");
  if (NA_RANK(rblapack_e) != 1)
    rb_raise(rb_eArgError, "rank of e (2th argument) must be %d", 1);
  if (NA_SHAPE0(rblapack_e) != (n-1))
    rb_raise(rb_eRuntimeError, "shape 0 of e must be %d", n-1);
  if (NA_TYPE(rblapack_e) != NA_SFLOAT)
    rblapack_e = na_change_type(rblapack_e, NA_SFLOAT);
  e = NA_PTR_TYPE(rblapack_e, real*);
  m = n;
  ldz = MAX(1,n);
  {
    int shape[2];
    shape[0] = ldz;
    shape[1] = m;
    rblapack_z = na_make_object(NA_SFLOAT, 2, shape, cNArray);
  }
  z = NA_PTR_TYPE(rblapack_z, real*);
  {
    int shape[1];
    shape[0] = m;
    rblapack_ifail = na_make_object(NA_LINT, 1, shape, cNArray);
  }
  ifail = NA_PTR_TYPE(rblapack_ifail, integer*);
  work = ALLOC_N(real, (5*n));
  iwork = ALLOC_N(integer, (n));

  sstein_(&n, d, e, &m, w, iblock, isplit, z, &ldz, work, iwork, ifail, &info);

  free(work);
  free(iwork);
  rblapack_info = INT2NUM(info);
  return rb_ary_new3(3, rblapack_z, rblapack_ifail, rblapack_info);
}

void
init_lapack_sstein(VALUE mLapack){
  rb_define_module_function(mLapack, "sstein", rblapack_sstein, -1);
  sHelp = ID2SYM(rb_intern("help"));
  sUsage = ID2SYM(rb_intern("usage"));
}
