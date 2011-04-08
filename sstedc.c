#include "rb_lapack.h"

extern VOID sstedc_(char *compz, integer *n, real *d, real *e, real *z, integer *ldz, real *work, integer *lwork, integer *iwork, integer *liwork, integer *info);

static VALUE sHelp, sUsage;

static VALUE
rblapack_sstedc(int argc, VALUE *argv, VALUE self){
  VALUE rblapack_compz;
  char compz; 
  VALUE rblapack_d;
  real *d; 
  VALUE rblapack_e;
  real *e; 
  VALUE rblapack_z;
  real *z; 
  VALUE rblapack_lwork;
  integer lwork; 
  VALUE rblapack_liwork;
  integer liwork; 
  VALUE rblapack_work;
  real *work; 
  VALUE rblapack_iwork;
  integer *iwork; 
  VALUE rblapack_info;
  integer info; 
  VALUE rblapack_d_out__;
  real *d_out__;
  VALUE rblapack_e_out__;
  real *e_out__;
  VALUE rblapack_z_out__;
  real *z_out__;

  integer n;
  integer ldz;

  VALUE rb_options;
  if (argc > 0 && TYPE(argv[argc-1]) == T_HASH) {
    argc--;
    rb_options = argv[argc];
    if (rb_hash_aref(rb_options, sHelp) == Qtrue) {
      printf("%s\n", "USAGE:\n  work, iwork, info, d, e, z = NumRu::Lapack.sstedc( compz, d, e, z, lwork, liwork, [:usage => usage, :help => help])\n\n\nFORTRAN MANUAL\n      SUBROUTINE SSTEDC( COMPZ, N, D, E, Z, LDZ, WORK, LWORK, IWORK, LIWORK, INFO )\n\n*  Purpose\n*  =======\n*\n*  SSTEDC computes all eigenvalues and, optionally, eigenvectors of a\n*  symmetric tridiagonal matrix using the divide and conquer method.\n*  The eigenvectors of a full or band real symmetric matrix can also be\n*  found if SSYTRD or SSPTRD or SSBTRD has been used to reduce this\n*  matrix to tridiagonal form.\n*\n*  This code makes very mild assumptions about floating point\n*  arithmetic. It will work on machines with a guard digit in\n*  add/subtract, or on those binary machines without guard digits\n*  which subtract like the Cray X-MP, Cray Y-MP, Cray C-90, or Cray-2.\n*  It could conceivably fail on hexadecimal or decimal machines\n*  without guard digits, but we know of none.  See SLAED3 for details.\n*\n\n*  Arguments\n*  =========\n*\n*  COMPZ   (input) CHARACTER*1\n*          = 'N':  Compute eigenvalues only.\n*          = 'I':  Compute eigenvectors of tridiagonal matrix also.\n*          = 'V':  Compute eigenvectors of original dense symmetric\n*                  matrix also.  On entry, Z contains the orthogonal\n*                  matrix used to reduce the original matrix to\n*                  tridiagonal form.\n*\n*  N       (input) INTEGER\n*          The dimension of the symmetric tridiagonal matrix.  N >= 0.\n*\n*  D       (input/output) REAL array, dimension (N)\n*          On entry, the diagonal elements of the tridiagonal matrix.\n*          On exit, if INFO = 0, the eigenvalues in ascending order.\n*\n*  E       (input/output) REAL array, dimension (N-1)\n*          On entry, the subdiagonal elements of the tridiagonal matrix.\n*          On exit, E has been destroyed.\n*\n*  Z       (input/output) REAL array, dimension (LDZ,N)\n*          On entry, if COMPZ = 'V', then Z contains the orthogonal\n*          matrix used in the reduction to tridiagonal form.\n*          On exit, if INFO = 0, then if COMPZ = 'V', Z contains the\n*          orthonormal eigenvectors of the original symmetric matrix,\n*          and if COMPZ = 'I', Z contains the orthonormal eigenvectors\n*          of the symmetric tridiagonal matrix.\n*          If  COMPZ = 'N', then Z is not referenced.\n*\n*  LDZ     (input) INTEGER\n*          The leading dimension of the array Z.  LDZ >= 1.\n*          If eigenvectors are desired, then LDZ >= max(1,N).\n*\n*  WORK    (workspace/output) REAL array, dimension (MAX(1,LWORK))\n*          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.\n*\n*  LWORK   (input) INTEGER\n*          The dimension of the array WORK.\n*          If COMPZ = 'N' or N <= 1 then LWORK must be at least 1.\n*          If COMPZ = 'V' and N > 1 then LWORK must be at least\n*                         ( 1 + 3*N + 2*N*lg N + 3*N**2 ),\n*                         where lg( N ) = smallest integer k such\n*                         that 2**k >= N.\n*          If COMPZ = 'I' and N > 1 then LWORK must be at least\n*                         ( 1 + 4*N + N**2 ).\n*          Note that for COMPZ = 'I' or 'V', then if N is less than or\n*          equal to the minimum divide size, usually 25, then LWORK need\n*          only be max(1,2*(N-1)).\n*\n*          If LWORK = -1, then a workspace query is assumed; the routine\n*          only calculates the optimal size of the WORK array, returns\n*          this value as the first entry of the WORK array, and no error\n*          message related to LWORK is issued by XERBLA.\n*\n*  IWORK   (workspace/output) INTEGER array, dimension (MAX(1,LIWORK))\n*          On exit, if INFO = 0, IWORK(1) returns the optimal LIWORK.\n*\n*  LIWORK  (input) INTEGER\n*          The dimension of the array IWORK.\n*          If COMPZ = 'N' or N <= 1 then LIWORK must be at least 1.\n*          If COMPZ = 'V' and N > 1 then LIWORK must be at least\n*                         ( 6 + 6*N + 5*N*lg N ).\n*          If COMPZ = 'I' and N > 1 then LIWORK must be at least\n*                         ( 3 + 5*N ).\n*          Note that for COMPZ = 'I' or 'V', then if N is less than or\n*          equal to the minimum divide size, usually 25, then LIWORK\n*          need only be 1.\n*\n*          If LIWORK = -1, then a workspace query is assumed; the\n*          routine only calculates the optimal size of the IWORK array,\n*          returns this value as the first entry of the IWORK array, and\n*          no error message related to LIWORK is issued by XERBLA.\n*\n*  INFO    (output) INTEGER\n*          = 0:  successful exit.\n*          < 0:  if INFO = -i, the i-th argument had an illegal value.\n*          > 0:  The algorithm failed to compute an eigenvalue while\n*                working on the submatrix lying in rows and columns\n*                INFO/(N+1) through mod(INFO,N+1).\n*\n\n*  Further Details\n*  ===============\n*\n*  Based on contributions by\n*     Jeff Rutter, Computer Science Division, University of California\n*     at Berkeley, USA\n*  Modified by Francoise Tisseur, University of Tennessee.\n*\n*  =====================================================================\n*\n\n");
      rb_exit(0);
    }
    if (rb_hash_aref(rb_options, sUsage) == Qtrue) {
      printf("%s\n", "USAGE:\n  work, iwork, info, d, e, z = NumRu::Lapack.sstedc( compz, d, e, z, lwork, liwork, [:usage => usage, :help => help])\n");
      rb_exit(0);
    } 
  } else
    rb_options = Qnil;
  if (argc != 6)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 6)", argc);
  rblapack_compz = argv[0];
  rblapack_d = argv[1];
  rblapack_e = argv[2];
  rblapack_z = argv[3];
  rblapack_lwork = argv[4];
  rblapack_liwork = argv[5];
  if (rb_options != Qnil) {
  }

  compz = StringValueCStr(rblapack_compz)[0];
  liwork = NUM2INT(rblapack_liwork);
  if (!NA_IsNArray(rblapack_z))
    rb_raise(rb_eArgError, "z (4th argument) must be NArray");
  if (NA_RANK(rblapack_z) != 2)
    rb_raise(rb_eArgError, "rank of z (4th argument) must be %d", 2);
  n = NA_SHAPE1(rblapack_z);
  ldz = NA_SHAPE0(rblapack_z);
  if (NA_TYPE(rblapack_z) != NA_SFLOAT)
    rblapack_z = na_change_type(rblapack_z, NA_SFLOAT);
  z = NA_PTR_TYPE(rblapack_z, real*);
  lwork = NUM2INT(rblapack_lwork);
  if (!NA_IsNArray(rblapack_d))
    rb_raise(rb_eArgError, "d (2th argument) must be NArray");
  if (NA_RANK(rblapack_d) != 1)
    rb_raise(rb_eArgError, "rank of d (2th argument) must be %d", 1);
  if (NA_SHAPE0(rblapack_d) != n)
    rb_raise(rb_eRuntimeError, "shape 0 of d must be the same as shape 1 of z");
  if (NA_TYPE(rblapack_d) != NA_SFLOAT)
    rblapack_d = na_change_type(rblapack_d, NA_SFLOAT);
  d = NA_PTR_TYPE(rblapack_d, real*);
  if (!NA_IsNArray(rblapack_e))
    rb_raise(rb_eArgError, "e (3th argument) must be NArray");
  if (NA_RANK(rblapack_e) != 1)
    rb_raise(rb_eArgError, "rank of e (3th argument) must be %d", 1);
  if (NA_SHAPE0(rblapack_e) != (n-1))
    rb_raise(rb_eRuntimeError, "shape 0 of e must be %d", n-1);
  if (NA_TYPE(rblapack_e) != NA_SFLOAT)
    rblapack_e = na_change_type(rblapack_e, NA_SFLOAT);
  e = NA_PTR_TYPE(rblapack_e, real*);
  {
    int shape[1];
    shape[0] = MAX(1,lwork);
    rblapack_work = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  work = NA_PTR_TYPE(rblapack_work, real*);
  {
    int shape[1];
    shape[0] = MAX(1,liwork);
    rblapack_iwork = na_make_object(NA_LINT, 1, shape, cNArray);
  }
  iwork = NA_PTR_TYPE(rblapack_iwork, integer*);
  {
    int shape[1];
    shape[0] = n;
    rblapack_d_out__ = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  d_out__ = NA_PTR_TYPE(rblapack_d_out__, real*);
  MEMCPY(d_out__, d, real, NA_TOTAL(rblapack_d));
  rblapack_d = rblapack_d_out__;
  d = d_out__;
  {
    int shape[1];
    shape[0] = n-1;
    rblapack_e_out__ = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  e_out__ = NA_PTR_TYPE(rblapack_e_out__, real*);
  MEMCPY(e_out__, e, real, NA_TOTAL(rblapack_e));
  rblapack_e = rblapack_e_out__;
  e = e_out__;
  {
    int shape[2];
    shape[0] = ldz;
    shape[1] = n;
    rblapack_z_out__ = na_make_object(NA_SFLOAT, 2, shape, cNArray);
  }
  z_out__ = NA_PTR_TYPE(rblapack_z_out__, real*);
  MEMCPY(z_out__, z, real, NA_TOTAL(rblapack_z));
  rblapack_z = rblapack_z_out__;
  z = z_out__;

  sstedc_(&compz, &n, d, e, z, &ldz, work, &lwork, iwork, &liwork, &info);

  rblapack_info = INT2NUM(info);
  return rb_ary_new3(6, rblapack_work, rblapack_iwork, rblapack_info, rblapack_d, rblapack_e, rblapack_z);
}

void
init_lapack_sstedc(VALUE mLapack){
  rb_define_module_function(mLapack, "sstedc", rblapack_sstedc, -1);
  sHelp = ID2SYM(rb_intern("help"));
  sUsage = ID2SYM(rb_intern("usage"));
}
