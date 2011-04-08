#include "rb_lapack.h"

extern VOID slasd3_(integer *nl, integer *nr, integer *sqre, integer *k, real *d, real *q, integer *ldq, real *dsigma, real *u, integer *ldu, real *u2, integer *ldu2, real *vt, integer *ldvt, real *vt2, integer *ldvt2, integer *idxc, integer *ctot, real *z, integer *info);

static VALUE sHelp, sUsage;

static VALUE
rblapack_slasd3(int argc, VALUE *argv, VALUE self){
  VALUE rblapack_nl;
  integer nl; 
  VALUE rblapack_nr;
  integer nr; 
  VALUE rblapack_sqre;
  integer sqre; 
  VALUE rblapack_dsigma;
  real *dsigma; 
  VALUE rblapack_u2;
  real *u2; 
  VALUE rblapack_vt2;
  real *vt2; 
  VALUE rblapack_idxc;
  integer *idxc; 
  VALUE rblapack_ctot;
  integer *ctot; 
  VALUE rblapack_z;
  real *z; 
  VALUE rblapack_d;
  real *d; 
  VALUE rblapack_u;
  real *u; 
  VALUE rblapack_vt;
  real *vt; 
  VALUE rblapack_info;
  integer info; 
  VALUE rblapack_dsigma_out__;
  real *dsigma_out__;
  VALUE rblapack_vt2_out__;
  real *vt2_out__;
  VALUE rblapack_z_out__;
  real *z_out__;
  real *q;

  integer k;
  integer ldu2;
  integer n;
  integer ldvt2;
  integer ldu;
  integer ldvt;
  integer m;
  integer ldq;

  VALUE rb_options;
  if (argc > 0 && TYPE(argv[argc-1]) == T_HASH) {
    argc--;
    rb_options = argv[argc];
    if (rb_hash_aref(rb_options, sHelp) == Qtrue) {
      printf("%s\n", "USAGE:\n  d, u, vt, info, dsigma, vt2, z = NumRu::Lapack.slasd3( nl, nr, sqre, dsigma, u2, vt2, idxc, ctot, z, [:usage => usage, :help => help])\n\n\nFORTRAN MANUAL\n      SUBROUTINE SLASD3( NL, NR, SQRE, K, D, Q, LDQ, DSIGMA, U, LDU, U2, LDU2, VT, LDVT, VT2, LDVT2, IDXC, CTOT, Z, INFO )\n\n*  Purpose\n*  =======\n*\n*  SLASD3 finds all the square roots of the roots of the secular\n*  equation, as defined by the values in D and Z.  It makes the\n*  appropriate calls to SLASD4 and then updates the singular\n*  vectors by matrix multiplication.\n*\n*  This code makes very mild assumptions about floating point\n*  arithmetic. It will work on machines with a guard digit in\n*  add/subtract, or on those binary machines without guard digits\n*  which subtract like the Cray XMP, Cray YMP, Cray C 90, or Cray 2.\n*  It could conceivably fail on hexadecimal or decimal machines\n*  without guard digits, but we know of none.\n*\n*  SLASD3 is called from SLASD1.\n*\n\n*  Arguments\n*  =========\n*\n*  NL     (input) INTEGER\n*         The row dimension of the upper block.  NL >= 1.\n*\n*  NR     (input) INTEGER\n*         The row dimension of the lower block.  NR >= 1.\n*\n*  SQRE   (input) INTEGER\n*         = 0: the lower block is an NR-by-NR square matrix.\n*         = 1: the lower block is an NR-by-(NR+1) rectangular matrix.\n*\n*         The bidiagonal matrix has N = NL + NR + 1 rows and\n*         M = N + SQRE >= N columns.\n*\n*  K      (input) INTEGER\n*         The size of the secular equation, 1 =< K = < N.\n*\n*  D      (output) REAL array, dimension(K)\n*         On exit the square roots of the roots of the secular equation,\n*         in ascending order.\n*\n*  Q      (workspace) REAL array,\n*                     dimension at least (LDQ,K).\n*\n*  LDQ    (input) INTEGER\n*         The leading dimension of the array Q.  LDQ >= K.\n*\n*  DSIGMA (input/output) REAL array, dimension(K)\n*         The first K elements of this array contain the old roots\n*         of the deflated updating problem.  These are the poles\n*         of the secular equation.\n*\n*  U      (output) REAL array, dimension (LDU, N)\n*         The last N - K columns of this matrix contain the deflated\n*         left singular vectors.\n*\n*  LDU    (input) INTEGER\n*         The leading dimension of the array U.  LDU >= N.\n*\n*  U2     (input) REAL array, dimension (LDU2, N)\n*         The first K columns of this matrix contain the non-deflated\n*         left singular vectors for the split problem.\n*\n*  LDU2   (input) INTEGER\n*         The leading dimension of the array U2.  LDU2 >= N.\n*\n*  VT     (output) REAL array, dimension (LDVT, M)\n*         The last M - K columns of VT' contain the deflated\n*         right singular vectors.\n*\n*  LDVT   (input) INTEGER\n*         The leading dimension of the array VT.  LDVT >= N.\n*\n*  VT2    (input/output) REAL array, dimension (LDVT2, N)\n*         The first K columns of VT2' contain the non-deflated\n*         right singular vectors for the split problem.\n*\n*  LDVT2  (input) INTEGER\n*         The leading dimension of the array VT2.  LDVT2 >= N.\n*\n*  IDXC   (input) INTEGER array, dimension (N)\n*         The permutation used to arrange the columns of U (and rows of\n*         VT) into three groups:  the first group contains non-zero\n*         entries only at and above (or before) NL +1; the second\n*         contains non-zero entries only at and below (or after) NL+2;\n*         and the third is dense. The first column of U and the row of\n*         VT are treated separately, however.\n*\n*         The rows of the singular vectors found by SLASD4\n*         must be likewise permuted before the matrix multiplies can\n*         take place.\n*\n*  CTOT   (input) INTEGER array, dimension (4)\n*         A count of the total number of the various types of columns\n*         in U (or rows in VT), as described in IDXC. The fourth column\n*         type is any column which has been deflated.\n*\n*  Z      (input/output) REAL array, dimension (K)\n*         The first K elements of this array contain the components\n*         of the deflation-adjusted updating row vector.\n*\n*  INFO   (output) INTEGER\n*         = 0:  successful exit.\n*         < 0:  if INFO = -i, the i-th argument had an illegal value.\n*         > 0:  if INFO = 1, a singular value did not converge\n*\n\n*  Further Details\n*  ===============\n*\n*  Based on contributions by\n*     Ming Gu and Huan Ren, Computer Science Division, University of\n*     California at Berkeley, USA\n*\n*  =====================================================================\n*\n\n");
      rb_exit(0);
    }
    if (rb_hash_aref(rb_options, sUsage) == Qtrue) {
      printf("%s\n", "USAGE:\n  d, u, vt, info, dsigma, vt2, z = NumRu::Lapack.slasd3( nl, nr, sqre, dsigma, u2, vt2, idxc, ctot, z, [:usage => usage, :help => help])\n");
      rb_exit(0);
    } 
  } else
    rb_options = Qnil;
  if (argc != 9)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 9)", argc);
  rblapack_nl = argv[0];
  rblapack_nr = argv[1];
  rblapack_sqre = argv[2];
  rblapack_dsigma = argv[3];
  rblapack_u2 = argv[4];
  rblapack_vt2 = argv[5];
  rblapack_idxc = argv[6];
  rblapack_ctot = argv[7];
  rblapack_z = argv[8];
  if (rb_options != Qnil) {
  }

  if (!NA_IsNArray(rblapack_ctot))
    rb_raise(rb_eArgError, "ctot (8th argument) must be NArray");
  if (NA_RANK(rblapack_ctot) != 1)
    rb_raise(rb_eArgError, "rank of ctot (8th argument) must be %d", 1);
  if (NA_SHAPE0(rblapack_ctot) != (4))
    rb_raise(rb_eRuntimeError, "shape 0 of ctot must be %d", 4);
  if (NA_TYPE(rblapack_ctot) != NA_LINT)
    rblapack_ctot = na_change_type(rblapack_ctot, NA_LINT);
  ctot = NA_PTR_TYPE(rblapack_ctot, integer*);
  if (!NA_IsNArray(rblapack_vt2))
    rb_raise(rb_eArgError, "vt2 (6th argument) must be NArray");
  if (NA_RANK(rblapack_vt2) != 2)
    rb_raise(rb_eArgError, "rank of vt2 (6th argument) must be %d", 2);
  n = NA_SHAPE1(rblapack_vt2);
  if (n != (nl + nr + 1))
    rb_raise(rb_eRuntimeError, "shape 1 of vt2 must be %d", nl + nr + 1);
  ldvt2 = NA_SHAPE0(rblapack_vt2);
  if (ldvt2 != (n))
    rb_raise(rb_eRuntimeError, "shape 0 of vt2 must be %d", n);
  n = ldvt2;
  if (NA_TYPE(rblapack_vt2) != NA_SFLOAT)
    rblapack_vt2 = na_change_type(rblapack_vt2, NA_SFLOAT);
  vt2 = NA_PTR_TYPE(rblapack_vt2, real*);
  nl = NUM2INT(rblapack_nl);
  if (!NA_IsNArray(rblapack_idxc))
    rb_raise(rb_eArgError, "idxc (7th argument) must be NArray");
  if (NA_RANK(rblapack_idxc) != 1)
    rb_raise(rb_eArgError, "rank of idxc (7th argument) must be %d", 1);
  if (NA_SHAPE0(rblapack_idxc) != n)
    rb_raise(rb_eRuntimeError, "shape 0 of idxc must be n");
  if (NA_TYPE(rblapack_idxc) != NA_LINT)
    rblapack_idxc = na_change_type(rblapack_idxc, NA_LINT);
  idxc = NA_PTR_TYPE(rblapack_idxc, integer*);
  if (!NA_IsNArray(rblapack_u2))
    rb_raise(rb_eArgError, "u2 (5th argument) must be NArray");
  if (NA_RANK(rblapack_u2) != 2)
    rb_raise(rb_eArgError, "rank of u2 (5th argument) must be %d", 2);
  if (NA_SHAPE1(rblapack_u2) != n)
    rb_raise(rb_eRuntimeError, "shape 1 of u2 must be n");
  ldu2 = NA_SHAPE0(rblapack_u2);
  if (ldu2 != (n))
    rb_raise(rb_eRuntimeError, "shape 0 of u2 must be %d", n);
  n = ldu2;
  if (NA_TYPE(rblapack_u2) != NA_SFLOAT)
    rblapack_u2 = na_change_type(rblapack_u2, NA_SFLOAT);
  u2 = NA_PTR_TYPE(rblapack_u2, real*);
  if (!NA_IsNArray(rblapack_dsigma))
    rb_raise(rb_eArgError, "dsigma (4th argument) must be NArray");
  if (NA_RANK(rblapack_dsigma) != 1)
    rb_raise(rb_eArgError, "rank of dsigma (4th argument) must be %d", 1);
  k = NA_SHAPE0(rblapack_dsigma);
  if (NA_TYPE(rblapack_dsigma) != NA_SFLOAT)
    rblapack_dsigma = na_change_type(rblapack_dsigma, NA_SFLOAT);
  dsigma = NA_PTR_TYPE(rblapack_dsigma, real*);
  sqre = NUM2INT(rblapack_sqre);
  if (!NA_IsNArray(rblapack_z))
    rb_raise(rb_eArgError, "z (9th argument) must be NArray");
  if (NA_RANK(rblapack_z) != 1)
    rb_raise(rb_eArgError, "rank of z (9th argument) must be %d", 1);
  if (NA_SHAPE0(rblapack_z) != k)
    rb_raise(rb_eRuntimeError, "shape 0 of z must be the same as shape 0 of dsigma");
  if (NA_TYPE(rblapack_z) != NA_SFLOAT)
    rblapack_z = na_change_type(rblapack_z, NA_SFLOAT);
  z = NA_PTR_TYPE(rblapack_z, real*);
  nr = NUM2INT(rblapack_nr);
  ldu2 = n;
  ldu = n;
  ldvt = n;
  ldvt2 = n;
  ldq = k;
  m = n+sqre;
  {
    int shape[1];
    shape[0] = k;
    rblapack_d = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  d = NA_PTR_TYPE(rblapack_d, real*);
  {
    int shape[2];
    shape[0] = ldu;
    shape[1] = n;
    rblapack_u = na_make_object(NA_SFLOAT, 2, shape, cNArray);
  }
  u = NA_PTR_TYPE(rblapack_u, real*);
  {
    int shape[2];
    shape[0] = ldvt;
    shape[1] = m;
    rblapack_vt = na_make_object(NA_SFLOAT, 2, shape, cNArray);
  }
  vt = NA_PTR_TYPE(rblapack_vt, real*);
  {
    int shape[1];
    shape[0] = k;
    rblapack_dsigma_out__ = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  dsigma_out__ = NA_PTR_TYPE(rblapack_dsigma_out__, real*);
  MEMCPY(dsigma_out__, dsigma, real, NA_TOTAL(rblapack_dsigma));
  rblapack_dsigma = rblapack_dsigma_out__;
  dsigma = dsigma_out__;
  {
    int shape[2];
    shape[0] = ldvt2;
    shape[1] = n;
    rblapack_vt2_out__ = na_make_object(NA_SFLOAT, 2, shape, cNArray);
  }
  vt2_out__ = NA_PTR_TYPE(rblapack_vt2_out__, real*);
  MEMCPY(vt2_out__, vt2, real, NA_TOTAL(rblapack_vt2));
  rblapack_vt2 = rblapack_vt2_out__;
  vt2 = vt2_out__;
  {
    int shape[1];
    shape[0] = k;
    rblapack_z_out__ = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  z_out__ = NA_PTR_TYPE(rblapack_z_out__, real*);
  MEMCPY(z_out__, z, real, NA_TOTAL(rblapack_z));
  rblapack_z = rblapack_z_out__;
  z = z_out__;
  q = ALLOC_N(real, (ldq)*(k));

  slasd3_(&nl, &nr, &sqre, &k, d, q, &ldq, dsigma, u, &ldu, u2, &ldu2, vt, &ldvt, vt2, &ldvt2, idxc, ctot, z, &info);

  free(q);
  rblapack_info = INT2NUM(info);
  return rb_ary_new3(7, rblapack_d, rblapack_u, rblapack_vt, rblapack_info, rblapack_dsigma, rblapack_vt2, rblapack_z);
}

void
init_lapack_slasd3(VALUE mLapack){
  rb_define_module_function(mLapack, "slasd3", rblapack_slasd3, -1);
  sHelp = ID2SYM(rb_intern("help"));
  sUsage = ID2SYM(rb_intern("usage"));
}
