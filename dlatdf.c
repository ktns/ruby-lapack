#include "rb_lapack.h"

extern VOID dlatdf_(integer *ijob, integer *n, doublereal *z, integer *ldz, doublereal *rhs, doublereal *rdsum, doublereal *rdscal, integer *ipiv, integer *jpiv);

static VALUE sHelp, sUsage;

static VALUE
rblapack_dlatdf(int argc, VALUE *argv, VALUE self){
  VALUE rblapack_ijob;
  integer ijob; 
  VALUE rblapack_z;
  doublereal *z; 
  VALUE rblapack_rhs;
  doublereal *rhs; 
  VALUE rblapack_rdsum;
  doublereal rdsum; 
  VALUE rblapack_rdscal;
  doublereal rdscal; 
  VALUE rblapack_ipiv;
  integer *ipiv; 
  VALUE rblapack_jpiv;
  integer *jpiv; 
  VALUE rblapack_rhs_out__;
  doublereal *rhs_out__;

  integer ldz;
  integer n;

  VALUE rb_options;
  if (argc > 0 && TYPE(argv[argc-1]) == T_HASH) {
    argc--;
    rb_options = argv[argc];
    if (rb_hash_aref(rb_options, sHelp) == Qtrue) {
      printf("%s\n", "USAGE:\n  rhs, rdsum, rdscal = NumRu::Lapack.dlatdf( ijob, z, rhs, rdsum, rdscal, ipiv, jpiv, [:usage => usage, :help => help])\n\n\nFORTRAN MANUAL\n      SUBROUTINE DLATDF( IJOB, N, Z, LDZ, RHS, RDSUM, RDSCAL, IPIV, JPIV )\n\n*  Purpose\n*  =======\n*\n*  DLATDF uses the LU factorization of the n-by-n matrix Z computed by\n*  DGETC2 and computes a contribution to the reciprocal Dif-estimate\n*  by solving Z * x = b for x, and choosing the r.h.s. b such that\n*  the norm of x is as large as possible. On entry RHS = b holds the\n*  contribution from earlier solved sub-systems, and on return RHS = x.\n*\n*  The factorization of Z returned by DGETC2 has the form Z = P*L*U*Q,\n*  where P and Q are permutation matrices. L is lower triangular with\n*  unit diagonal elements and U is upper triangular.\n*\n\n*  Arguments\n*  =========\n*\n*  IJOB    (input) INTEGER\n*          IJOB = 2: First compute an approximative null-vector e\n*              of Z using DGECON, e is normalized and solve for\n*              Zx = +-e - f with the sign giving the greater value\n*              of 2-norm(x). About 5 times as expensive as Default.\n*          IJOB .ne. 2: Local look ahead strategy where all entries of\n*              the r.h.s. b is choosen as either +1 or -1 (Default).\n*\n*  N       (input) INTEGER\n*          The number of columns of the matrix Z.\n*\n*  Z       (input) DOUBLE PRECISION array, dimension (LDZ, N)\n*          On entry, the LU part of the factorization of the n-by-n\n*          matrix Z computed by DGETC2:  Z = P * L * U * Q\n*\n*  LDZ     (input) INTEGER\n*          The leading dimension of the array Z.  LDA >= max(1, N).\n*\n*  RHS     (input/output) DOUBLE PRECISION array, dimension (N)\n*          On entry, RHS contains contributions from other subsystems.\n*          On exit, RHS contains the solution of the subsystem with\n*          entries acoording to the value of IJOB (see above).\n*\n*  RDSUM   (input/output) DOUBLE PRECISION\n*          On entry, the sum of squares of computed contributions to\n*          the Dif-estimate under computation by DTGSYL, where the\n*          scaling factor RDSCAL (see below) has been factored out.\n*          On exit, the corresponding sum of squares updated with the\n*          contributions from the current sub-system.\n*          If TRANS = 'T' RDSUM is not touched.\n*          NOTE: RDSUM only makes sense when DTGSY2 is called by STGSYL.\n*\n*  RDSCAL  (input/output) DOUBLE PRECISION\n*          On entry, scaling factor used to prevent overflow in RDSUM.\n*          On exit, RDSCAL is updated w.r.t. the current contributions\n*          in RDSUM.\n*          If TRANS = 'T', RDSCAL is not touched.\n*          NOTE: RDSCAL only makes sense when DTGSY2 is called by\n*                DTGSYL.\n*\n*  IPIV    (input) INTEGER array, dimension (N).\n*          The pivot indices; for 1 <= i <= N, row i of the\n*          matrix has been interchanged with row IPIV(i).\n*\n*  JPIV    (input) INTEGER array, dimension (N).\n*          The pivot indices; for 1 <= j <= N, column j of the\n*          matrix has been interchanged with column JPIV(j).\n*\n\n*  Further Details\n*  ===============\n*\n*  Based on contributions by\n*     Bo Kagstrom and Peter Poromaa, Department of Computing Science,\n*     Umea University, S-901 87 Umea, Sweden.\n*\n*  This routine is a further developed implementation of algorithm\n*  BSOLVE in [1] using complete pivoting in the LU factorization.\n*\n*  [1] Bo Kagstrom and Lars Westin,\n*      Generalized Schur Methods with Condition Estimators for\n*      Solving the Generalized Sylvester Equation, IEEE Transactions\n*      on Automatic Control, Vol. 34, No. 7, July 1989, pp 745-751.\n*\n*  [2] Peter Poromaa,\n*      On Efficient and Robust Estimators for the Separation\n*      between two Regular Matrix Pairs with Applications in\n*      Condition Estimation. Report IMINF-95.05, Departement of\n*      Computing Science, Umea University, S-901 87 Umea, Sweden, 1995.\n*\n*  =====================================================================\n*\n\n");
      rb_exit(0);
    }
    if (rb_hash_aref(rb_options, sUsage) == Qtrue) {
      printf("%s\n", "USAGE:\n  rhs, rdsum, rdscal = NumRu::Lapack.dlatdf( ijob, z, rhs, rdsum, rdscal, ipiv, jpiv, [:usage => usage, :help => help])\n");
      rb_exit(0);
    } 
  } else
    rb_options = Qnil;
  if (argc != 7)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 7)", argc);
  rblapack_ijob = argv[0];
  rblapack_z = argv[1];
  rblapack_rhs = argv[2];
  rblapack_rdsum = argv[3];
  rblapack_rdscal = argv[4];
  rblapack_ipiv = argv[5];
  rblapack_jpiv = argv[6];
  if (rb_options != Qnil) {
  }

  if (!NA_IsNArray(rblapack_ipiv))
    rb_raise(rb_eArgError, "ipiv (6th argument) must be NArray");
  if (NA_RANK(rblapack_ipiv) != 1)
    rb_raise(rb_eArgError, "rank of ipiv (6th argument) must be %d", 1);
  n = NA_SHAPE0(rblapack_ipiv);
  if (NA_TYPE(rblapack_ipiv) != NA_LINT)
    rblapack_ipiv = na_change_type(rblapack_ipiv, NA_LINT);
  ipiv = NA_PTR_TYPE(rblapack_ipiv, integer*);
  rdscal = NUM2DBL(rblapack_rdscal);
  ijob = NUM2INT(rblapack_ijob);
  if (!NA_IsNArray(rblapack_z))
    rb_raise(rb_eArgError, "z (2th argument) must be NArray");
  if (NA_RANK(rblapack_z) != 2)
    rb_raise(rb_eArgError, "rank of z (2th argument) must be %d", 2);
  if (NA_SHAPE1(rblapack_z) != n)
    rb_raise(rb_eRuntimeError, "shape 1 of z must be the same as shape 0 of ipiv");
  ldz = NA_SHAPE0(rblapack_z);
  if (NA_TYPE(rblapack_z) != NA_DFLOAT)
    rblapack_z = na_change_type(rblapack_z, NA_DFLOAT);
  z = NA_PTR_TYPE(rblapack_z, doublereal*);
  rdsum = NUM2DBL(rblapack_rdsum);
  if (!NA_IsNArray(rblapack_rhs))
    rb_raise(rb_eArgError, "rhs (3th argument) must be NArray");
  if (NA_RANK(rblapack_rhs) != 1)
    rb_raise(rb_eArgError, "rank of rhs (3th argument) must be %d", 1);
  if (NA_SHAPE0(rblapack_rhs) != n)
    rb_raise(rb_eRuntimeError, "shape 0 of rhs must be the same as shape 0 of ipiv");
  if (NA_TYPE(rblapack_rhs) != NA_DFLOAT)
    rblapack_rhs = na_change_type(rblapack_rhs, NA_DFLOAT);
  rhs = NA_PTR_TYPE(rblapack_rhs, doublereal*);
  if (!NA_IsNArray(rblapack_jpiv))
    rb_raise(rb_eArgError, "jpiv (7th argument) must be NArray");
  if (NA_RANK(rblapack_jpiv) != 1)
    rb_raise(rb_eArgError, "rank of jpiv (7th argument) must be %d", 1);
  if (NA_SHAPE0(rblapack_jpiv) != n)
    rb_raise(rb_eRuntimeError, "shape 0 of jpiv must be the same as shape 0 of ipiv");
  if (NA_TYPE(rblapack_jpiv) != NA_LINT)
    rblapack_jpiv = na_change_type(rblapack_jpiv, NA_LINT);
  jpiv = NA_PTR_TYPE(rblapack_jpiv, integer*);
  {
    int shape[1];
    shape[0] = n;
    rblapack_rhs_out__ = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  rhs_out__ = NA_PTR_TYPE(rblapack_rhs_out__, doublereal*);
  MEMCPY(rhs_out__, rhs, doublereal, NA_TOTAL(rblapack_rhs));
  rblapack_rhs = rblapack_rhs_out__;
  rhs = rhs_out__;

  dlatdf_(&ijob, &n, z, &ldz, rhs, &rdsum, &rdscal, ipiv, jpiv);

  rblapack_rdsum = rb_float_new((double)rdsum);
  rblapack_rdscal = rb_float_new((double)rdscal);
  return rb_ary_new3(3, rblapack_rhs, rblapack_rdsum, rblapack_rdscal);
}

void
init_lapack_dlatdf(VALUE mLapack){
  rb_define_module_function(mLapack, "dlatdf", rblapack_dlatdf, -1);
  sHelp = ID2SYM(rb_intern("help"));
  sUsage = ID2SYM(rb_intern("usage"));
}
