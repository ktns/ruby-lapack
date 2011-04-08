#include "rb_lapack.h"

extern VOID sgghrd_(char *compq, char *compz, integer *n, integer *ilo, integer *ihi, real *a, integer *lda, real *b, integer *ldb, real *q, integer *ldq, real *z, integer *ldz, integer *info);

static VALUE sHelp, sUsage;

static VALUE
rblapack_sgghrd(int argc, VALUE *argv, VALUE self){
  VALUE rblapack_compq;
  char compq; 
  VALUE rblapack_compz;
  char compz; 
  VALUE rblapack_ilo;
  integer ilo; 
  VALUE rblapack_ihi;
  integer ihi; 
  VALUE rblapack_a;
  real *a; 
  VALUE rblapack_b;
  real *b; 
  VALUE rblapack_q;
  real *q; 
  VALUE rblapack_z;
  real *z; 
  VALUE rblapack_info;
  integer info; 
  VALUE rblapack_a_out__;
  real *a_out__;
  VALUE rblapack_b_out__;
  real *b_out__;
  VALUE rblapack_q_out__;
  real *q_out__;
  VALUE rblapack_z_out__;
  real *z_out__;

  integer lda;
  integer n;
  integer ldb;
  integer ldq;
  integer ldz;

  VALUE rb_options;
  if (argc > 0 && TYPE(argv[argc-1]) == T_HASH) {
    argc--;
    rb_options = argv[argc];
    if (rb_hash_aref(rb_options, sHelp) == Qtrue) {
      printf("%s\n", "USAGE:\n  info, a, b, q, z = NumRu::Lapack.sgghrd( compq, compz, ilo, ihi, a, b, q, z, [:usage => usage, :help => help])\n\n\nFORTRAN MANUAL\n      SUBROUTINE SGGHRD( COMPQ, COMPZ, N, ILO, IHI, A, LDA, B, LDB, Q, LDQ, Z, LDZ, INFO )\n\n*  Purpose\n*  =======\n*\n*  SGGHRD reduces a pair of real matrices (A,B) to generalized upper\n*  Hessenberg form using orthogonal transformations, where A is a\n*  general matrix and B is upper triangular.  The form of the\n*  generalized eigenvalue problem is\n*     A*x = lambda*B*x,\n*  and B is typically made upper triangular by computing its QR\n*  factorization and moving the orthogonal matrix Q to the left side\n*  of the equation.\n*\n*  This subroutine simultaneously reduces A to a Hessenberg matrix H:\n*     Q**T*A*Z = H\n*  and transforms B to another upper triangular matrix T:\n*     Q**T*B*Z = T\n*  in order to reduce the problem to its standard form\n*     H*y = lambda*T*y\n*  where y = Z**T*x.\n*\n*  The orthogonal matrices Q and Z are determined as products of Givens\n*  rotations.  They may either be formed explicitly, or they may be\n*  postmultiplied into input matrices Q1 and Z1, so that\n*\n*       Q1 * A * Z1**T = (Q1*Q) * H * (Z1*Z)**T\n*\n*       Q1 * B * Z1**T = (Q1*Q) * T * (Z1*Z)**T\n*\n*  If Q1 is the orthogonal matrix from the QR factorization of B in the\n*  original equation A*x = lambda*B*x, then SGGHRD reduces the original\n*  problem to generalized Hessenberg form.\n*\n\n*  Arguments\n*  =========\n*\n*  COMPQ   (input) CHARACTER*1\n*          = 'N': do not compute Q;\n*          = 'I': Q is initialized to the unit matrix, and the\n*                 orthogonal matrix Q is returned;\n*          = 'V': Q must contain an orthogonal matrix Q1 on entry,\n*                 and the product Q1*Q is returned.\n*\n*  COMPZ   (input) CHARACTER*1\n*          = 'N': do not compute Z;\n*          = 'I': Z is initialized to the unit matrix, and the\n*                 orthogonal matrix Z is returned;\n*          = 'V': Z must contain an orthogonal matrix Z1 on entry,\n*                 and the product Z1*Z is returned.\n*\n*  N       (input) INTEGER\n*          The order of the matrices A and B.  N >= 0.\n*\n*  ILO     (input) INTEGER\n*  IHI     (input) INTEGER\n*          ILO and IHI mark the rows and columns of A which are to be\n*          reduced.  It is assumed that A is already upper triangular\n*          in rows and columns 1:ILO-1 and IHI+1:N.  ILO and IHI are\n*          normally set by a previous call to SGGBAL; otherwise they\n*          should be set to 1 and N respectively.\n*          1 <= ILO <= IHI <= N, if N > 0; ILO=1 and IHI=0, if N=0.\n*\n*  A       (input/output) REAL array, dimension (LDA, N)\n*          On entry, the N-by-N general matrix to be reduced.\n*          On exit, the upper triangle and the first subdiagonal of A\n*          are overwritten with the upper Hessenberg matrix H, and the\n*          rest is set to zero.\n*\n*  LDA     (input) INTEGER\n*          The leading dimension of the array A.  LDA >= max(1,N).\n*\n*  B       (input/output) REAL array, dimension (LDB, N)\n*          On entry, the N-by-N upper triangular matrix B.\n*          On exit, the upper triangular matrix T = Q**T B Z.  The\n*          elements below the diagonal are set to zero.\n*\n*  LDB     (input) INTEGER\n*          The leading dimension of the array B.  LDB >= max(1,N).\n*\n*  Q       (input/output) REAL array, dimension (LDQ, N)\n*          On entry, if COMPQ = 'V', the orthogonal matrix Q1,\n*          typically from the QR factorization of B.\n*          On exit, if COMPQ='I', the orthogonal matrix Q, and if\n*          COMPQ = 'V', the product Q1*Q.\n*          Not referenced if COMPQ='N'.\n*\n*  LDQ     (input) INTEGER\n*          The leading dimension of the array Q.\n*          LDQ >= N if COMPQ='V' or 'I'; LDQ >= 1 otherwise.\n*\n*  Z       (input/output) REAL array, dimension (LDZ, N)\n*          On entry, if COMPZ = 'V', the orthogonal matrix Z1.\n*          On exit, if COMPZ='I', the orthogonal matrix Z, and if\n*          COMPZ = 'V', the product Z1*Z.\n*          Not referenced if COMPZ='N'.\n*\n*  LDZ     (input) INTEGER\n*          The leading dimension of the array Z.\n*          LDZ >= N if COMPZ='V' or 'I'; LDZ >= 1 otherwise.\n*\n*  INFO    (output) INTEGER\n*          = 0:  successful exit.\n*          < 0:  if INFO = -i, the i-th argument had an illegal value.\n*\n\n*  Further Details\n*  ===============\n*\n*  This routine reduces A to Hessenberg and B to triangular form by\n*  an unblocked reduction, as described in _Matrix_Computations_,\n*  by Golub and Van Loan (Johns Hopkins Press.)\n*\n*  =====================================================================\n*\n\n");
      rb_exit(0);
    }
    if (rb_hash_aref(rb_options, sUsage) == Qtrue) {
      printf("%s\n", "USAGE:\n  info, a, b, q, z = NumRu::Lapack.sgghrd( compq, compz, ilo, ihi, a, b, q, z, [:usage => usage, :help => help])\n");
      rb_exit(0);
    } 
  } else
    rb_options = Qnil;
  if (argc != 8)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 8)", argc);
  rblapack_compq = argv[0];
  rblapack_compz = argv[1];
  rblapack_ilo = argv[2];
  rblapack_ihi = argv[3];
  rblapack_a = argv[4];
  rblapack_b = argv[5];
  rblapack_q = argv[6];
  rblapack_z = argv[7];
  if (rb_options != Qnil) {
  }

  if (!NA_IsNArray(rblapack_a))
    rb_raise(rb_eArgError, "a (5th argument) must be NArray");
  if (NA_RANK(rblapack_a) != 2)
    rb_raise(rb_eArgError, "rank of a (5th argument) must be %d", 2);
  n = NA_SHAPE1(rblapack_a);
  lda = NA_SHAPE0(rblapack_a);
  if (NA_TYPE(rblapack_a) != NA_SFLOAT)
    rblapack_a = na_change_type(rblapack_a, NA_SFLOAT);
  a = NA_PTR_TYPE(rblapack_a, real*);
  ilo = NUM2INT(rblapack_ilo);
  if (!NA_IsNArray(rblapack_b))
    rb_raise(rb_eArgError, "b (6th argument) must be NArray");
  if (NA_RANK(rblapack_b) != 2)
    rb_raise(rb_eArgError, "rank of b (6th argument) must be %d", 2);
  if (NA_SHAPE1(rblapack_b) != n)
    rb_raise(rb_eRuntimeError, "shape 1 of b must be the same as shape 1 of a");
  ldb = NA_SHAPE0(rblapack_b);
  if (NA_TYPE(rblapack_b) != NA_SFLOAT)
    rblapack_b = na_change_type(rblapack_b, NA_SFLOAT);
  b = NA_PTR_TYPE(rblapack_b, real*);
  compz = StringValueCStr(rblapack_compz)[0];
  if (!NA_IsNArray(rblapack_z))
    rb_raise(rb_eArgError, "z (8th argument) must be NArray");
  if (NA_RANK(rblapack_z) != 2)
    rb_raise(rb_eArgError, "rank of z (8th argument) must be %d", 2);
  if (NA_SHAPE1(rblapack_z) != n)
    rb_raise(rb_eRuntimeError, "shape 1 of z must be the same as shape 1 of a");
  ldz = NA_SHAPE0(rblapack_z);
  if (NA_TYPE(rblapack_z) != NA_SFLOAT)
    rblapack_z = na_change_type(rblapack_z, NA_SFLOAT);
  z = NA_PTR_TYPE(rblapack_z, real*);
  compq = StringValueCStr(rblapack_compq)[0];
  if (!NA_IsNArray(rblapack_q))
    rb_raise(rb_eArgError, "q (7th argument) must be NArray");
  if (NA_RANK(rblapack_q) != 2)
    rb_raise(rb_eArgError, "rank of q (7th argument) must be %d", 2);
  if (NA_SHAPE1(rblapack_q) != n)
    rb_raise(rb_eRuntimeError, "shape 1 of q must be the same as shape 1 of a");
  ldq = NA_SHAPE0(rblapack_q);
  if (NA_TYPE(rblapack_q) != NA_SFLOAT)
    rblapack_q = na_change_type(rblapack_q, NA_SFLOAT);
  q = NA_PTR_TYPE(rblapack_q, real*);
  ihi = NUM2INT(rblapack_ihi);
  {
    int shape[2];
    shape[0] = lda;
    shape[1] = n;
    rblapack_a_out__ = na_make_object(NA_SFLOAT, 2, shape, cNArray);
  }
  a_out__ = NA_PTR_TYPE(rblapack_a_out__, real*);
  MEMCPY(a_out__, a, real, NA_TOTAL(rblapack_a));
  rblapack_a = rblapack_a_out__;
  a = a_out__;
  {
    int shape[2];
    shape[0] = ldb;
    shape[1] = n;
    rblapack_b_out__ = na_make_object(NA_SFLOAT, 2, shape, cNArray);
  }
  b_out__ = NA_PTR_TYPE(rblapack_b_out__, real*);
  MEMCPY(b_out__, b, real, NA_TOTAL(rblapack_b));
  rblapack_b = rblapack_b_out__;
  b = b_out__;
  {
    int shape[2];
    shape[0] = ldq;
    shape[1] = n;
    rblapack_q_out__ = na_make_object(NA_SFLOAT, 2, shape, cNArray);
  }
  q_out__ = NA_PTR_TYPE(rblapack_q_out__, real*);
  MEMCPY(q_out__, q, real, NA_TOTAL(rblapack_q));
  rblapack_q = rblapack_q_out__;
  q = q_out__;
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

  sgghrd_(&compq, &compz, &n, &ilo, &ihi, a, &lda, b, &ldb, q, &ldq, z, &ldz, &info);

  rblapack_info = INT2NUM(info);
  return rb_ary_new3(5, rblapack_info, rblapack_a, rblapack_b, rblapack_q, rblapack_z);
}

void
init_lapack_sgghrd(VALUE mLapack){
  rb_define_module_function(mLapack, "sgghrd", rblapack_sgghrd, -1);
  sHelp = ID2SYM(rb_intern("help"));
  sUsage = ID2SYM(rb_intern("usage"));
}
