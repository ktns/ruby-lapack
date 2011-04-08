#include "rb_lapack.h"

extern VOID ctgexc_(logical *wantq, logical *wantz, integer *n, complex *a, integer *lda, complex *b, integer *ldb, complex *q, integer *ldq, complex *z, integer *ldz, integer *ifst, integer *ilst, integer *info);

static VALUE sHelp, sUsage;

static VALUE
rblapack_ctgexc(int argc, VALUE *argv, VALUE self){
  VALUE rblapack_wantq;
  logical wantq; 
  VALUE rblapack_wantz;
  logical wantz; 
  VALUE rblapack_a;
  complex *a; 
  VALUE rblapack_b;
  complex *b; 
  VALUE rblapack_q;
  complex *q; 
  VALUE rblapack_ldq;
  integer ldq; 
  VALUE rblapack_z;
  complex *z; 
  VALUE rblapack_ifst;
  integer ifst; 
  VALUE rblapack_ilst;
  integer ilst; 
  VALUE rblapack_info;
  integer info; 
  VALUE rblapack_a_out__;
  complex *a_out__;
  VALUE rblapack_b_out__;
  complex *b_out__;
  VALUE rblapack_q_out__;
  complex *q_out__;
  VALUE rblapack_z_out__;
  complex *z_out__;

  integer lda;
  integer n;
  integer ldb;
  integer ldz;

  VALUE rb_options;
  if (argc > 0 && TYPE(argv[argc-1]) == T_HASH) {
    argc--;
    rb_options = argv[argc];
    if (rb_hash_aref(rb_options, sHelp) == Qtrue) {
      printf("%s\n", "USAGE:\n  info, a, b, q, z, ilst = NumRu::Lapack.ctgexc( wantq, wantz, a, b, q, ldq, z, ifst, ilst, [:usage => usage, :help => help])\n\n\nFORTRAN MANUAL\n      SUBROUTINE CTGEXC( WANTQ, WANTZ, N, A, LDA, B, LDB, Q, LDQ, Z, LDZ, IFST, ILST, INFO )\n\n*  Purpose\n*  =======\n*\n*  CTGEXC reorders the generalized Schur decomposition of a complex\n*  matrix pair (A,B), using an unitary equivalence transformation\n*  (A, B) := Q * (A, B) * Z', so that the diagonal block of (A, B) with\n*  row index IFST is moved to row ILST.\n*\n*  (A, B) must be in generalized Schur canonical form, that is, A and\n*  B are both upper triangular.\n*\n*  Optionally, the matrices Q and Z of generalized Schur vectors are\n*  updated.\n*\n*         Q(in) * A(in) * Z(in)' = Q(out) * A(out) * Z(out)'\n*         Q(in) * B(in) * Z(in)' = Q(out) * B(out) * Z(out)'\n*\n\n*  Arguments\n*  =========\n*\n*  WANTQ   (input) LOGICAL\n*          .TRUE. : update the left transformation matrix Q;\n*          .FALSE.: do not update Q.\n*\n*  WANTZ   (input) LOGICAL\n*          .TRUE. : update the right transformation matrix Z;\n*          .FALSE.: do not update Z.\n*\n*  N       (input) INTEGER\n*          The order of the matrices A and B. N >= 0.\n*\n*  A       (input/output) COMPLEX array, dimension (LDA,N)\n*          On entry, the upper triangular matrix A in the pair (A, B).\n*          On exit, the updated matrix A.\n*\n*  LDA     (input)  INTEGER\n*          The leading dimension of the array A. LDA >= max(1,N).\n*\n*  B       (input/output) COMPLEX array, dimension (LDB,N)\n*          On entry, the upper triangular matrix B in the pair (A, B).\n*          On exit, the updated matrix B.\n*\n*  LDB     (input)  INTEGER\n*          The leading dimension of the array B. LDB >= max(1,N).\n*\n*  Q       (input/output) COMPLEX array, dimension (LDZ,N)\n*          On entry, if WANTQ = .TRUE., the unitary matrix Q.\n*          On exit, the updated matrix Q.\n*          If WANTQ = .FALSE., Q is not referenced.\n*\n*  LDQ     (input) INTEGER\n*          The leading dimension of the array Q. LDQ >= 1;\n*          If WANTQ = .TRUE., LDQ >= N.\n*\n*  Z       (input/output) COMPLEX array, dimension (LDZ,N)\n*          On entry, if WANTZ = .TRUE., the unitary matrix Z.\n*          On exit, the updated matrix Z.\n*          If WANTZ = .FALSE., Z is not referenced.\n*\n*  LDZ     (input) INTEGER\n*          The leading dimension of the array Z. LDZ >= 1;\n*          If WANTZ = .TRUE., LDZ >= N.\n*\n*  IFST    (input) INTEGER\n*  ILST    (input/output) INTEGER\n*          Specify the reordering of the diagonal blocks of (A, B).\n*          The block with row index IFST is moved to row ILST, by a\n*          sequence of swapping between adjacent blocks.\n*\n*  INFO    (output) INTEGER\n*           =0:  Successful exit.\n*           <0:  if INFO = -i, the i-th argument had an illegal value.\n*           =1:  The transformed matrix pair (A, B) would be too far\n*                from generalized Schur form; the problem is ill-\n*                conditioned. (A, B) may have been partially reordered,\n*                and ILST points to the first row of the current\n*                position of the block being moved.\n*\n*\n\n*  Further Details\n*  ===============\n*\n*  Based on contributions by\n*     Bo Kagstrom and Peter Poromaa, Department of Computing Science,\n*     Umea University, S-901 87 Umea, Sweden.\n*\n*  [1] B. Kagstrom; A Direct Method for Reordering Eigenvalues in the\n*      Generalized Real Schur Form of a Regular Matrix Pair (A, B), in\n*      M.S. Moonen et al (eds), Linear Algebra for Large Scale and\n*      Real-Time Applications, Kluwer Academic Publ. 1993, pp 195-218.\n*\n*  [2] B. Kagstrom and P. Poromaa; Computing Eigenspaces with Specified\n*      Eigenvalues of a Regular Matrix Pair (A, B) and Condition\n*      Estimation: Theory, Algorithms and Software, Report\n*      UMINF - 94.04, Department of Computing Science, Umea University,\n*      S-901 87 Umea, Sweden, 1994. Also as LAPACK Working Note 87.\n*      To appear in Numerical Algorithms, 1996.\n*\n*  [3] B. Kagstrom and P. Poromaa, LAPACK-Style Algorithms and Software\n*      for Solving the Generalized Sylvester Equation and Estimating the\n*      Separation between Regular Matrix Pairs, Report UMINF - 93.23,\n*      Department of Computing Science, Umea University, S-901 87 Umea,\n*      Sweden, December 1993, Revised April 1994, Also as LAPACK working\n*      Note 75. To appear in ACM Trans. on Math. Software, Vol 22, No 1,\n*      1996.\n*\n*  =====================================================================\n*\n*     .. Local Scalars ..\n      INTEGER            HERE\n*     ..\n*     .. External Subroutines ..\n      EXTERNAL           CTGEX2, XERBLA\n*     ..\n*     .. Intrinsic Functions ..\n      INTRINSIC          MAX\n*     ..\n\n");
      rb_exit(0);
    }
    if (rb_hash_aref(rb_options, sUsage) == Qtrue) {
      printf("%s\n", "USAGE:\n  info, a, b, q, z, ilst = NumRu::Lapack.ctgexc( wantq, wantz, a, b, q, ldq, z, ifst, ilst, [:usage => usage, :help => help])\n");
      rb_exit(0);
    } 
  } else
    rb_options = Qnil;
  if (argc != 9)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 9)", argc);
  rblapack_wantq = argv[0];
  rblapack_wantz = argv[1];
  rblapack_a = argv[2];
  rblapack_b = argv[3];
  rblapack_q = argv[4];
  rblapack_ldq = argv[5];
  rblapack_z = argv[6];
  rblapack_ifst = argv[7];
  rblapack_ilst = argv[8];
  if (rb_options != Qnil) {
  }

  if (!NA_IsNArray(rblapack_a))
    rb_raise(rb_eArgError, "a (3th argument) must be NArray");
  if (NA_RANK(rblapack_a) != 2)
    rb_raise(rb_eArgError, "rank of a (3th argument) must be %d", 2);
  n = NA_SHAPE1(rblapack_a);
  lda = NA_SHAPE0(rblapack_a);
  if (NA_TYPE(rblapack_a) != NA_SCOMPLEX)
    rblapack_a = na_change_type(rblapack_a, NA_SCOMPLEX);
  a = NA_PTR_TYPE(rblapack_a, complex*);
  wantz = (rblapack_wantz == Qtrue);
  if (!NA_IsNArray(rblapack_b))
    rb_raise(rb_eArgError, "b (4th argument) must be NArray");
  if (NA_RANK(rblapack_b) != 2)
    rb_raise(rb_eArgError, "rank of b (4th argument) must be %d", 2);
  if (NA_SHAPE1(rblapack_b) != n)
    rb_raise(rb_eRuntimeError, "shape 1 of b must be the same as shape 1 of a");
  ldb = NA_SHAPE0(rblapack_b);
  if (NA_TYPE(rblapack_b) != NA_SCOMPLEX)
    rblapack_b = na_change_type(rblapack_b, NA_SCOMPLEX);
  b = NA_PTR_TYPE(rblapack_b, complex*);
  ldq = NUM2INT(rblapack_ldq);
  wantq = (rblapack_wantq == Qtrue);
  if (!NA_IsNArray(rblapack_z))
    rb_raise(rb_eArgError, "z (7th argument) must be NArray");
  if (NA_RANK(rblapack_z) != 2)
    rb_raise(rb_eArgError, "rank of z (7th argument) must be %d", 2);
  if (NA_SHAPE1(rblapack_z) != n)
    rb_raise(rb_eRuntimeError, "shape 1 of z must be the same as shape 1 of a");
  ldz = NA_SHAPE0(rblapack_z);
  if (NA_TYPE(rblapack_z) != NA_SCOMPLEX)
    rblapack_z = na_change_type(rblapack_z, NA_SCOMPLEX);
  z = NA_PTR_TYPE(rblapack_z, complex*);
  ilst = NUM2INT(rblapack_ilst);
  if (!NA_IsNArray(rblapack_q))
    rb_raise(rb_eArgError, "q (5th argument) must be NArray");
  if (NA_RANK(rblapack_q) != 2)
    rb_raise(rb_eArgError, "rank of q (5th argument) must be %d", 2);
  if (NA_SHAPE1(rblapack_q) != n)
    rb_raise(rb_eRuntimeError, "shape 1 of q must be the same as shape 1 of a");
  if (NA_SHAPE0(rblapack_q) != ldz)
    rb_raise(rb_eRuntimeError, "shape 0 of q must be the same as shape 0 of z");
  if (NA_TYPE(rblapack_q) != NA_SCOMPLEX)
    rblapack_q = na_change_type(rblapack_q, NA_SCOMPLEX);
  q = NA_PTR_TYPE(rblapack_q, complex*);
  ifst = NUM2INT(rblapack_ifst);
  {
    int shape[2];
    shape[0] = lda;
    shape[1] = n;
    rblapack_a_out__ = na_make_object(NA_SCOMPLEX, 2, shape, cNArray);
  }
  a_out__ = NA_PTR_TYPE(rblapack_a_out__, complex*);
  MEMCPY(a_out__, a, complex, NA_TOTAL(rblapack_a));
  rblapack_a = rblapack_a_out__;
  a = a_out__;
  {
    int shape[2];
    shape[0] = ldb;
    shape[1] = n;
    rblapack_b_out__ = na_make_object(NA_SCOMPLEX, 2, shape, cNArray);
  }
  b_out__ = NA_PTR_TYPE(rblapack_b_out__, complex*);
  MEMCPY(b_out__, b, complex, NA_TOTAL(rblapack_b));
  rblapack_b = rblapack_b_out__;
  b = b_out__;
  {
    int shape[2];
    shape[0] = ldz;
    shape[1] = n;
    rblapack_q_out__ = na_make_object(NA_SCOMPLEX, 2, shape, cNArray);
  }
  q_out__ = NA_PTR_TYPE(rblapack_q_out__, complex*);
  MEMCPY(q_out__, q, complex, NA_TOTAL(rblapack_q));
  rblapack_q = rblapack_q_out__;
  q = q_out__;
  {
    int shape[2];
    shape[0] = ldz;
    shape[1] = n;
    rblapack_z_out__ = na_make_object(NA_SCOMPLEX, 2, shape, cNArray);
  }
  z_out__ = NA_PTR_TYPE(rblapack_z_out__, complex*);
  MEMCPY(z_out__, z, complex, NA_TOTAL(rblapack_z));
  rblapack_z = rblapack_z_out__;
  z = z_out__;

  ctgexc_(&wantq, &wantz, &n, a, &lda, b, &ldb, q, &ldq, z, &ldz, &ifst, &ilst, &info);

  rblapack_info = INT2NUM(info);
  rblapack_ilst = INT2NUM(ilst);
  return rb_ary_new3(6, rblapack_info, rblapack_a, rblapack_b, rblapack_q, rblapack_z, rblapack_ilst);
}

void
init_lapack_ctgexc(VALUE mLapack){
  rb_define_module_function(mLapack, "ctgexc", rblapack_ctgexc, -1);
  sHelp = ID2SYM(rb_intern("help"));
  sUsage = ID2SYM(rb_intern("usage"));
}
