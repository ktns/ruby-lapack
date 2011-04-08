#include "rb_lapack.h"

extern VOID ztgex2_(logical *wantq, logical *wantz, integer *n, doublecomplex *a, integer *lda, doublecomplex *b, integer *ldb, doublecomplex *q, integer *ldq, doublecomplex *z, integer *ldz, integer *j1, integer *info);

static VALUE sHelp, sUsage;

static VALUE
rblapack_ztgex2(int argc, VALUE *argv, VALUE self){
  VALUE rblapack_wantq;
  logical wantq; 
  VALUE rblapack_wantz;
  logical wantz; 
  VALUE rblapack_a;
  doublecomplex *a; 
  VALUE rblapack_b;
  doublecomplex *b; 
  VALUE rblapack_q;
  doublecomplex *q; 
  VALUE rblapack_ldq;
  integer ldq; 
  VALUE rblapack_z;
  doublecomplex *z; 
  VALUE rblapack_ldz;
  integer ldz; 
  VALUE rblapack_j1;
  integer j1; 
  VALUE rblapack_info;
  integer info; 
  VALUE rblapack_a_out__;
  doublecomplex *a_out__;
  VALUE rblapack_b_out__;
  doublecomplex *b_out__;
  VALUE rblapack_q_out__;
  doublecomplex *q_out__;
  VALUE rblapack_z_out__;
  doublecomplex *z_out__;

  integer lda;
  integer n;
  integer ldb;

  VALUE rb_options;
  if (argc > 0 && TYPE(argv[argc-1]) == T_HASH) {
    argc--;
    rb_options = argv[argc];
    if (rb_hash_aref(rb_options, sHelp) == Qtrue) {
      printf("%s\n", "USAGE:\n  info, a, b, q, z = NumRu::Lapack.ztgex2( wantq, wantz, a, b, q, ldq, z, ldz, j1, [:usage => usage, :help => help])\n\n\nFORTRAN MANUAL\n      SUBROUTINE ZTGEX2( WANTQ, WANTZ, N, A, LDA, B, LDB, Q, LDQ, Z, LDZ, J1, INFO )\n\n*  Purpose\n*  =======\n*\n*  ZTGEX2 swaps adjacent diagonal 1 by 1 blocks (A11,B11) and (A22,B22)\n*  in an upper triangular matrix pair (A, B) by an unitary equivalence\n*  transformation.\n*\n*  (A, B) must be in generalized Schur canonical form, that is, A and\n*  B are both upper triangular.\n*\n*  Optionally, the matrices Q and Z of generalized Schur vectors are\n*  updated.\n*\n*         Q(in) * A(in) * Z(in)' = Q(out) * A(out) * Z(out)'\n*         Q(in) * B(in) * Z(in)' = Q(out) * B(out) * Z(out)'\n*\n*\n\n*  Arguments\n*  =========\n*\n*  WANTQ   (input) LOGICAL\n*          .TRUE. : update the left transformation matrix Q;\n*          .FALSE.: do not update Q.\n*\n*  WANTZ   (input) LOGICAL\n*          .TRUE. : update the right transformation matrix Z;\n*          .FALSE.: do not update Z.\n*\n*  N       (input) INTEGER\n*          The order of the matrices A and B. N >= 0.\n*\n*  A       (input/output) COMPLEX*16 arrays, dimensions (LDA,N)\n*          On entry, the matrix A in the pair (A, B).\n*          On exit, the updated matrix A.\n*\n*  LDA     (input)  INTEGER\n*          The leading dimension of the array A. LDA >= max(1,N).\n*\n*  B       (input/output) COMPLEX*16 arrays, dimensions (LDB,N)\n*          On entry, the matrix B in the pair (A, B).\n*          On exit, the updated matrix B.\n*\n*  LDB     (input)  INTEGER\n*          The leading dimension of the array B. LDB >= max(1,N).\n*\n*  Q       (input/output) COMPLEX*16 array, dimension (LDZ,N)\n*          If WANTQ = .TRUE, on entry, the unitary matrix Q. On exit,\n*          the updated matrix Q.\n*          Not referenced if WANTQ = .FALSE..\n*\n*  LDQ     (input) INTEGER\n*          The leading dimension of the array Q. LDQ >= 1;\n*          If WANTQ = .TRUE., LDQ >= N.\n*\n*  Z       (input/output) COMPLEX*16 array, dimension (LDZ,N)\n*          If WANTZ = .TRUE, on entry, the unitary matrix Z. On exit,\n*          the updated matrix Z.\n*          Not referenced if WANTZ = .FALSE..\n*\n*  LDZ     (input) INTEGER\n*          The leading dimension of the array Z. LDZ >= 1;\n*          If WANTZ = .TRUE., LDZ >= N.\n*\n*  J1      (input) INTEGER\n*          The index to the first block (A11, B11).\n*\n*  INFO    (output) INTEGER\n*           =0:  Successful exit.\n*           =1:  The transformed matrix pair (A, B) would be too far\n*                from generalized Schur form; the problem is ill-\n*                conditioned. \n*\n*\n\n*  Further Details\n*  ===============\n*\n*  Based on contributions by\n*     Bo Kagstrom and Peter Poromaa, Department of Computing Science,\n*     Umea University, S-901 87 Umea, Sweden.\n*\n*  In the current code both weak and strong stability tests are\n*  performed. The user can omit the strong stability test by changing\n*  the internal logical parameter WANDS to .FALSE.. See ref. [2] for\n*  details.\n*\n*  [1] B. Kagstrom; A Direct Method for Reordering Eigenvalues in the\n*      Generalized Real Schur Form of a Regular Matrix Pair (A, B), in\n*      M.S. Moonen et al (eds), Linear Algebra for Large Scale and\n*      Real-Time Applications, Kluwer Academic Publ. 1993, pp 195-218.\n*\n*  [2] B. Kagstrom and P. Poromaa; Computing Eigenspaces with Specified\n*      Eigenvalues of a Regular Matrix Pair (A, B) and Condition\n*      Estimation: Theory, Algorithms and Software, Report UMINF-94.04,\n*      Department of Computing Science, Umea University, S-901 87 Umea,\n*      Sweden, 1994. Also as LAPACK Working Note 87. To appear in\n*      Numerical Algorithms, 1996.\n*\n*  =====================================================================\n*\n\n");
      rb_exit(0);
    }
    if (rb_hash_aref(rb_options, sUsage) == Qtrue) {
      printf("%s\n", "USAGE:\n  info, a, b, q, z = NumRu::Lapack.ztgex2( wantq, wantz, a, b, q, ldq, z, ldz, j1, [:usage => usage, :help => help])\n");
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
  rblapack_ldz = argv[7];
  rblapack_j1 = argv[8];
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
  wantz = (rblapack_wantz == Qtrue);
  ldz = NUM2INT(rblapack_ldz);
  if (!NA_IsNArray(rblapack_b))
    rb_raise(rb_eArgError, "b (4th argument) must be NArray");
  if (NA_RANK(rblapack_b) != 2)
    rb_raise(rb_eArgError, "rank of b (4th argument) must be %d", 2);
  if (NA_SHAPE1(rblapack_b) != n)
    rb_raise(rb_eRuntimeError, "shape 1 of b must be the same as shape 1 of a");
  ldb = NA_SHAPE0(rblapack_b);
  if (NA_TYPE(rblapack_b) != NA_DCOMPLEX)
    rblapack_b = na_change_type(rblapack_b, NA_DCOMPLEX);
  b = NA_PTR_TYPE(rblapack_b, doublecomplex*);
  j1 = NUM2INT(rblapack_j1);
  wantq = (rblapack_wantq == Qtrue);
  ldq = NUM2INT(rblapack_ldq);
  if (!NA_IsNArray(rblapack_q))
    rb_raise(rb_eArgError, "q (5th argument) must be NArray");
  if (NA_RANK(rblapack_q) != 2)
    rb_raise(rb_eArgError, "rank of q (5th argument) must be %d", 2);
  if (NA_SHAPE1(rblapack_q) != (wantq ? n : 0))
    rb_raise(rb_eRuntimeError, "shape 1 of q must be %d", wantq ? n : 0);
  if (NA_SHAPE0(rblapack_q) != (wantq ? ldq : 0))
    rb_raise(rb_eRuntimeError, "shape 0 of q must be %d", wantq ? ldq : 0);
  if (NA_TYPE(rblapack_q) != NA_DCOMPLEX)
    rblapack_q = na_change_type(rblapack_q, NA_DCOMPLEX);
  q = NA_PTR_TYPE(rblapack_q, doublecomplex*);
  if (!NA_IsNArray(rblapack_z))
    rb_raise(rb_eArgError, "z (7th argument) must be NArray");
  if (NA_RANK(rblapack_z) != 2)
    rb_raise(rb_eArgError, "rank of z (7th argument) must be %d", 2);
  if (NA_SHAPE1(rblapack_z) != (wantq ? n : 0))
    rb_raise(rb_eRuntimeError, "shape 1 of z must be %d", wantq ? n : 0);
  if (NA_SHAPE0(rblapack_z) != (wantq ? ldz : 0))
    rb_raise(rb_eRuntimeError, "shape 0 of z must be %d", wantq ? ldz : 0);
  if (NA_TYPE(rblapack_z) != NA_DCOMPLEX)
    rblapack_z = na_change_type(rblapack_z, NA_DCOMPLEX);
  z = NA_PTR_TYPE(rblapack_z, doublecomplex*);
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
  {
    int shape[2];
    shape[0] = ldb;
    shape[1] = n;
    rblapack_b_out__ = na_make_object(NA_DCOMPLEX, 2, shape, cNArray);
  }
  b_out__ = NA_PTR_TYPE(rblapack_b_out__, doublecomplex*);
  MEMCPY(b_out__, b, doublecomplex, NA_TOTAL(rblapack_b));
  rblapack_b = rblapack_b_out__;
  b = b_out__;
  {
    int shape[2];
    shape[0] = wantq ? ldq : 0;
    shape[1] = wantq ? n : 0;
    rblapack_q_out__ = na_make_object(NA_DCOMPLEX, 2, shape, cNArray);
  }
  q_out__ = NA_PTR_TYPE(rblapack_q_out__, doublecomplex*);
  MEMCPY(q_out__, q, doublecomplex, NA_TOTAL(rblapack_q));
  rblapack_q = rblapack_q_out__;
  q = q_out__;
  {
    int shape[2];
    shape[0] = wantq ? ldz : 0;
    shape[1] = wantq ? n : 0;
    rblapack_z_out__ = na_make_object(NA_DCOMPLEX, 2, shape, cNArray);
  }
  z_out__ = NA_PTR_TYPE(rblapack_z_out__, doublecomplex*);
  MEMCPY(z_out__, z, doublecomplex, NA_TOTAL(rblapack_z));
  rblapack_z = rblapack_z_out__;
  z = z_out__;

  ztgex2_(&wantq, &wantz, &n, a, &lda, b, &ldb, q, &ldq, z, &ldz, &j1, &info);

  rblapack_info = INT2NUM(info);
  return rb_ary_new3(5, rblapack_info, rblapack_a, rblapack_b, rblapack_q, rblapack_z);
}

void
init_lapack_ztgex2(VALUE mLapack){
  rb_define_module_function(mLapack, "ztgex2", rblapack_ztgex2, -1);
  sHelp = ID2SYM(rb_intern("help"));
  sUsage = ID2SYM(rb_intern("usage"));
}
