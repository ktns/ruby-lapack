#include "rb_lapack.h"

static VALUE
rb_ztgex2(int argc, VALUE *argv, VALUE self){
  VALUE rb_wantq;
  logical wantq; 
  VALUE rb_wantz;
  logical wantz; 
  VALUE rb_a;
  doublecomplex *a; 
  VALUE rb_b;
  doublecomplex *b; 
  VALUE rb_q;
  doublecomplex *q; 
  VALUE rb_ldq;
  integer ldq; 
  VALUE rb_z;
  doublecomplex *z; 
  VALUE rb_ldz;
  integer ldz; 
  VALUE rb_j1;
  integer j1; 
  VALUE rb_info;
  integer info; 
  VALUE rb_a_out__;
  doublecomplex *a_out__;
  VALUE rb_b_out__;
  doublecomplex *b_out__;
  VALUE rb_q_out__;
  doublecomplex *q_out__;
  VALUE rb_z_out__;
  doublecomplex *z_out__;

  integer lda;
  integer n;
  integer ldb;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  info, a, b, q, z = NumRu::Lapack.ztgex2( wantq, wantz, a, b, q, ldq, z, ldz, j1)\n    or\n  NumRu::Lapack.ztgex2  # print help\n\n\nFORTRAN MANUAL\n      SUBROUTINE ZTGEX2( WANTQ, WANTZ, N, A, LDA, B, LDB, Q, LDQ, Z, LDZ, J1, INFO )\n\n*  Purpose\n*  =======\n*\n*  ZTGEX2 swaps adjacent diagonal 1 by 1 blocks (A11,B11) and (A22,B22)\n*  in an upper triangular matrix pair (A, B) by an unitary equivalence\n*  transformation.\n*\n*  (A, B) must be in generalized Schur canonical form, that is, A and\n*  B are both upper triangular.\n*\n*  Optionally, the matrices Q and Z of generalized Schur vectors are\n*  updated.\n*\n*         Q(in) * A(in) * Z(in)' = Q(out) * A(out) * Z(out)'\n*         Q(in) * B(in) * Z(in)' = Q(out) * B(out) * Z(out)'\n*\n*\n\n*  Arguments\n*  =========\n*\n*  WANTQ   (input) LOGICAL\n*          .TRUE. : update the left transformation matrix Q;\n*          .FALSE.: do not update Q.\n*\n*  WANTZ   (input) LOGICAL\n*          .TRUE. : update the right transformation matrix Z;\n*          .FALSE.: do not update Z.\n*\n*  N       (input) INTEGER\n*          The order of the matrices A and B. N >= 0.\n*\n*  A       (input/output) COMPLEX*16 arrays, dimensions (LDA,N)\n*          On entry, the matrix A in the pair (A, B).\n*          On exit, the updated matrix A.\n*\n*  LDA     (input)  INTEGER\n*          The leading dimension of the array A. LDA >= max(1,N).\n*\n*  B       (input/output) COMPLEX*16 arrays, dimensions (LDB,N)\n*          On entry, the matrix B in the pair (A, B).\n*          On exit, the updated matrix B.\n*\n*  LDB     (input)  INTEGER\n*          The leading dimension of the array B. LDB >= max(1,N).\n*\n*  Q       (input/output) COMPLEX*16 array, dimension (LDZ,N)\n*          If WANTQ = .TRUE, on entry, the unitary matrix Q. On exit,\n*          the updated matrix Q.\n*          Not referenced if WANTQ = .FALSE..\n*\n*  LDQ     (input) INTEGER\n*          The leading dimension of the array Q. LDQ >= 1;\n*          If WANTQ = .TRUE., LDQ >= N.\n*\n*  Z       (input/output) COMPLEX*16 array, dimension (LDZ,N)\n*          If WANTZ = .TRUE, on entry, the unitary matrix Z. On exit,\n*          the updated matrix Z.\n*          Not referenced if WANTZ = .FALSE..\n*\n*  LDZ     (input) INTEGER\n*          The leading dimension of the array Z. LDZ >= 1;\n*          If WANTZ = .TRUE., LDZ >= N.\n*\n*  J1      (input) INTEGER\n*          The index to the first block (A11, B11).\n*\n*  INFO    (output) INTEGER\n*           =0:  Successful exit.\n*           =1:  The transformed matrix pair (A, B) would be too far\n*                from generalized Schur form; the problem is ill-\n*                conditioned. \n*\n*\n\n*  Further Details\n*  ===============\n*\n*  Based on contributions by\n*     Bo Kagstrom and Peter Poromaa, Department of Computing Science,\n*     Umea University, S-901 87 Umea, Sweden.\n*\n*  In the current code both weak and strong stability tests are\n*  performed. The user can omit the strong stability test by changing\n*  the internal logical parameter WANDS to .FALSE.. See ref. [2] for\n*  details.\n*\n*  [1] B. Kagstrom; A Direct Method for Reordering Eigenvalues in the\n*      Generalized Real Schur Form of a Regular Matrix Pair (A, B), in\n*      M.S. Moonen et al (eds), Linear Algebra for Large Scale and\n*      Real-Time Applications, Kluwer Academic Publ. 1993, pp 195-218.\n*\n*  [2] B. Kagstrom and P. Poromaa; Computing Eigenspaces with Specified\n*      Eigenvalues of a Regular Matrix Pair (A, B) and Condition\n*      Estimation: Theory, Algorithms and Software, Report UMINF-94.04,\n*      Department of Computing Science, Umea University, S-901 87 Umea,\n*      Sweden, 1994. Also as LAPACK Working Note 87. To appear in\n*      Numerical Algorithms, 1996.\n*\n*  =====================================================================\n*\n\n");
    return Qnil;
  }
  if (argc != 9)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 9)", argc);
  rb_wantq = argv[0];
  rb_wantz = argv[1];
  rb_a = argv[2];
  rb_b = argv[3];
  rb_q = argv[4];
  rb_ldq = argv[5];
  rb_z = argv[6];
  rb_ldz = argv[7];
  rb_j1 = argv[8];

  wantq = (rb_wantq == Qtrue);
  wantz = (rb_wantz == Qtrue);
  ldq = NUM2INT(rb_ldq);
  ldz = NUM2INT(rb_ldz);
  j1 = NUM2INT(rb_j1);
  if (!NA_IsNArray(rb_a))
    rb_raise(rb_eArgError, "a (3th argument) must be NArray");
  if (NA_RANK(rb_a) != 2)
    rb_raise(rb_eArgError, "rank of a (3th argument) must be %d", 2);
  lda = NA_SHAPE0(rb_a);
  n = NA_SHAPE1(rb_a);
  if (NA_TYPE(rb_a) != NA_DCOMPLEX)
    rb_a = na_change_type(rb_a, NA_DCOMPLEX);
  a = NA_PTR_TYPE(rb_a, doublecomplex*);
  if (!NA_IsNArray(rb_b))
    rb_raise(rb_eArgError, "b (4th argument) must be NArray");
  if (NA_RANK(rb_b) != 2)
    rb_raise(rb_eArgError, "rank of b (4th argument) must be %d", 2);
  ldb = NA_SHAPE0(rb_b);
  if (NA_SHAPE1(rb_b) != n)
    rb_raise(rb_eRuntimeError, "shape 1 of b must be the same as shape 1 of a");
  if (NA_TYPE(rb_b) != NA_DCOMPLEX)
    rb_b = na_change_type(rb_b, NA_DCOMPLEX);
  b = NA_PTR_TYPE(rb_b, doublecomplex*);
  if (!NA_IsNArray(rb_q))
    rb_raise(rb_eArgError, "q (5th argument) must be NArray");
  if (NA_RANK(rb_q) != 2)
    rb_raise(rb_eArgError, "rank of q (5th argument) must be %d", 2);
  if (NA_SHAPE0(rb_q) != (wantq ? ldq : 0))
    rb_raise(rb_eRuntimeError, "shape 0 of q must be %d", wantq ? ldq : 0);
  if (NA_SHAPE1(rb_q) != (wantq ? n : 0))
    rb_raise(rb_eRuntimeError, "shape 1 of q must be %d", wantq ? n : 0);
  if (NA_TYPE(rb_q) != NA_DCOMPLEX)
    rb_q = na_change_type(rb_q, NA_DCOMPLEX);
  q = NA_PTR_TYPE(rb_q, doublecomplex*);
  if (!NA_IsNArray(rb_z))
    rb_raise(rb_eArgError, "z (7th argument) must be NArray");
  if (NA_RANK(rb_z) != 2)
    rb_raise(rb_eArgError, "rank of z (7th argument) must be %d", 2);
  if (NA_SHAPE0(rb_z) != (wantq ? ldz : 0))
    rb_raise(rb_eRuntimeError, "shape 0 of z must be %d", wantq ? ldz : 0);
  if (NA_SHAPE1(rb_z) != (wantq ? n : 0))
    rb_raise(rb_eRuntimeError, "shape 1 of z must be %d", wantq ? n : 0);
  if (NA_TYPE(rb_z) != NA_DCOMPLEX)
    rb_z = na_change_type(rb_z, NA_DCOMPLEX);
  z = NA_PTR_TYPE(rb_z, doublecomplex*);
  {
    int shape[2];
    shape[0] = lda;
    shape[1] = n;
    rb_a_out__ = na_make_object(NA_DCOMPLEX, 2, shape, cNArray);
  }
  a_out__ = NA_PTR_TYPE(rb_a_out__, doublecomplex*);
  MEMCPY(a_out__, a, doublecomplex, NA_TOTAL(rb_a));
  rb_a = rb_a_out__;
  a = a_out__;
  {
    int shape[2];
    shape[0] = ldb;
    shape[1] = n;
    rb_b_out__ = na_make_object(NA_DCOMPLEX, 2, shape, cNArray);
  }
  b_out__ = NA_PTR_TYPE(rb_b_out__, doublecomplex*);
  MEMCPY(b_out__, b, doublecomplex, NA_TOTAL(rb_b));
  rb_b = rb_b_out__;
  b = b_out__;
  {
    int shape[2];
    shape[0] = wantq ? ldq : 0;
    shape[1] = wantq ? n : 0;
    rb_q_out__ = na_make_object(NA_DCOMPLEX, 2, shape, cNArray);
  }
  q_out__ = NA_PTR_TYPE(rb_q_out__, doublecomplex*);
  MEMCPY(q_out__, q, doublecomplex, NA_TOTAL(rb_q));
  rb_q = rb_q_out__;
  q = q_out__;
  {
    int shape[2];
    shape[0] = wantq ? ldz : 0;
    shape[1] = wantq ? n : 0;
    rb_z_out__ = na_make_object(NA_DCOMPLEX, 2, shape, cNArray);
  }
  z_out__ = NA_PTR_TYPE(rb_z_out__, doublecomplex*);
  MEMCPY(z_out__, z, doublecomplex, NA_TOTAL(rb_z));
  rb_z = rb_z_out__;
  z = z_out__;

  ztgex2_(&wantq, &wantz, &n, a, &lda, b, &ldb, q, &ldq, z, &ldz, &j1, &info);

  rb_info = INT2NUM(info);
  return rb_ary_new3(5, rb_info, rb_a, rb_b, rb_q, rb_z);
}

void
init_lapack_ztgex2(VALUE mLapack){
  rb_define_module_function(mLapack, "ztgex2", rb_ztgex2, -1);
}
