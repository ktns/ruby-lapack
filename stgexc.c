#include "rb_lapack.h"

extern VOID stgexc_(logical *wantq, logical *wantz, integer *n, real *a, integer *lda, real *b, integer *ldb, real *q, integer *ldq, real *z, integer *ldz, integer *ifst, integer *ilst, real *work, integer *lwork, integer *info);

static VALUE sHelp, sUsage;

static VALUE
rblapack_stgexc(int argc, VALUE *argv, VALUE self){
  VALUE rblapack_wantq;
  logical wantq; 
  VALUE rblapack_wantz;
  logical wantz; 
  VALUE rblapack_a;
  real *a; 
  VALUE rblapack_b;
  real *b; 
  VALUE rblapack_q;
  real *q; 
  VALUE rblapack_ldq;
  integer ldq; 
  VALUE rblapack_z;
  real *z; 
  VALUE rblapack_ifst;
  integer ifst; 
  VALUE rblapack_ilst;
  integer ilst; 
  VALUE rblapack_lwork;
  integer lwork; 
  VALUE rblapack_work;
  real *work; 
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
  integer ldz;

  VALUE rb_options;
  if (argc > 0 && TYPE(argv[argc-1]) == T_HASH) {
    argc--;
    rb_options = argv[argc];
    if (rb_hash_aref(rb_options, sHelp) == Qtrue) {
      printf("%s\n", "USAGE:\n  work, info, a, b, q, z, ifst, ilst = NumRu::Lapack.stgexc( wantq, wantz, a, b, q, ldq, z, ifst, ilst, lwork, [:usage => usage, :help => help])\n\n\nFORTRAN MANUAL\n      SUBROUTINE STGEXC( WANTQ, WANTZ, N, A, LDA, B, LDB, Q, LDQ, Z, LDZ, IFST, ILST, WORK, LWORK, INFO )\n\n*  Purpose\n*  =======\n*\n*  STGEXC reorders the generalized real Schur decomposition of a real\n*  matrix pair (A,B) using an orthogonal equivalence transformation\n*\n*                 (A, B) = Q * (A, B) * Z',\n*\n*  so that the diagonal block of (A, B) with row index IFST is moved\n*  to row ILST.\n*\n*  (A, B) must be in generalized real Schur canonical form (as returned\n*  by SGGES), i.e. A is block upper triangular with 1-by-1 and 2-by-2\n*  diagonal blocks. B is upper triangular.\n*\n*  Optionally, the matrices Q and Z of generalized Schur vectors are\n*  updated.\n*\n*         Q(in) * A(in) * Z(in)' = Q(out) * A(out) * Z(out)'\n*         Q(in) * B(in) * Z(in)' = Q(out) * B(out) * Z(out)'\n*\n*\n\n*  Arguments\n*  =========\n*\n*  WANTQ   (input) LOGICAL\n*          .TRUE. : update the left transformation matrix Q;\n*          .FALSE.: do not update Q.\n*\n*  WANTZ   (input) LOGICAL\n*          .TRUE. : update the right transformation matrix Z;\n*          .FALSE.: do not update Z.\n*\n*  N       (input) INTEGER\n*          The order of the matrices A and B. N >= 0.\n*\n*  A       (input/output) REAL array, dimension (LDA,N)\n*          On entry, the matrix A in generalized real Schur canonical\n*          form.\n*          On exit, the updated matrix A, again in generalized\n*          real Schur canonical form.\n*\n*  LDA     (input)  INTEGER\n*          The leading dimension of the array A. LDA >= max(1,N).\n*\n*  B       (input/output) REAL array, dimension (LDB,N)\n*          On entry, the matrix B in generalized real Schur canonical\n*          form (A,B).\n*          On exit, the updated matrix B, again in generalized\n*          real Schur canonical form (A,B).\n*\n*  LDB     (input)  INTEGER\n*          The leading dimension of the array B. LDB >= max(1,N).\n*\n*  Q       (input/output) REAL array, dimension (LDZ,N)\n*          On entry, if WANTQ = .TRUE., the orthogonal matrix Q.\n*          On exit, the updated matrix Q.\n*          If WANTQ = .FALSE., Q is not referenced.\n*\n*  LDQ     (input) INTEGER\n*          The leading dimension of the array Q. LDQ >= 1.\n*          If WANTQ = .TRUE., LDQ >= N.\n*\n*  Z       (input/output) REAL array, dimension (LDZ,N)\n*          On entry, if WANTZ = .TRUE., the orthogonal matrix Z.\n*          On exit, the updated matrix Z.\n*          If WANTZ = .FALSE., Z is not referenced.\n*\n*  LDZ     (input) INTEGER\n*          The leading dimension of the array Z. LDZ >= 1.\n*          If WANTZ = .TRUE., LDZ >= N.\n*\n*  IFST    (input/output) INTEGER\n*  ILST    (input/output) INTEGER\n*          Specify the reordering of the diagonal blocks of (A, B).\n*          The block with row index IFST is moved to row ILST, by a\n*          sequence of swapping between adjacent blocks.\n*          On exit, if IFST pointed on entry to the second row of\n*          a 2-by-2 block, it is changed to point to the first row;\n*          ILST always points to the first row of the block in its\n*          final position (which may differ from its input value by\n*          +1 or -1). 1 <= IFST, ILST <= N.\n*\n*  WORK    (workspace/output) REAL array, dimension (MAX(1,LWORK))\n*          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.\n*\n*  LWORK   (input) INTEGER\n*          The dimension of the array WORK.\n*          LWORK >= 1 when N <= 1, otherwise LWORK >= 4*N + 16.\n*\n*          If LWORK = -1, then a workspace query is assumed; the routine\n*          only calculates the optimal size of the WORK array, returns\n*          this value as the first entry of the WORK array, and no error\n*          message related to LWORK is issued by XERBLA.\n*\n*  INFO    (output) INTEGER\n*           =0:  successful exit.\n*           <0:  if INFO = -i, the i-th argument had an illegal value.\n*           =1:  The transformed matrix pair (A, B) would be too far\n*                from generalized Schur form; the problem is ill-\n*                conditioned. (A, B) may have been partially reordered,\n*                and ILST points to the first row of the current\n*                position of the block being moved.\n*\n\n*  Further Details\n*  ===============\n*\n*  Based on contributions by\n*     Bo Kagstrom and Peter Poromaa, Department of Computing Science,\n*     Umea University, S-901 87 Umea, Sweden.\n*\n*  [1] B. Kagstrom; A Direct Method for Reordering Eigenvalues in the\n*      Generalized Real Schur Form of a Regular Matrix Pair (A, B), in\n*      M.S. Moonen et al (eds), Linear Algebra for Large Scale and\n*      Real-Time Applications, Kluwer Academic Publ. 1993, pp 195-218.\n*\n*  =====================================================================\n*\n\n");
      rb_exit(0);
    }
    if (rb_hash_aref(rb_options, sUsage) == Qtrue) {
      printf("%s\n", "USAGE:\n  work, info, a, b, q, z, ifst, ilst = NumRu::Lapack.stgexc( wantq, wantz, a, b, q, ldq, z, ifst, ilst, lwork, [:usage => usage, :help => help])\n");
      rb_exit(0);
    } 
  } else
    rb_options = Qnil;
  if (argc != 10)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 10)", argc);
  rblapack_wantq = argv[0];
  rblapack_wantz = argv[1];
  rblapack_a = argv[2];
  rblapack_b = argv[3];
  rblapack_q = argv[4];
  rblapack_ldq = argv[5];
  rblapack_z = argv[6];
  rblapack_ifst = argv[7];
  rblapack_ilst = argv[8];
  rblapack_lwork = argv[9];
  if (rb_options != Qnil) {
  }

  if (!NA_IsNArray(rblapack_a))
    rb_raise(rb_eArgError, "a (3th argument) must be NArray");
  if (NA_RANK(rblapack_a) != 2)
    rb_raise(rb_eArgError, "rank of a (3th argument) must be %d", 2);
  n = NA_SHAPE1(rblapack_a);
  lda = NA_SHAPE0(rblapack_a);
  if (NA_TYPE(rblapack_a) != NA_SFLOAT)
    rblapack_a = na_change_type(rblapack_a, NA_SFLOAT);
  a = NA_PTR_TYPE(rblapack_a, real*);
  wantz = (rblapack_wantz == Qtrue);
  if (!NA_IsNArray(rblapack_b))
    rb_raise(rb_eArgError, "b (4th argument) must be NArray");
  if (NA_RANK(rblapack_b) != 2)
    rb_raise(rb_eArgError, "rank of b (4th argument) must be %d", 2);
  if (NA_SHAPE1(rblapack_b) != n)
    rb_raise(rb_eRuntimeError, "shape 1 of b must be the same as shape 1 of a");
  ldb = NA_SHAPE0(rblapack_b);
  if (NA_TYPE(rblapack_b) != NA_SFLOAT)
    rblapack_b = na_change_type(rblapack_b, NA_SFLOAT);
  b = NA_PTR_TYPE(rblapack_b, real*);
  ldq = NUM2INT(rblapack_ldq);
  wantq = (rblapack_wantq == Qtrue);
  if (!NA_IsNArray(rblapack_z))
    rb_raise(rb_eArgError, "z (7th argument) must be NArray");
  if (NA_RANK(rblapack_z) != 2)
    rb_raise(rb_eArgError, "rank of z (7th argument) must be %d", 2);
  if (NA_SHAPE1(rblapack_z) != n)
    rb_raise(rb_eRuntimeError, "shape 1 of z must be the same as shape 1 of a");
  ldz = NA_SHAPE0(rblapack_z);
  if (NA_TYPE(rblapack_z) != NA_SFLOAT)
    rblapack_z = na_change_type(rblapack_z, NA_SFLOAT);
  z = NA_PTR_TYPE(rblapack_z, real*);
  lwork = NUM2INT(rblapack_lwork);
  ilst = NUM2INT(rblapack_ilst);
  if (!NA_IsNArray(rblapack_q))
    rb_raise(rb_eArgError, "q (5th argument) must be NArray");
  if (NA_RANK(rblapack_q) != 2)
    rb_raise(rb_eArgError, "rank of q (5th argument) must be %d", 2);
  if (NA_SHAPE1(rblapack_q) != n)
    rb_raise(rb_eRuntimeError, "shape 1 of q must be the same as shape 1 of a");
  if (NA_SHAPE0(rblapack_q) != ldz)
    rb_raise(rb_eRuntimeError, "shape 0 of q must be the same as shape 0 of z");
  if (NA_TYPE(rblapack_q) != NA_SFLOAT)
    rblapack_q = na_change_type(rblapack_q, NA_SFLOAT);
  q = NA_PTR_TYPE(rblapack_q, real*);
  ifst = NUM2INT(rblapack_ifst);
  {
    int shape[1];
    shape[0] = MAX(1,lwork);
    rblapack_work = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  work = NA_PTR_TYPE(rblapack_work, real*);
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
    shape[0] = ldz;
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

  stgexc_(&wantq, &wantz, &n, a, &lda, b, &ldb, q, &ldq, z, &ldz, &ifst, &ilst, work, &lwork, &info);

  rblapack_info = INT2NUM(info);
  rblapack_ifst = INT2NUM(ifst);
  rblapack_ilst = INT2NUM(ilst);
  return rb_ary_new3(8, rblapack_work, rblapack_info, rblapack_a, rblapack_b, rblapack_q, rblapack_z, rblapack_ifst, rblapack_ilst);
}

void
init_lapack_stgexc(VALUE mLapack){
  rb_define_module_function(mLapack, "stgexc", rblapack_stgexc, -1);
  sHelp = ID2SYM(rb_intern("help"));
  sUsage = ID2SYM(rb_intern("usage"));
}
