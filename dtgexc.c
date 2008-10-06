#include "rb_lapack.h"

static VALUE
rb_dtgexc(int argc, VALUE *argv, VALUE self){
  VALUE rb_wantq;
  logical wantq; 
  VALUE rb_wantz;
  logical wantz; 
  VALUE rb_a;
  doublereal *a; 
  VALUE rb_b;
  doublereal *b; 
  VALUE rb_q;
  doublereal *q; 
  VALUE rb_ldq;
  integer ldq; 
  VALUE rb_z;
  doublereal *z; 
  VALUE rb_ifst;
  integer ifst; 
  VALUE rb_ilst;
  integer ilst; 
  VALUE rb_lwork;
  integer lwork; 
  VALUE rb_work;
  doublereal *work; 
  VALUE rb_info;
  integer info; 
  VALUE rb_a_out__;
  doublereal *a_out__;
  VALUE rb_b_out__;
  doublereal *b_out__;
  VALUE rb_q_out__;
  doublereal *q_out__;
  VALUE rb_z_out__;
  doublereal *z_out__;

  integer lda;
  integer n;
  integer ldb;
  integer ldz;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  work, info, a, b, q, z, ifst, ilst = NumRu::Lapack.dtgexc( wantq, wantz, a, b, q, ldq, z, ifst, ilst, lwork)\n    or\n  NumRu::Lapack.dtgexc  # print help\n\n\nFORTRAN MANUAL\n      SUBROUTINE DTGEXC( WANTQ, WANTZ, N, A, LDA, B, LDB, Q, LDQ, Z, LDZ, IFST, ILST, WORK, LWORK, INFO )\n\n*  Purpose\n*  =======\n*\n*  DTGEXC reorders the generalized real Schur decomposition of a real\n*  matrix pair (A,B) using an orthogonal equivalence transformation\n*\n*                 (A, B) = Q * (A, B) * Z',\n*\n*  so that the diagonal block of (A, B) with row index IFST is moved\n*  to row ILST.\n*\n*  (A, B) must be in generalized real Schur canonical form (as returned\n*  by DGGES), i.e. A is block upper triangular with 1-by-1 and 2-by-2\n*  diagonal blocks. B is upper triangular.\n*\n*  Optionally, the matrices Q and Z of generalized Schur vectors are\n*  updated.\n*\n*         Q(in) * A(in) * Z(in)' = Q(out) * A(out) * Z(out)'\n*         Q(in) * B(in) * Z(in)' = Q(out) * B(out) * Z(out)'\n*\n*\n\n*  Arguments\n*  =========\n*\n*  WANTQ   (input) LOGICAL\n*          .TRUE. : update the left transformation matrix Q;\n*          .FALSE.: do not update Q.\n*\n*  WANTZ   (input) LOGICAL\n*          .TRUE. : update the right transformation matrix Z;\n*          .FALSE.: do not update Z.\n*\n*  N       (input) INTEGER\n*          The order of the matrices A and B. N >= 0.\n*\n*  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)\n*          On entry, the matrix A in generalized real Schur canonical\n*          form.\n*          On exit, the updated matrix A, again in generalized\n*          real Schur canonical form.\n*\n*  LDA     (input)  INTEGER\n*          The leading dimension of the array A. LDA >= max(1,N).\n*\n*  B       (input/output) DOUBLE PRECISION array, dimension (LDB,N)\n*          On entry, the matrix B in generalized real Schur canonical\n*          form (A,B).\n*          On exit, the updated matrix B, again in generalized\n*          real Schur canonical form (A,B).\n*\n*  LDB     (input)  INTEGER\n*          The leading dimension of the array B. LDB >= max(1,N).\n*\n*  Q       (input/output) DOUBLE PRECISION array, dimension (LDZ,N)\n*          On entry, if WANTQ = .TRUE., the orthogonal matrix Q.\n*          On exit, the updated matrix Q.\n*          If WANTQ = .FALSE., Q is not referenced.\n*\n*  LDQ     (input) INTEGER\n*          The leading dimension of the array Q. LDQ >= 1.\n*          If WANTQ = .TRUE., LDQ >= N.\n*\n*  Z       (input/output) DOUBLE PRECISION array, dimension (LDZ,N)\n*          On entry, if WANTZ = .TRUE., the orthogonal matrix Z.\n*          On exit, the updated matrix Z.\n*          If WANTZ = .FALSE., Z is not referenced.\n*\n*  LDZ     (input) INTEGER\n*          The leading dimension of the array Z. LDZ >= 1.\n*          If WANTZ = .TRUE., LDZ >= N.\n*\n*  IFST    (input/output) INTEGER\n*  ILST    (input/output) INTEGER\n*          Specify the reordering of the diagonal blocks of (A, B).\n*          The block with row index IFST is moved to row ILST, by a\n*          sequence of swapping between adjacent blocks.\n*          On exit, if IFST pointed on entry to the second row of\n*          a 2-by-2 block, it is changed to point to the first row;\n*          ILST always points to the first row of the block in its\n*          final position (which may differ from its input value by\n*          +1 or -1). 1 <= IFST, ILST <= N.\n*\n*  WORK    (workspace/output) DOUBLE PRECISION array, dimension (MAX(1,LWORK))\n*          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.\n*\n*  LWORK   (input) INTEGER\n*          The dimension of the array WORK.\n*          LWORK >= 1 when N <= 1, otherwise LWORK >= 4*N + 16.\n*\n*          If LWORK = -1, then a workspace query is assumed; the routine\n*          only calculates the optimal size of the WORK array, returns\n*          this value as the first entry of the WORK array, and no error\n*          message related to LWORK is issued by XERBLA.\n*\n*  INFO    (output) INTEGER\n*           =0:  successful exit.\n*           <0:  if INFO = -i, the i-th argument had an illegal value.\n*           =1:  The transformed matrix pair (A, B) would be too far\n*                from generalized Schur form; the problem is ill-\n*                conditioned. (A, B) may have been partially reordered,\n*                and ILST points to the first row of the current\n*                position of the block being moved.\n*\n\n*  Further Details\n*  ===============\n*\n*  Based on contributions by\n*     Bo Kagstrom and Peter Poromaa, Department of Computing Science,\n*     Umea University, S-901 87 Umea, Sweden.\n*\n*  [1] B. Kagstrom; A Direct Method for Reordering Eigenvalues in the\n*      Generalized Real Schur Form of a Regular Matrix Pair (A, B), in\n*      M.S. Moonen et al (eds), Linear Algebra for Large Scale and\n*      Real-Time Applications, Kluwer Academic Publ. 1993, pp 195-218.\n*\n*  =====================================================================\n*\n\n");
    return Qnil;
  }
  if (argc != 10)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 10)", argc);
  rb_wantq = argv[0];
  rb_wantz = argv[1];
  rb_a = argv[2];
  rb_b = argv[3];
  rb_q = argv[4];
  rb_ldq = argv[5];
  rb_z = argv[6];
  rb_ifst = argv[7];
  rb_ilst = argv[8];
  rb_lwork = argv[9];

  wantq = (rb_wantq == Qtrue);
  wantz = (rb_wantz == Qtrue);
  ldq = NUM2INT(rb_ldq);
  ifst = NUM2INT(rb_ifst);
  ilst = NUM2INT(rb_ilst);
  lwork = NUM2INT(rb_lwork);
  if (!NA_IsNArray(rb_a))
    rb_raise(rb_eArgError, "a (3th argument) must be NArray");
  if (NA_RANK(rb_a) != 2)
    rb_raise(rb_eArgError, "rank of a (3th argument) must be %d", 2);
  lda = NA_SHAPE0(rb_a);
  n = NA_SHAPE1(rb_a);
  if (NA_TYPE(rb_a) != NA_DFLOAT)
    rb_a = na_change_type(rb_a, NA_DFLOAT);
  a = NA_PTR_TYPE(rb_a, doublereal*);
  if (!NA_IsNArray(rb_b))
    rb_raise(rb_eArgError, "b (4th argument) must be NArray");
  if (NA_RANK(rb_b) != 2)
    rb_raise(rb_eArgError, "rank of b (4th argument) must be %d", 2);
  ldb = NA_SHAPE0(rb_b);
  if (NA_SHAPE1(rb_b) != n)
    rb_raise(rb_eRuntimeError, "shape 1 of b must be the same as shape 1 of a");
  if (NA_TYPE(rb_b) != NA_DFLOAT)
    rb_b = na_change_type(rb_b, NA_DFLOAT);
  b = NA_PTR_TYPE(rb_b, doublereal*);
  if (!NA_IsNArray(rb_q))
    rb_raise(rb_eArgError, "q (5th argument) must be NArray");
  if (NA_RANK(rb_q) != 2)
    rb_raise(rb_eArgError, "rank of q (5th argument) must be %d", 2);
  ldz = NA_SHAPE0(rb_q);
  if (NA_SHAPE1(rb_q) != n)
    rb_raise(rb_eRuntimeError, "shape 1 of q must be the same as shape 1 of a");
  if (NA_TYPE(rb_q) != NA_DFLOAT)
    rb_q = na_change_type(rb_q, NA_DFLOAT);
  q = NA_PTR_TYPE(rb_q, doublereal*);
  if (!NA_IsNArray(rb_z))
    rb_raise(rb_eArgError, "z (7th argument) must be NArray");
  if (NA_RANK(rb_z) != 2)
    rb_raise(rb_eArgError, "rank of z (7th argument) must be %d", 2);
  if (NA_SHAPE0(rb_z) != ldz)
    rb_raise(rb_eRuntimeError, "shape 0 of z must be the same as shape 0 of q");
  if (NA_SHAPE1(rb_z) != n)
    rb_raise(rb_eRuntimeError, "shape 1 of z must be the same as shape 1 of a");
  if (NA_TYPE(rb_z) != NA_DFLOAT)
    rb_z = na_change_type(rb_z, NA_DFLOAT);
  z = NA_PTR_TYPE(rb_z, doublereal*);
  {
    int shape[1];
    shape[0] = MAX(1,lwork);
    rb_work = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  work = NA_PTR_TYPE(rb_work, doublereal*);
  {
    int shape[2];
    shape[0] = lda;
    shape[1] = n;
    rb_a_out__ = na_make_object(NA_DFLOAT, 2, shape, cNArray);
  }
  a_out__ = NA_PTR_TYPE(rb_a_out__, doublereal*);
  MEMCPY(a_out__, a, doublereal, NA_TOTAL(rb_a));
  rb_a = rb_a_out__;
  a = a_out__;
  {
    int shape[2];
    shape[0] = ldb;
    shape[1] = n;
    rb_b_out__ = na_make_object(NA_DFLOAT, 2, shape, cNArray);
  }
  b_out__ = NA_PTR_TYPE(rb_b_out__, doublereal*);
  MEMCPY(b_out__, b, doublereal, NA_TOTAL(rb_b));
  rb_b = rb_b_out__;
  b = b_out__;
  {
    int shape[2];
    shape[0] = ldz;
    shape[1] = n;
    rb_q_out__ = na_make_object(NA_DFLOAT, 2, shape, cNArray);
  }
  q_out__ = NA_PTR_TYPE(rb_q_out__, doublereal*);
  MEMCPY(q_out__, q, doublereal, NA_TOTAL(rb_q));
  rb_q = rb_q_out__;
  q = q_out__;
  {
    int shape[2];
    shape[0] = ldz;
    shape[1] = n;
    rb_z_out__ = na_make_object(NA_DFLOAT, 2, shape, cNArray);
  }
  z_out__ = NA_PTR_TYPE(rb_z_out__, doublereal*);
  MEMCPY(z_out__, z, doublereal, NA_TOTAL(rb_z));
  rb_z = rb_z_out__;
  z = z_out__;

  dtgexc_(&wantq, &wantz, &n, a, &lda, b, &ldb, q, &ldq, z, &ldz, &ifst, &ilst, work, &lwork, &info);

  rb_info = INT2NUM(info);
  rb_ifst = INT2NUM(ifst);
  rb_ilst = INT2NUM(ilst);
  return rb_ary_new3(8, rb_work, rb_info, rb_a, rb_b, rb_q, rb_z, rb_ifst, rb_ilst);
}

void
init_lapack_dtgexc(VALUE mLapack){
  rb_define_module_function(mLapack, "dtgexc", rb_dtgexc, -1);
}
