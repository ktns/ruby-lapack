#include "rb_lapack.h"

static VALUE
rb_dgghrd(int argc, VALUE *argv, VALUE self){
  VALUE rb_compq;
  char compq; 
  VALUE rb_compz;
  char compz; 
  VALUE rb_ilo;
  integer ilo; 
  VALUE rb_ihi;
  integer ihi; 
  VALUE rb_a;
  doublereal *a; 
  VALUE rb_b;
  doublereal *b; 
  VALUE rb_q;
  doublereal *q; 
  VALUE rb_z;
  doublereal *z; 
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
  integer ldq;
  integer ldz;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  info, a, b, q, z = NumRu::Lapack.dgghrd( compq, compz, ilo, ihi, a, b, q, z)\n    or\n  NumRu::Lapack.dgghrd  # print help\n\n\nFORTRAN MANUAL\n      SUBROUTINE DGGHRD( COMPQ, COMPZ, N, ILO, IHI, A, LDA, B, LDB, Q, LDQ, Z, LDZ, INFO )\n\n*  Purpose\n*  =======\n*\n*  DGGHRD reduces a pair of real matrices (A,B) to generalized upper\n*  Hessenberg form using orthogonal transformations, where A is a\n*  general matrix and B is upper triangular.  The form of the\n*  generalized eigenvalue problem is\n*     A*x = lambda*B*x,\n*  and B is typically made upper triangular by computing its QR\n*  factorization and moving the orthogonal matrix Q to the left side\n*  of the equation.\n*\n*  This subroutine simultaneously reduces A to a Hessenberg matrix H:\n*     Q**T*A*Z = H\n*  and transforms B to another upper triangular matrix T:\n*     Q**T*B*Z = T\n*  in order to reduce the problem to its standard form\n*     H*y = lambda*T*y\n*  where y = Z**T*x.\n*\n*  The orthogonal matrices Q and Z are determined as products of Givens\n*  rotations.  They may either be formed explicitly, or they may be\n*  postmultiplied into input matrices Q1 and Z1, so that\n*\n*       Q1 * A * Z1**T = (Q1*Q) * H * (Z1*Z)**T\n*\n*       Q1 * B * Z1**T = (Q1*Q) * T * (Z1*Z)**T\n*\n*  If Q1 is the orthogonal matrix from the QR factorization of B in the\n*  original equation A*x = lambda*B*x, then DGGHRD reduces the original\n*  problem to generalized Hessenberg form.\n*\n\n*  Arguments\n*  =========\n*\n*  COMPQ   (input) CHARACTER*1\n*          = 'N': do not compute Q;\n*          = 'I': Q is initialized to the unit matrix, and the\n*                 orthogonal matrix Q is returned;\n*          = 'V': Q must contain an orthogonal matrix Q1 on entry,\n*                 and the product Q1*Q is returned.\n*\n*  COMPZ   (input) CHARACTER*1\n*          = 'N': do not compute Z;\n*          = 'I': Z is initialized to the unit matrix, and the\n*                 orthogonal matrix Z is returned;\n*          = 'V': Z must contain an orthogonal matrix Z1 on entry,\n*                 and the product Z1*Z is returned.\n*\n*  N       (input) INTEGER\n*          The order of the matrices A and B.  N >= 0.\n*\n*  ILO     (input) INTEGER\n*  IHI     (input) INTEGER\n*          ILO and IHI mark the rows and columns of A which are to be\n*          reduced.  It is assumed that A is already upper triangular\n*          in rows and columns 1:ILO-1 and IHI+1:N.  ILO and IHI are\n*          normally set by a previous call to SGGBAL; otherwise they\n*          should be set to 1 and N respectively.\n*          1 <= ILO <= IHI <= N, if N > 0; ILO=1 and IHI=0, if N=0.\n*\n*  A       (input/output) DOUBLE PRECISION array, dimension (LDA, N)\n*          On entry, the N-by-N general matrix to be reduced.\n*          On exit, the upper triangle and the first subdiagonal of A\n*          are overwritten with the upper Hessenberg matrix H, and the\n*          rest is set to zero.\n*\n*  LDA     (input) INTEGER\n*          The leading dimension of the array A.  LDA >= max(1,N).\n*\n*  B       (input/output) DOUBLE PRECISION array, dimension (LDB, N)\n*          On entry, the N-by-N upper triangular matrix B.\n*          On exit, the upper triangular matrix T = Q**T B Z.  The\n*          elements below the diagonal are set to zero.\n*\n*  LDB     (input) INTEGER\n*          The leading dimension of the array B.  LDB >= max(1,N).\n*\n*  Q       (input/output) DOUBLE PRECISION array, dimension (LDQ, N)\n*          On entry, if COMPQ = 'V', the orthogonal matrix Q1,\n*          typically from the QR factorization of B.\n*          On exit, if COMPQ='I', the orthogonal matrix Q, and if\n*          COMPQ = 'V', the product Q1*Q.\n*          Not referenced if COMPQ='N'.\n*\n*  LDQ     (input) INTEGER\n*          The leading dimension of the array Q.\n*          LDQ >= N if COMPQ='V' or 'I'; LDQ >= 1 otherwise.\n*\n*  Z       (input/output) DOUBLE PRECISION array, dimension (LDZ, N)\n*          On entry, if COMPZ = 'V', the orthogonal matrix Z1.\n*          On exit, if COMPZ='I', the orthogonal matrix Z, and if\n*          COMPZ = 'V', the product Z1*Z.\n*          Not referenced if COMPZ='N'.\n*\n*  LDZ     (input) INTEGER\n*          The leading dimension of the array Z.\n*          LDZ >= N if COMPZ='V' or 'I'; LDZ >= 1 otherwise.\n*\n*  INFO    (output) INTEGER\n*          = 0:  successful exit.\n*          < 0:  if INFO = -i, the i-th argument had an illegal value.\n*\n\n*  Further Details\n*  ===============\n*\n*  This routine reduces A to Hessenberg and B to triangular form by\n*  an unblocked reduction, as described in _Matrix_Computations_,\n*  by Golub and Van Loan (Johns Hopkins Press.)\n*\n*  =====================================================================\n*\n\n");
    return Qnil;
  }
  if (argc != 8)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 8)", argc);
  rb_compq = argv[0];
  rb_compz = argv[1];
  rb_ilo = argv[2];
  rb_ihi = argv[3];
  rb_a = argv[4];
  rb_b = argv[5];
  rb_q = argv[6];
  rb_z = argv[7];

  compq = StringValueCStr(rb_compq)[0];
  compz = StringValueCStr(rb_compz)[0];
  ilo = NUM2INT(rb_ilo);
  ihi = NUM2INT(rb_ihi);
  if (!NA_IsNArray(rb_a))
    rb_raise(rb_eArgError, "a (5th argument) must be NArray");
  if (NA_RANK(rb_a) != 2)
    rb_raise(rb_eArgError, "rank of a (5th argument) must be %d", 2);
  lda = NA_SHAPE0(rb_a);
  n = NA_SHAPE1(rb_a);
  if (NA_TYPE(rb_a) != NA_DFLOAT)
    rb_a = na_change_type(rb_a, NA_DFLOAT);
  a = NA_PTR_TYPE(rb_a, doublereal*);
  if (!NA_IsNArray(rb_b))
    rb_raise(rb_eArgError, "b (6th argument) must be NArray");
  if (NA_RANK(rb_b) != 2)
    rb_raise(rb_eArgError, "rank of b (6th argument) must be %d", 2);
  ldb = NA_SHAPE0(rb_b);
  if (NA_SHAPE1(rb_b) != n)
    rb_raise(rb_eRuntimeError, "shape 1 of b must be the same as shape 1 of a");
  if (NA_TYPE(rb_b) != NA_DFLOAT)
    rb_b = na_change_type(rb_b, NA_DFLOAT);
  b = NA_PTR_TYPE(rb_b, doublereal*);
  if (!NA_IsNArray(rb_q))
    rb_raise(rb_eArgError, "q (7th argument) must be NArray");
  if (NA_RANK(rb_q) != 2)
    rb_raise(rb_eArgError, "rank of q (7th argument) must be %d", 2);
  ldq = NA_SHAPE0(rb_q);
  if (NA_SHAPE1(rb_q) != n)
    rb_raise(rb_eRuntimeError, "shape 1 of q must be the same as shape 1 of a");
  if (NA_TYPE(rb_q) != NA_DFLOAT)
    rb_q = na_change_type(rb_q, NA_DFLOAT);
  q = NA_PTR_TYPE(rb_q, doublereal*);
  if (!NA_IsNArray(rb_z))
    rb_raise(rb_eArgError, "z (8th argument) must be NArray");
  if (NA_RANK(rb_z) != 2)
    rb_raise(rb_eArgError, "rank of z (8th argument) must be %d", 2);
  ldz = NA_SHAPE0(rb_z);
  if (NA_SHAPE1(rb_z) != n)
    rb_raise(rb_eRuntimeError, "shape 1 of z must be the same as shape 1 of a");
  if (NA_TYPE(rb_z) != NA_DFLOAT)
    rb_z = na_change_type(rb_z, NA_DFLOAT);
  z = NA_PTR_TYPE(rb_z, doublereal*);
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
    shape[0] = ldq;
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

  dgghrd_(&compq, &compz, &n, &ilo, &ihi, a, &lda, b, &ldb, q, &ldq, z, &ldz, &info);

  rb_info = INT2NUM(info);
  return rb_ary_new3(5, rb_info, rb_a, rb_b, rb_q, rb_z);
}

void
init_lapack_dgghrd(VALUE mLapack){
  rb_define_module_function(mLapack, "dgghrd", rb_dgghrd, -1);
}
