#include "rb_lapack.h"

extern doublereal zlanhf_(char *norm, char *transr, char *uplo, integer *n, doublecomplex *a, doublereal *work);

static VALUE sHelp, sUsage;

static VALUE
rblapack_zlanhf(int argc, VALUE *argv, VALUE self){
  VALUE rblapack_norm;
  char norm; 
  VALUE rblapack_transr;
  char transr; 
  VALUE rblapack_uplo;
  char uplo; 
  VALUE rblapack_n;
  integer n; 
  VALUE rblapack_a;
  doublecomplex *a; 
  VALUE rblapack___out__;
  doublereal __out__; 
  doublereal *work;

  integer lwork;

  VALUE rb_options;
  if (argc > 0 && TYPE(argv[argc-1]) == T_HASH) {
    argc--;
    rb_options = argv[argc];
    if (rb_hash_aref(rb_options, sHelp) == Qtrue) {
      printf("%s\n", "USAGE:\n  __out__ = NumRu::Lapack.zlanhf( norm, transr, uplo, n, a, [:usage => usage, :help => help])\n\n\nFORTRAN MANUAL\n      DOUBLE PRECISION FUNCTION ZLANHF( NORM, TRANSR, UPLO, N, A, WORK )\n\n*  Purpose\n*  =======\n*\n*  ZLANHF  returns the value of the one norm,  or the Frobenius norm, or\n*  the  infinity norm,  or the  element of  largest absolute value  of a\n*  complex Hermitian matrix A in RFP format.\n*\n*  Description\n*  ===========\n*\n*  ZLANHF returns the value\n*\n*     ZLANHF = ( max(abs(A(i,j))), NORM = 'M' or 'm'\n*              (\n*              ( norm1(A),         NORM = '1', 'O' or 'o'\n*              (\n*              ( normI(A),         NORM = 'I' or 'i'\n*              (\n*              ( normF(A),         NORM = 'F', 'f', 'E' or 'e'\n*\n*  where  norm1  denotes the  one norm of a matrix (maximum column sum),\n*  normI  denotes the  infinity norm  of a matrix  (maximum row sum) and\n*  normF  denotes the  Frobenius norm of a matrix (square root of sum of\n*  squares).  Note that  max(abs(A(i,j)))  is not a  matrix norm.\n*\n\n*  Arguments\n*  =========\n*\n*  NORM      (input) CHARACTER\n*            Specifies the value to be returned in ZLANHF as described\n*            above.\n*\n*  TRANSR    (input) CHARACTER\n*            Specifies whether the RFP format of A is normal or\n*            conjugate-transposed format.\n*            = 'N':  RFP format is Normal\n*            = 'C':  RFP format is Conjugate-transposed\n*\n*  UPLO      (input) CHARACTER\n*            On entry, UPLO specifies whether the RFP matrix A came from\n*            an upper or lower triangular matrix as follows:\n*\n*            UPLO = 'U' or 'u' RFP A came from an upper triangular\n*            matrix\n*\n*            UPLO = 'L' or 'l' RFP A came from a  lower triangular\n*            matrix\n*\n*  N         (input) INTEGER\n*            The order of the matrix A.  N >= 0.  When N = 0, ZLANHF is\n*            set to zero.\n*\n*   A        (input) COMPLEX*16 array, dimension ( N*(N+1)/2 );\n*            On entry, the matrix A in RFP Format.\n*            RFP Format is described by TRANSR, UPLO and N as follows:\n*            If TRANSR='N' then RFP A is (0:N,0:K-1) when N is even;\n*            K=N/2. RFP A is (0:N-1,0:K) when N is odd; K=N/2. If\n*            TRANSR = 'C' then RFP is the Conjugate-transpose of RFP A\n*            as defined when TRANSR = 'N'. The contents of RFP A are\n*            defined by UPLO as follows: If UPLO = 'U' the RFP A\n*            contains the ( N*(N+1)/2 ) elements of upper packed A\n*            either in normal or conjugate-transpose Format. If\n*            UPLO = 'L' the RFP A contains the ( N*(N+1) /2 ) elements\n*            of lower packed A either in normal or conjugate-transpose\n*            Format. The LDA of RFP A is (N+1)/2 when TRANSR = 'C'. When\n*            TRANSR is 'N' the LDA is N+1 when N is even and is N when\n*            is odd. See the Note below for more details.\n*            Unchanged on exit.\n*\n*  WORK      (workspace) DOUBLE PRECISION array, dimension (LWORK),\n*            where LWORK >= N when NORM = 'I' or '1' or 'O'; otherwise,\n*            WORK is not referenced.\n*\n\n*  Further Details\n*  ===============\n*\n*  We first consider Standard Packed Format when N is even.\n*  We give an example where N = 6.\n*\n*      AP is Upper             AP is Lower\n*\n*   00 01 02 03 04 05       00\n*      11 12 13 14 15       10 11\n*         22 23 24 25       20 21 22\n*            33 34 35       30 31 32 33\n*               44 45       40 41 42 43 44\n*                  55       50 51 52 53 54 55\n*\n*\n*  Let TRANSR = 'N'. RFP holds AP as follows:\n*  For UPLO = 'U' the upper trapezoid A(0:5,0:2) consists of the last\n*  three columns of AP upper. The lower triangle A(4:6,0:2) consists of\n*  conjugate-transpose of the first three columns of AP upper.\n*  For UPLO = 'L' the lower trapezoid A(1:6,0:2) consists of the first\n*  three columns of AP lower. The upper triangle A(0:2,0:2) consists of\n*  conjugate-transpose of the last three columns of AP lower.\n*  To denote conjugate we place -- above the element. This covers the\n*  case N even and TRANSR = 'N'.\n*\n*         RFP A                   RFP A\n*\n*                                -- -- --\n*        03 04 05                33 43 53\n*                                   -- --\n*        13 14 15                00 44 54\n*                                      --\n*        23 24 25                10 11 55\n*\n*        33 34 35                20 21 22\n*        --\n*        00 44 45                30 31 32\n*        -- --\n*        01 11 55                40 41 42\n*        -- -- --\n*        02 12 22                50 51 52\n*\n*  Now let TRANSR = 'C'. RFP A in both UPLO cases is just the conjugate-\n*  transpose of RFP A above. One therefore gets:\n*\n*\n*           RFP A                   RFP A\n*\n*     -- -- -- --                -- -- -- -- -- --\n*     03 13 23 33 00 01 02    33 00 10 20 30 40 50\n*     -- -- -- -- --                -- -- -- -- --\n*     04 14 24 34 44 11 12    43 44 11 21 31 41 51\n*     -- -- -- -- -- --                -- -- -- --\n*     05 15 25 35 45 55 22    53 54 55 22 32 42 52\n*\n*\n*  We next  consider Standard Packed Format when N is odd.\n*  We give an example where N = 5.\n*\n*     AP is Upper                 AP is Lower\n*\n*   00 01 02 03 04              00\n*      11 12 13 14              10 11\n*         22 23 24              20 21 22\n*            33 34              30 31 32 33\n*               44              40 41 42 43 44\n*\n*\n*  Let TRANSR = 'N'. RFP holds AP as follows:\n*  For UPLO = 'U' the upper trapezoid A(0:4,0:2) consists of the last\n*  three columns of AP upper. The lower triangle A(3:4,0:1) consists of\n*  conjugate-transpose of the first two   columns of AP upper.\n*  For UPLO = 'L' the lower trapezoid A(0:4,0:2) consists of the first\n*  three columns of AP lower. The upper triangle A(0:1,1:2) consists of\n*  conjugate-transpose of the last two   columns of AP lower.\n*  To denote conjugate we place -- above the element. This covers the\n*  case N odd  and TRANSR = 'N'.\n*\n*         RFP A                   RFP A\n*\n*                                   -- --\n*        02 03 04                00 33 43\n*                                      --\n*        12 13 14                10 11 44\n*\n*        22 23 24                20 21 22\n*        --\n*        00 33 34                30 31 32\n*        -- --\n*        01 11 44                40 41 42\n*\n*  Now let TRANSR = 'C'. RFP A in both UPLO cases is just the conjugate-\n*  transpose of RFP A above. One therefore gets:\n*\n*\n*           RFP A                   RFP A\n*\n*     -- -- --                   -- -- -- -- -- --\n*     02 12 22 00 01             00 10 20 30 40 50\n*     -- -- -- --                   -- -- -- -- --\n*     03 13 23 33 11             33 11 21 31 41 51\n*     -- -- -- -- --                   -- -- -- --\n*     04 14 24 34 44             43 44 22 32 42 52\n*\n*  =====================================================================\n*\n\n");
      rb_exit(0);
    }
    if (rb_hash_aref(rb_options, sUsage) == Qtrue) {
      printf("%s\n", "USAGE:\n  __out__ = NumRu::Lapack.zlanhf( norm, transr, uplo, n, a, [:usage => usage, :help => help])\n");
      rb_exit(0);
    } 
  } else
    rb_options = Qnil;
  if (argc != 5)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 5)", argc);
  rblapack_norm = argv[0];
  rblapack_transr = argv[1];
  rblapack_uplo = argv[2];
  rblapack_n = argv[3];
  rblapack_a = argv[4];
  if (rb_options != Qnil) {
  }

  transr = StringValueCStr(rblapack_transr)[0];
  n = NUM2INT(rblapack_n);
  uplo = StringValueCStr(rblapack_uplo)[0];
  norm = StringValueCStr(rblapack_norm)[0];
  lwork = ((lsame_(&norm,"I")) || ((('1') || ('o')))) ? n : 0;
  if (!NA_IsNArray(rblapack_a))
    rb_raise(rb_eArgError, "a (5th argument) must be NArray");
  if (NA_RANK(rblapack_a) != 1)
    rb_raise(rb_eArgError, "rank of a (5th argument) must be %d", 1);
  if (NA_SHAPE0(rblapack_a) != (n*(n+1)/2))
    rb_raise(rb_eRuntimeError, "shape 0 of a must be %d", n*(n+1)/2);
  if (NA_TYPE(rblapack_a) != NA_DCOMPLEX)
    rblapack_a = na_change_type(rblapack_a, NA_DCOMPLEX);
  a = NA_PTR_TYPE(rblapack_a, doublecomplex*);
  work = ALLOC_N(doublereal, (lwork));

  __out__ = zlanhf_(&norm, &transr, &uplo, &n, a, work);

  free(work);
  rblapack___out__ = rb_float_new((double)__out__);
  return rblapack___out__;
}

void
init_lapack_zlanhf(VALUE mLapack){
  rb_define_module_function(mLapack, "zlanhf", rblapack_zlanhf, -1);
  sHelp = ID2SYM(rb_intern("help"));
  sUsage = ID2SYM(rb_intern("usage"));
}
