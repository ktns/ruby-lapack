#include "rb_lapack.h"

extern real slansf_(char *norm, char *transr, char *uplo, integer *n, real *a, real *work);

static VALUE
rb_slansf(int argc, VALUE *argv, VALUE self){
  VALUE rb_norm;
  char norm; 
  VALUE rb_transr;
  char transr; 
  VALUE rb_uplo;
  char uplo; 
  VALUE rb_n;
  integer n; 
  VALUE rb_a;
  real *a; 
  VALUE rb___out__;
  real __out__; 
  real *work;


  if (argc == 0) {
    printf("%s\n", "USAGE:\n  __out__ = NumRu::Lapack.slansf( norm, transr, uplo, n, a)\n    or\n  NumRu::Lapack.slansf  # print help\n\n\nFORTRAN MANUAL\n      REAL FUNCTION SLANSF( NORM, TRANSR, UPLO, N, A, WORK )\n\n*  Purpose\n*  =======\n*\n*  SLANSF returns the value of the one norm, or the Frobenius norm, or\n*  the infinity norm, or the element of largest absolute value of a\n*  real symmetric matrix A in RFP format.\n*\n*  Description\n*  ===========\n*\n*  SLANSF returns the value\n*\n*     SLANSF = ( max(abs(A(i,j))), NORM = 'M' or 'm'\n*              (\n*              ( norm1(A),         NORM = '1', 'O' or 'o'\n*              (\n*              ( normI(A),         NORM = 'I' or 'i'\n*              (\n*              ( normF(A),         NORM = 'F', 'f', 'E' or 'e'\n*\n*  where  norm1  denotes the  one norm of a matrix (maximum column sum),\n*  normI  denotes the  infinity norm  of a matrix  (maximum row sum) and\n*  normF  denotes the  Frobenius norm of a matrix (square root of sum of\n*  squares).  Note that  max(abs(A(i,j)))  is not a  matrix norm.\n*\n\n*  Arguments\n*  =========\n*\n*  NORM    (input) CHARACTER*1\n*          Specifies the value to be returned in SLANSF as described\n*          above.\n*\n*  TRANSR  (input) CHARACTER*1\n*          Specifies whether the RFP format of A is normal or\n*          transposed format.\n*          = 'N':  RFP format is Normal;\n*          = 'T':  RFP format is Transpose.\n*\n*  UPLO    (input) CHARACTER*1\n*           On entry, UPLO specifies whether the RFP matrix A came from\n*           an upper or lower triangular matrix as follows:\n*           = 'U': RFP A came from an upper triangular matrix;\n*           = 'L': RFP A came from a lower triangular matrix.\n*\n*  N       (input) INTEGER\n*          The order of the matrix A. N >= 0. When N = 0, SLANSF is\n*          set to zero.\n*\n*  A       (input) REAL array, dimension ( N*(N+1)/2 );\n*          On entry, the upper (if UPLO = 'U') or lower (if UPLO = 'L')\n*          part of the symmetric matrix A stored in RFP format. See the\n*          \"Notes\" below for more details.\n*          Unchanged on exit.\n*\n*  WORK    (workspace) REAL array, dimension (MAX(1,LWORK)),\n*          where LWORK >= N when NORM = 'I' or '1' or 'O'; otherwise,\n*          WORK is not referenced.\n*\n\n*  Further Details\n*  ===============\n*\n*  We first consider Rectangular Full Packed (RFP) Format when N is\n*  even. We give an example where N = 6.\n*\n*      AP is Upper             AP is Lower\n*\n*   00 01 02 03 04 05       00\n*      11 12 13 14 15       10 11\n*         22 23 24 25       20 21 22\n*            33 34 35       30 31 32 33\n*               44 45       40 41 42 43 44\n*                  55       50 51 52 53 54 55\n*\n*\n*  Let TRANSR = 'N'. RFP holds AP as follows:\n*  For UPLO = 'U' the upper trapezoid A(0:5,0:2) consists of the last\n*  three columns of AP upper. The lower triangle A(4:6,0:2) consists of\n*  the transpose of the first three columns of AP upper.\n*  For UPLO = 'L' the lower trapezoid A(1:6,0:2) consists of the first\n*  three columns of AP lower. The upper triangle A(0:2,0:2) consists of\n*  the transpose of the last three columns of AP lower.\n*  This covers the case N even and TRANSR = 'N'.\n*\n*         RFP A                   RFP A\n*\n*        03 04 05                33 43 53\n*        13 14 15                00 44 54\n*        23 24 25                10 11 55\n*        33 34 35                20 21 22\n*        00 44 45                30 31 32\n*        01 11 55                40 41 42\n*        02 12 22                50 51 52\n*\n*  Now let TRANSR = 'T'. RFP A in both UPLO cases is just the\n*  transpose of RFP A above. One therefore gets:\n*\n*\n*           RFP A                   RFP A\n*\n*     03 13 23 33 00 01 02    33 00 10 20 30 40 50\n*     04 14 24 34 44 11 12    43 44 11 21 31 41 51\n*     05 15 25 35 45 55 22    53 54 55 22 32 42 52\n*\n*\n*  We then consider Rectangular Full Packed (RFP) Format when N is\n*  odd. We give an example where N = 5.\n*\n*     AP is Upper                 AP is Lower\n*\n*   00 01 02 03 04              00\n*      11 12 13 14              10 11\n*         22 23 24              20 21 22\n*            33 34              30 31 32 33\n*               44              40 41 42 43 44\n*\n*\n*  Let TRANSR = 'N'. RFP holds AP as follows:\n*  For UPLO = 'U' the upper trapezoid A(0:4,0:2) consists of the last\n*  three columns of AP upper. The lower triangle A(3:4,0:1) consists of\n*  the transpose of the first two columns of AP upper.\n*  For UPLO = 'L' the lower trapezoid A(0:4,0:2) consists of the first\n*  three columns of AP lower. The upper triangle A(0:1,1:2) consists of\n*  the transpose of the last two columns of AP lower.\n*  This covers the case N odd and TRANSR = 'N'.\n*\n*         RFP A                   RFP A\n*\n*        02 03 04                00 33 43\n*        12 13 14                10 11 44\n*        22 23 24                20 21 22\n*        00 33 34                30 31 32\n*        01 11 44                40 41 42\n*\n*  Now let TRANSR = 'T'. RFP A in both UPLO cases is just the\n*  transpose of RFP A above. One therefore gets:\n*\n*           RFP A                   RFP A\n*\n*     02 12 22 00 01             00 10 20 30 40 50\n*     03 13 23 33 11             33 11 21 31 41 51\n*     04 14 24 34 44             43 44 22 32 42 52\n*\n*  Reference\n*  =========\n*\n*  =====================================================================\n*\n*     ..\n\n");
    return Qnil;
  }
  if (argc != 5)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 5)", argc);
  rb_norm = argv[0];
  rb_transr = argv[1];
  rb_uplo = argv[2];
  rb_n = argv[3];
  rb_a = argv[4];

  transr = StringValueCStr(rb_transr)[0];
  n = NUM2INT(rb_n);
  uplo = StringValueCStr(rb_uplo)[0];
  norm = StringValueCStr(rb_norm)[0];
  if (!NA_IsNArray(rb_a))
    rb_raise(rb_eArgError, "a (5th argument) must be NArray");
  if (NA_RANK(rb_a) != 1)
    rb_raise(rb_eArgError, "rank of a (5th argument) must be %d", 1);
  if (NA_SHAPE0(rb_a) != (n*(n+1)/2))
    rb_raise(rb_eRuntimeError, "shape 0 of a must be %d", n*(n+1)/2);
  if (NA_TYPE(rb_a) != NA_SFLOAT)
    rb_a = na_change_type(rb_a, NA_SFLOAT);
  a = NA_PTR_TYPE(rb_a, real*);
  work = ALLOC_N(real, (MAX(1,(lsame_(&norm,"I")||lsame_(&norm,"1")||lsame_(&norm,"o")) ? n : 0)));

  __out__ = slansf_(&norm, &transr, &uplo, &n, a, work);

  free(work);
  rb___out__ = rb_float_new((double)__out__);
  return rb___out__;
}

void
init_lapack_slansf(VALUE mLapack){
  rb_define_module_function(mLapack, "slansf", rb_slansf, -1);
}
