#include "rb_lapack.h"

static VALUE
rb_ztfttp(int argc, VALUE *argv, VALUE self){
  VALUE rb_transr;
  char transr; 
  VALUE rb_uplo;
  char uplo; 
  VALUE rb_n;
  integer n; 
  VALUE rb_arf;
  doublecomplex *arf; 
  VALUE rb_ap;
  doublecomplex *ap; 
  VALUE rb_info;
  integer info; 


  if (argc == 0) {
    printf("%s\n", "USAGE:\n  ap, info = NumRu::Lapack.ztfttp( transr, uplo, n, arf)\n    or\n  NumRu::Lapack.ztfttp  # print help\n\n\nFORTRAN MANUAL\n      SUBROUTINE ZTFTTP( TRANSR, UPLO, N, ARF, AP, INFO )\n\n*  Purpose\n*  =======\n*\n*  ZTFTTP copies a triangular matrix A from rectangular full packed\n*  format (TF) to standard packed format (TP).\n*\n\n*  Arguments\n*  =========\n*\n*  TRANSR   (input) CHARACTER\n*          = 'N':  ARF is in Normal format;\n*          = 'C':  ARF is in Conjugate-transpose format;\n*\n*  UPLO    (input) CHARACTER\n*          = 'U':  A is upper triangular;\n*          = 'L':  A is lower triangular.\n*\n*  N       (input) INTEGER\n*          The order of the matrix A. N >= 0.\n*\n*  ARF     (input) COMPLEX*16 array, dimension ( N*(N+1)/2 ),\n*          On entry, the upper or lower triangular matrix A stored in\n*          RFP format. For a further discussion see Notes below.\n*\n*  AP      (output) COMPLEX*16 array, dimension ( N*(N+1)/2 ),\n*          On exit, the upper or lower triangular matrix A, packed\n*          columnwise in a linear array. The j-th column of A is stored\n*          in the array AP as follows:\n*          if UPLO = 'U', AP(i + (j-1)*j/2) = A(i,j) for 1<=i<=j;\n*          if UPLO = 'L', AP(i + (j-1)*(2n-j)/2) = A(i,j) for j<=i<=n.\n*\n*  INFO    (output) INTEGER\n*          = 0:  successful exit\n*          < 0:  if INFO = -i, the i-th argument had an illegal value\n*\n\n*  Further Details\n*  ===============\n*\n*  We first consider Standard Packed Format when N is even.\n*  We give an example where N = 6.\n*\n*      AP is Upper             AP is Lower\n*\n*   00 01 02 03 04 05       00\n*      11 12 13 14 15       10 11\n*         22 23 24 25       20 21 22\n*            33 34 35       30 31 32 33\n*               44 45       40 41 42 43 44\n*                  55       50 51 52 53 54 55\n*\n*\n*  Let TRANSR = 'N'. RFP holds AP as follows:\n*  For UPLO = 'U' the upper trapezoid A(0:5,0:2) consists of the last\n*  three columns of AP upper. The lower triangle A(4:6,0:2) consists of\n*  conjugate-transpose of the first three columns of AP upper.\n*  For UPLO = 'L' the lower trapezoid A(1:6,0:2) consists of the first\n*  three columns of AP lower. The upper triangle A(0:2,0:2) consists of\n*  conjugate-transpose of the last three columns of AP lower.\n*  To denote conjugate we place -- above the element. This covers the\n*  case N even and TRANSR = 'N'.\n*\n*         RFP A                   RFP A\n*\n*                                -- -- --\n*        03 04 05                33 43 53\n*                                   -- --\n*        13 14 15                00 44 54\n*                                      --\n*        23 24 25                10 11 55\n*\n*        33 34 35                20 21 22\n*        --\n*        00 44 45                30 31 32\n*        -- --\n*        01 11 55                40 41 42\n*        -- -- --\n*        02 12 22                50 51 52\n*\n*  Now let TRANSR = 'C'. RFP A in both UPLO cases is just the conjugate-\n*  transpose of RFP A above. One therefore gets:\n*\n*\n*           RFP A                   RFP A\n*\n*     -- -- -- --                -- -- -- -- -- --\n*     03 13 23 33 00 01 02    33 00 10 20 30 40 50\n*     -- -- -- -- --                -- -- -- -- --\n*     04 14 24 34 44 11 12    43 44 11 21 31 41 51\n*     -- -- -- -- -- --                -- -- -- --\n*     05 15 25 35 45 55 22    53 54 55 22 32 42 52\n*\n*\n*  We next consider Standard Packed Format when N is odd.\n*  We give an example where N = 5.\n*\n*     AP is Upper                 AP is Lower\n*\n*   00 01 02 03 04              00\n*      11 12 13 14              10 11\n*         22 23 24              20 21 22\n*            33 34              30 31 32 33\n*               44              40 41 42 43 44\n*\n*\n*  Let TRANSR = 'N'. RFP holds AP as follows:\n*  For UPLO = 'U' the upper trapezoid A(0:4,0:2) consists of the last\n*  three columns of AP upper. The lower triangle A(3:4,0:1) consists of\n*  conjugate-transpose of the first two   columns of AP upper.\n*  For UPLO = 'L' the lower trapezoid A(0:4,0:2) consists of the first\n*  three columns of AP lower. The upper triangle A(0:1,1:2) consists of\n*  conjugate-transpose of the last two   columns of AP lower.\n*  To denote conjugate we place -- above the element. This covers the\n*  case N odd  and TRANSR = 'N'.\n*\n*         RFP A                   RFP A\n*\n*                                   -- --\n*        02 03 04                00 33 43\n*                                      --\n*        12 13 14                10 11 44\n*\n*        22 23 24                20 21 22\n*        --\n*        00 33 34                30 31 32\n*        -- --\n*        01 11 44                40 41 42\n*\n*  Now let TRANSR = 'C'. RFP A in both UPLO cases is just the conjugate-\n*  transpose of RFP A above. One therefore gets:\n*\n*\n*           RFP A                   RFP A\n*\n*     -- -- --                   -- -- -- -- -- --\n*     02 12 22 00 01             00 10 20 30 40 50\n*     -- -- -- --                   -- -- -- -- --\n*     03 13 23 33 11             33 11 21 31 41 51\n*     -- -- -- -- --                   -- -- -- --\n*     04 14 24 34 44             43 44 22 32 42 52\n*\n*  =====================================================================\n*\n\n");
    return Qnil;
  }
  if (argc != 4)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 4)", argc);
  rb_transr = argv[0];
  rb_uplo = argv[1];
  rb_n = argv[2];
  rb_arf = argv[3];

  transr = StringValueCStr(rb_transr)[0];
  uplo = StringValueCStr(rb_uplo)[0];
  n = NUM2INT(rb_n);
  if (!NA_IsNArray(rb_arf))
    rb_raise(rb_eArgError, "arf (4th argument) must be NArray");
  if (NA_RANK(rb_arf) != 1)
    rb_raise(rb_eArgError, "rank of arf (4th argument) must be %d", 1);
  if (NA_SHAPE0(rb_arf) != (( n*(n+1)/2 )))
    rb_raise(rb_eRuntimeError, "shape 0 of arf must be %d", ( n*(n+1)/2 ));
  if (NA_TYPE(rb_arf) != NA_DCOMPLEX)
    rb_arf = na_change_type(rb_arf, NA_DCOMPLEX);
  arf = NA_PTR_TYPE(rb_arf, doublecomplex*);
  {
    int shape[1];
    shape[0] = ( n*(n+1)/2 );
    rb_ap = na_make_object(NA_DCOMPLEX, 1, shape, cNArray);
  }
  ap = NA_PTR_TYPE(rb_ap, doublecomplex*);

  ztfttp_(&transr, &uplo, &n, arf, ap, &info);

  rb_info = INT2NUM(info);
  return rb_ary_new3(2, rb_ap, rb_info);
}

void
init_lapack_ztfttp(VALUE mLapack){
  rb_define_module_function(mLapack, "ztfttp", rb_ztfttp, -1);
}
