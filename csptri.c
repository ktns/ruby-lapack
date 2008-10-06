#include "rb_lapack.h"

static VALUE
rb_csptri(int argc, VALUE *argv, VALUE self){
  VALUE rb_uplo;
  char uplo; 
  VALUE rb_ap;
  complex *ap; 
  VALUE rb_ipiv;
  integer *ipiv; 
  VALUE rb_info;
  integer info; 
  VALUE rb_ap_out__;
  complex *ap_out__;
  complex *work;

  integer n;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  info, ap = NumRu::Lapack.csptri( uplo, ap, ipiv)\n    or\n  NumRu::Lapack.csptri  # print help\n\n\nFORTRAN MANUAL\n      SUBROUTINE CSPTRI( UPLO, N, AP, IPIV, WORK, INFO )\n\n*  Purpose\n*  =======\n*\n*  CSPTRI computes the inverse of a complex symmetric indefinite matrix\n*  A in packed storage using the factorization A = U*D*U**T or\n*  A = L*D*L**T computed by CSPTRF.\n*\n\n*  Arguments\n*  =========\n*\n*  UPLO    (input) CHARACTER*1\n*          Specifies whether the details of the factorization are stored\n*          as an upper or lower triangular matrix.\n*          = 'U':  Upper triangular, form is A = U*D*U**T;\n*          = 'L':  Lower triangular, form is A = L*D*L**T.\n*\n*  N       (input) INTEGER\n*          The order of the matrix A.  N >= 0.\n*\n*  AP      (input/output) COMPLEX array, dimension (N*(N+1)/2)\n*          On entry, the block diagonal matrix D and the multipliers\n*          used to obtain the factor U or L as computed by CSPTRF,\n*          stored as a packed triangular matrix.\n*\n*          On exit, if INFO = 0, the (symmetric) inverse of the original\n*          matrix, stored as a packed triangular matrix. The j-th column\n*          of inv(A) is stored in the array AP as follows:\n*          if UPLO = 'U', AP(i + (j-1)*j/2) = inv(A)(i,j) for 1<=i<=j;\n*          if UPLO = 'L',\n*             AP(i + (j-1)*(2n-j)/2) = inv(A)(i,j) for j<=i<=n.\n*\n*  IPIV    (input) INTEGER array, dimension (N)\n*          Details of the interchanges and the block structure of D\n*          as determined by CSPTRF.\n*\n*  WORK    (workspace) COMPLEX array, dimension (N)\n*\n*  INFO    (output) INTEGER\n*          = 0: successful exit\n*          < 0: if INFO = -i, the i-th argument had an illegal value\n*          > 0: if INFO = i, D(i,i) = 0; the matrix is singular and its\n*               inverse could not be computed.\n*\n\n*  =====================================================================\n*\n\n");
    return Qnil;
  }
  if (argc != 3)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 3)", argc);
  rb_uplo = argv[0];
  rb_ap = argv[1];
  rb_ipiv = argv[2];

  uplo = StringValueCStr(rb_uplo)[0];
  if (!NA_IsNArray(rb_ipiv))
    rb_raise(rb_eArgError, "ipiv (2th argument) must be NArray");
  if (NA_RANK(rb_ipiv) != 1)
    rb_raise(rb_eArgError, "rank of ipiv (2th argument) must be %d", 1);
  n = NA_SHAPE0(rb_ipiv);
  if (NA_TYPE(rb_ipiv) != NA_LINT)
    rb_ipiv = na_change_type(rb_ipiv, NA_LINT);
  ipiv = NA_PTR_TYPE(rb_ipiv, integer*);
  if (!NA_IsNArray(rb_ap))
    rb_raise(rb_eArgError, "ap (3th argument) must be NArray");
  if (NA_RANK(rb_ap) != 1)
    rb_raise(rb_eArgError, "rank of ap (3th argument) must be %d", 1);
  if (NA_SHAPE0(rb_ap) != (n*(n+1)/2))
    rb_raise(rb_eRuntimeError, "shape 0 of ap must be %d", n*(n+1)/2);
  if (NA_TYPE(rb_ap) != NA_SCOMPLEX)
    rb_ap = na_change_type(rb_ap, NA_SCOMPLEX);
  ap = NA_PTR_TYPE(rb_ap, complex*);
  {
    int shape[1];
    shape[0] = n*(n+1)/2;
    rb_ap_out__ = na_make_object(NA_SCOMPLEX, 1, shape, cNArray);
  }
  ap_out__ = NA_PTR_TYPE(rb_ap_out__, complex*);
  MEMCPY(ap_out__, ap, complex, NA_TOTAL(rb_ap));
  rb_ap = rb_ap_out__;
  ap = ap_out__;
  work = ALLOC_N(complex, (n));

  csptri_(&uplo, &n, ap, ipiv, work, &info);

  free(work);
  rb_info = INT2NUM(info);
  return rb_ary_new3(2, rb_info, rb_ap);
}

void
init_lapack_csptri(VALUE mLapack){
  rb_define_module_function(mLapack, "csptri", rb_csptri, -1);
}
