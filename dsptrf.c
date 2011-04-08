#include "rb_lapack.h"

extern VOID dsptrf_(char *uplo, integer *n, doublereal *ap, integer *ipiv, integer *info);

static VALUE sHelp, sUsage;

static VALUE
rblapack_dsptrf(int argc, VALUE *argv, VALUE self){
  VALUE rblapack_uplo;
  char uplo; 
  VALUE rblapack_ap;
  doublereal *ap; 
  VALUE rblapack_ipiv;
  integer *ipiv; 
  VALUE rblapack_info;
  integer info; 
  VALUE rblapack_ap_out__;
  doublereal *ap_out__;

  integer ldap;
  integer n;

  VALUE rb_options;
  if (argc > 0 && TYPE(argv[argc-1]) == T_HASH) {
    argc--;
    rb_options = argv[argc];
    if (rb_hash_aref(rb_options, sHelp) == Qtrue) {
      printf("%s\n", "USAGE:\n  ipiv, info, ap = NumRu::Lapack.dsptrf( uplo, ap, [:usage => usage, :help => help])\n\n\nFORTRAN MANUAL\n      SUBROUTINE DSPTRF( UPLO, N, AP, IPIV, INFO )\n\n*  Purpose\n*  =======\n*\n*  DSPTRF computes the factorization of a real symmetric matrix A stored\n*  in packed format using the Bunch-Kaufman diagonal pivoting method:\n*\n*     A = U*D*U**T  or  A = L*D*L**T\n*\n*  where U (or L) is a product of permutation and unit upper (lower)\n*  triangular matrices, and D is symmetric and block diagonal with\n*  1-by-1 and 2-by-2 diagonal blocks.\n*\n\n*  Arguments\n*  =========\n*\n*  UPLO    (input) CHARACTER*1\n*          = 'U':  Upper triangle of A is stored;\n*          = 'L':  Lower triangle of A is stored.\n*\n*  N       (input) INTEGER\n*          The order of the matrix A.  N >= 0.\n*\n*  AP      (input/output) DOUBLE PRECISION array, dimension (N*(N+1)/2)\n*          On entry, the upper or lower triangle of the symmetric matrix\n*          A, packed columnwise in a linear array.  The j-th column of A\n*          is stored in the array AP as follows:\n*          if UPLO = 'U', AP(i + (j-1)*j/2) = A(i,j) for 1<=i<=j;\n*          if UPLO = 'L', AP(i + (j-1)*(2n-j)/2) = A(i,j) for j<=i<=n.\n*\n*          On exit, the block diagonal matrix D and the multipliers used\n*          to obtain the factor U or L, stored as a packed triangular\n*          matrix overwriting A (see below for further details).\n*\n*  IPIV    (output) INTEGER array, dimension (N)\n*          Details of the interchanges and the block structure of D.\n*          If IPIV(k) > 0, then rows and columns k and IPIV(k) were\n*          interchanged and D(k,k) is a 1-by-1 diagonal block.\n*          If UPLO = 'U' and IPIV(k) = IPIV(k-1) < 0, then rows and\n*          columns k-1 and -IPIV(k) were interchanged and D(k-1:k,k-1:k)\n*          is a 2-by-2 diagonal block.  If UPLO = 'L' and IPIV(k) =\n*          IPIV(k+1) < 0, then rows and columns k+1 and -IPIV(k) were\n*          interchanged and D(k:k+1,k:k+1) is a 2-by-2 diagonal block.\n*\n*  INFO    (output) INTEGER\n*          = 0: successful exit\n*          < 0: if INFO = -i, the i-th argument had an illegal value\n*          > 0: if INFO = i, D(i,i) is exactly zero.  The factorization\n*               has been completed, but the block diagonal matrix D is\n*               exactly singular, and division by zero will occur if it\n*               is used to solve a system of equations.\n*\n\n*  Further Details\n*  ===============\n*\n*  5-96 - Based on modifications by J. Lewis, Boeing Computer Services\n*         Company\n*\n*  If UPLO = 'U', then A = U*D*U', where\n*     U = P(n)*U(n)* ... *P(k)U(k)* ...,\n*  i.e., U is a product of terms P(k)*U(k), where k decreases from n to\n*  1 in steps of 1 or 2, and D is a block diagonal matrix with 1-by-1\n*  and 2-by-2 diagonal blocks D(k).  P(k) is a permutation matrix as\n*  defined by IPIV(k), and U(k) is a unit upper triangular matrix, such\n*  that if the diagonal block D(k) is of order s (s = 1 or 2), then\n*\n*             (   I    v    0   )   k-s\n*     U(k) =  (   0    I    0   )   s\n*             (   0    0    I   )   n-k\n*                k-s   s   n-k\n*\n*  If s = 1, D(k) overwrites A(k,k), and v overwrites A(1:k-1,k).\n*  If s = 2, the upper triangle of D(k) overwrites A(k-1,k-1), A(k-1,k),\n*  and A(k,k), and v overwrites A(1:k-2,k-1:k).\n*\n*  If UPLO = 'L', then A = L*D*L', where\n*     L = P(1)*L(1)* ... *P(k)*L(k)* ...,\n*  i.e., L is a product of terms P(k)*L(k), where k increases from 1 to\n*  n in steps of 1 or 2, and D is a block diagonal matrix with 1-by-1\n*  and 2-by-2 diagonal blocks D(k).  P(k) is a permutation matrix as\n*  defined by IPIV(k), and L(k) is a unit lower triangular matrix, such\n*  that if the diagonal block D(k) is of order s (s = 1 or 2), then\n*\n*             (   I    0     0   )  k-1\n*     L(k) =  (   0    I     0   )  s\n*             (   0    v     I   )  n-k-s+1\n*                k-1   s  n-k-s+1\n*\n*  If s = 1, D(k) overwrites A(k,k), and v overwrites A(k+1:n,k).\n*  If s = 2, the lower triangle of D(k) overwrites A(k,k), A(k+1,k),\n*  and A(k+1,k+1), and v overwrites A(k+2:n,k:k+1).\n*\n*  =====================================================================\n*\n\n");
      rb_exit(0);
    }
    if (rb_hash_aref(rb_options, sUsage) == Qtrue) {
      printf("%s\n", "USAGE:\n  ipiv, info, ap = NumRu::Lapack.dsptrf( uplo, ap, [:usage => usage, :help => help])\n");
      rb_exit(0);
    } 
  } else
    rb_options = Qnil;
  if (argc != 2)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 2)", argc);
  rblapack_uplo = argv[0];
  rblapack_ap = argv[1];
  if (rb_options != Qnil) {
  }

  uplo = StringValueCStr(rblapack_uplo)[0];
  if (!NA_IsNArray(rblapack_ap))
    rb_raise(rb_eArgError, "ap (2th argument) must be NArray");
  if (NA_RANK(rblapack_ap) != 1)
    rb_raise(rb_eArgError, "rank of ap (2th argument) must be %d", 1);
  ldap = NA_SHAPE0(rblapack_ap);
  if (NA_TYPE(rblapack_ap) != NA_DFLOAT)
    rblapack_ap = na_change_type(rblapack_ap, NA_DFLOAT);
  ap = NA_PTR_TYPE(rblapack_ap, doublereal*);
  n = ((int)sqrtf(ldap*8+1.0f)-1)/2;
  {
    int shape[1];
    shape[0] = n;
    rblapack_ipiv = na_make_object(NA_LINT, 1, shape, cNArray);
  }
  ipiv = NA_PTR_TYPE(rblapack_ipiv, integer*);
  {
    int shape[1];
    shape[0] = ldap;
    rblapack_ap_out__ = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  ap_out__ = NA_PTR_TYPE(rblapack_ap_out__, doublereal*);
  MEMCPY(ap_out__, ap, doublereal, NA_TOTAL(rblapack_ap));
  rblapack_ap = rblapack_ap_out__;
  ap = ap_out__;

  dsptrf_(&uplo, &n, ap, ipiv, &info);

  rblapack_info = INT2NUM(info);
  return rb_ary_new3(3, rblapack_ipiv, rblapack_info, rblapack_ap);
}

void
init_lapack_dsptrf(VALUE mLapack){
  rb_define_module_function(mLapack, "dsptrf", rblapack_dsptrf, -1);
  sHelp = ID2SYM(rb_intern("help"));
  sUsage = ID2SYM(rb_intern("usage"));
}
