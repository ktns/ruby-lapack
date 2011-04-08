#include "rb_lapack.h"

extern real clangb_(char *norm, integer *n, integer *kl, integer *ku, complex *ab, integer *ldab, real *work);

static VALUE sHelp, sUsage;

static VALUE
rblapack_clangb(int argc, VALUE *argv, VALUE self){
  VALUE rblapack_norm;
  char norm; 
  VALUE rblapack_kl;
  integer kl; 
  VALUE rblapack_ku;
  integer ku; 
  VALUE rblapack_ab;
  complex *ab; 
  VALUE rblapack___out__;
  real __out__; 
  real *work;

  integer ldab;
  integer n;

  VALUE rb_options;
  if (argc > 0 && TYPE(argv[argc-1]) == T_HASH) {
    argc--;
    rb_options = argv[argc];
    if (rb_hash_aref(rb_options, sHelp) == Qtrue) {
      printf("%s\n", "USAGE:\n  __out__ = NumRu::Lapack.clangb( norm, kl, ku, ab, [:usage => usage, :help => help])\n\n\nFORTRAN MANUAL\n      REAL             FUNCTION CLANGB( NORM, N, KL, KU, AB, LDAB, WORK )\n\n*  Purpose\n*  =======\n*\n*  CLANGB  returns the value of the one norm,  or the Frobenius norm, or\n*  the  infinity norm,  or the element of  largest absolute value  of an\n*  n by n band matrix  A,  with kl sub-diagonals and ku super-diagonals.\n*\n*  Description\n*  ===========\n*\n*  CLANGB returns the value\n*\n*     CLANGB = ( max(abs(A(i,j))), NORM = 'M' or 'm'\n*              (\n*              ( norm1(A),         NORM = '1', 'O' or 'o'\n*              (\n*              ( normI(A),         NORM = 'I' or 'i'\n*              (\n*              ( normF(A),         NORM = 'F', 'f', 'E' or 'e'\n*\n*  where  norm1  denotes the  one norm of a matrix (maximum column sum),\n*  normI  denotes the  infinity norm  of a matrix  (maximum row sum) and\n*  normF  denotes the  Frobenius norm of a matrix (square root of sum of\n*  squares).  Note that  max(abs(A(i,j)))  is not a consistent matrix norm.\n*\n\n*  Arguments\n*  =========\n*\n*  NORM    (input) CHARACTER*1\n*          Specifies the value to be returned in CLANGB as described\n*          above.\n*\n*  N       (input) INTEGER\n*          The order of the matrix A.  N >= 0.  When N = 0, CLANGB is\n*          set to zero.\n*\n*  KL      (input) INTEGER\n*          The number of sub-diagonals of the matrix A.  KL >= 0.\n*\n*  KU      (input) INTEGER\n*          The number of super-diagonals of the matrix A.  KU >= 0.\n*\n*  AB      (input) COMPLEX array, dimension (LDAB,N)\n*          The band matrix A, stored in rows 1 to KL+KU+1.  The j-th\n*          column of A is stored in the j-th column of the array AB as\n*          follows:\n*          AB(ku+1+i-j,j) = A(i,j) for max(1,j-ku)<=i<=min(n,j+kl).\n*\n*  LDAB    (input) INTEGER\n*          The leading dimension of the array AB.  LDAB >= KL+KU+1.\n*\n*  WORK    (workspace) REAL array, dimension (MAX(1,LWORK)),\n*          where LWORK >= N when NORM = 'I'; otherwise, WORK is not\n*          referenced.\n*\n\n* =====================================================================\n*\n\n");
      rb_exit(0);
    }
    if (rb_hash_aref(rb_options, sUsage) == Qtrue) {
      printf("%s\n", "USAGE:\n  __out__ = NumRu::Lapack.clangb( norm, kl, ku, ab, [:usage => usage, :help => help])\n");
      rb_exit(0);
    } 
  } else
    rb_options = Qnil;
  if (argc != 4)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 4)", argc);
  rblapack_norm = argv[0];
  rblapack_kl = argv[1];
  rblapack_ku = argv[2];
  rblapack_ab = argv[3];
  if (rb_options != Qnil) {
  }

  if (!NA_IsNArray(rblapack_ab))
    rb_raise(rb_eArgError, "ab (4th argument) must be NArray");
  if (NA_RANK(rblapack_ab) != 2)
    rb_raise(rb_eArgError, "rank of ab (4th argument) must be %d", 2);
  n = NA_SHAPE1(rblapack_ab);
  ldab = NA_SHAPE0(rblapack_ab);
  if (NA_TYPE(rblapack_ab) != NA_SCOMPLEX)
    rblapack_ab = na_change_type(rblapack_ab, NA_SCOMPLEX);
  ab = NA_PTR_TYPE(rblapack_ab, complex*);
  kl = NUM2INT(rblapack_kl);
  norm = StringValueCStr(rblapack_norm)[0];
  ku = NUM2INT(rblapack_ku);
  work = ALLOC_N(real, (MAX(1,lsame_(&norm,"I") ? n : 0)));

  __out__ = clangb_(&norm, &n, &kl, &ku, ab, &ldab, work);

  free(work);
  rblapack___out__ = rb_float_new((double)__out__);
  return rblapack___out__;
}

void
init_lapack_clangb(VALUE mLapack){
  rb_define_module_function(mLapack, "clangb", rblapack_clangb, -1);
  sHelp = ID2SYM(rb_intern("help"));
  sUsage = ID2SYM(rb_intern("usage"));
}
