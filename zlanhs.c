#include "rb_lapack.h"

extern doublereal zlanhs_(char *norm, integer *n, doublecomplex *a, integer *lda, doublereal *work);

static VALUE sHelp, sUsage;

static VALUE
rblapack_zlanhs(int argc, VALUE *argv, VALUE self){
  VALUE rblapack_norm;
  char norm; 
  VALUE rblapack_a;
  doublecomplex *a; 
  VALUE rblapack___out__;
  doublereal __out__; 
  doublereal *work;

  integer lda;
  integer n;
  integer lwork;

  VALUE rb_options;
  if (argc > 0 && TYPE(argv[argc-1]) == T_HASH) {
    argc--;
    rb_options = argv[argc];
    if (rb_hash_aref(rb_options, sHelp) == Qtrue) {
      printf("%s\n", "USAGE:\n  __out__ = NumRu::Lapack.zlanhs( norm, a, [:usage => usage, :help => help])\n\n\nFORTRAN MANUAL\n      DOUBLE PRECISION FUNCTION ZLANHS( NORM, N, A, LDA, WORK )\n\n*  Purpose\n*  =======\n*\n*  ZLANHS  returns the value of the one norm,  or the Frobenius norm, or\n*  the  infinity norm,  or the  element of  largest absolute value  of a\n*  Hessenberg matrix A.\n*\n*  Description\n*  ===========\n*\n*  ZLANHS returns the value\n*\n*     ZLANHS = ( max(abs(A(i,j))), NORM = 'M' or 'm'\n*              (\n*              ( norm1(A),         NORM = '1', 'O' or 'o'\n*              (\n*              ( normI(A),         NORM = 'I' or 'i'\n*              (\n*              ( normF(A),         NORM = 'F', 'f', 'E' or 'e'\n*\n*  where  norm1  denotes the  one norm of a matrix (maximum column sum),\n*  normI  denotes the  infinity norm  of a matrix  (maximum row sum) and\n*  normF  denotes the  Frobenius norm of a matrix (square root of sum of\n*  squares).  Note that  max(abs(A(i,j)))  is not a consistent matrix norm.\n*\n\n*  Arguments\n*  =========\n*\n*  NORM    (input) CHARACTER*1\n*          Specifies the value to be returned in ZLANHS as described\n*          above.\n*\n*  N       (input) INTEGER\n*          The order of the matrix A.  N >= 0.  When N = 0, ZLANHS is\n*          set to zero.\n*\n*  A       (input) COMPLEX*16 array, dimension (LDA,N)\n*          The n by n upper Hessenberg matrix A; the part of A below the\n*          first sub-diagonal is not referenced.\n*\n*  LDA     (input) INTEGER\n*          The leading dimension of the array A.  LDA >= max(N,1).\n*\n*  WORK    (workspace) DOUBLE PRECISION array, dimension (MAX(1,LWORK)),\n*          where LWORK >= N when NORM = 'I'; otherwise, WORK is not\n*          referenced.\n*\n\n* =====================================================================\n*\n\n");
      rb_exit(0);
    }
    if (rb_hash_aref(rb_options, sUsage) == Qtrue) {
      printf("%s\n", "USAGE:\n  __out__ = NumRu::Lapack.zlanhs( norm, a, [:usage => usage, :help => help])\n");
      rb_exit(0);
    } 
  } else
    rb_options = Qnil;
  if (argc != 2)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 2)", argc);
  rblapack_norm = argv[0];
  rblapack_a = argv[1];
  if (rb_options != Qnil) {
  }

  if (!NA_IsNArray(rblapack_a))
    rb_raise(rb_eArgError, "a (2th argument) must be NArray");
  if (NA_RANK(rblapack_a) != 2)
    rb_raise(rb_eArgError, "rank of a (2th argument) must be %d", 2);
  n = NA_SHAPE1(rblapack_a);
  lda = NA_SHAPE0(rblapack_a);
  if (NA_TYPE(rblapack_a) != NA_DCOMPLEX)
    rblapack_a = na_change_type(rblapack_a, NA_DCOMPLEX);
  a = NA_PTR_TYPE(rblapack_a, doublecomplex*);
  norm = StringValueCStr(rblapack_norm)[0];
  lwork = lsame_(&norm,"I") ? n : 0;
  work = ALLOC_N(doublereal, (MAX(1,lwork)));

  __out__ = zlanhs_(&norm, &n, a, &lda, work);

  free(work);
  rblapack___out__ = rb_float_new((double)__out__);
  return rblapack___out__;
}

void
init_lapack_zlanhs(VALUE mLapack){
  rb_define_module_function(mLapack, "zlanhs", rblapack_zlanhs, -1);
  sHelp = ID2SYM(rb_intern("help"));
  sUsage = ID2SYM(rb_intern("usage"));
}
