#include "rb_lapack.h"

extern doublereal dlange_(char *norm, integer *m, integer *n, doublereal *a, integer *lda, doublereal *work);

static VALUE sHelp, sUsage;

static VALUE
rblapack_dlange(int argc, VALUE *argv, VALUE self){
  VALUE rblapack_norm;
  char norm; 
  VALUE rblapack_m;
  integer m; 
  VALUE rblapack_a;
  doublereal *a; 
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
      printf("%s\n", "USAGE:\n  __out__ = NumRu::Lapack.dlange( norm, m, a, [:usage => usage, :help => help])\n\n\nFORTRAN MANUAL\n      DOUBLE PRECISION FUNCTION DLANGE( NORM, M, N, A, LDA, WORK )\n\n*  Purpose\n*  =======\n*\n*  DLANGE  returns the value of the one norm,  or the Frobenius norm, or\n*  the  infinity norm,  or the  element of  largest absolute value  of a\n*  real matrix A.\n*\n*  Description\n*  ===========\n*\n*  DLANGE returns the value\n*\n*     DLANGE = ( max(abs(A(i,j))), NORM = 'M' or 'm'\n*              (\n*              ( norm1(A),         NORM = '1', 'O' or 'o'\n*              (\n*              ( normI(A),         NORM = 'I' or 'i'\n*              (\n*              ( normF(A),         NORM = 'F', 'f', 'E' or 'e'\n*\n*  where  norm1  denotes the  one norm of a matrix (maximum column sum),\n*  normI  denotes the  infinity norm  of a matrix  (maximum row sum) and\n*  normF  denotes the  Frobenius norm of a matrix (square root of sum of\n*  squares).  Note that  max(abs(A(i,j)))  is not a consistent matrix norm.\n*\n\n*  Arguments\n*  =========\n*\n*  NORM    (input) CHARACTER*1\n*          Specifies the value to be returned in DLANGE as described\n*          above.\n*\n*  M       (input) INTEGER\n*          The number of rows of the matrix A.  M >= 0.  When M = 0,\n*          DLANGE is set to zero.\n*\n*  N       (input) INTEGER\n*          The number of columns of the matrix A.  N >= 0.  When N = 0,\n*          DLANGE is set to zero.\n*\n*  A       (input) DOUBLE PRECISION array, dimension (LDA,N)\n*          The m by n matrix A.\n*\n*  LDA     (input) INTEGER\n*          The leading dimension of the array A.  LDA >= max(M,1).\n*\n*  WORK    (workspace) DOUBLE PRECISION array, dimension (MAX(1,LWORK)),\n*          where LWORK >= M when NORM = 'I'; otherwise, WORK is not\n*          referenced.\n*\n\n* =====================================================================\n*\n\n");
      rb_exit(0);
    }
    if (rb_hash_aref(rb_options, sUsage) == Qtrue) {
      printf("%s\n", "USAGE:\n  __out__ = NumRu::Lapack.dlange( norm, m, a, [:usage => usage, :help => help])\n");
      rb_exit(0);
    } 
  } else
    rb_options = Qnil;
  if (argc != 3)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 3)", argc);
  rblapack_norm = argv[0];
  rblapack_m = argv[1];
  rblapack_a = argv[2];
  if (rb_options != Qnil) {
  }

  if (!NA_IsNArray(rblapack_a))
    rb_raise(rb_eArgError, "a (3th argument) must be NArray");
  if (NA_RANK(rblapack_a) != 2)
    rb_raise(rb_eArgError, "rank of a (3th argument) must be %d", 2);
  n = NA_SHAPE1(rblapack_a);
  lda = NA_SHAPE0(rblapack_a);
  if (NA_TYPE(rblapack_a) != NA_DFLOAT)
    rblapack_a = na_change_type(rblapack_a, NA_DFLOAT);
  a = NA_PTR_TYPE(rblapack_a, doublereal*);
  m = NUM2INT(rblapack_m);
  norm = StringValueCStr(rblapack_norm)[0];
  lwork = lsame_(&norm,"I") ? m : 0;
  work = ALLOC_N(doublereal, (MAX(1,lwork)));

  __out__ = dlange_(&norm, &m, &n, a, &lda, work);

  free(work);
  rblapack___out__ = rb_float_new((double)__out__);
  return rblapack___out__;
}

void
init_lapack_dlange(VALUE mLapack){
  rb_define_module_function(mLapack, "dlange", rblapack_dlange, -1);
  sHelp = ID2SYM(rb_intern("help"));
  sUsage = ID2SYM(rb_intern("usage"));
}
