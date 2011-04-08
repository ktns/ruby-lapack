#include "rb_lapack.h"

extern VOID dlasrt_(char *id, integer *n, doublereal *d, integer *info);

static VALUE sHelp, sUsage;

static VALUE
rblapack_dlasrt(int argc, VALUE *argv, VALUE self){
  VALUE rblapack_id;
  char id; 
  VALUE rblapack_d;
  doublereal *d; 
  VALUE rblapack_info;
  integer info; 
  VALUE rblapack_d_out__;
  doublereal *d_out__;

  integer n;

  VALUE rb_options;
  if (argc > 0 && TYPE(argv[argc-1]) == T_HASH) {
    argc--;
    rb_options = argv[argc];
    if (rb_hash_aref(rb_options, sHelp) == Qtrue) {
      printf("%s\n", "USAGE:\n  info, d = NumRu::Lapack.dlasrt( id, d, [:usage => usage, :help => help])\n\n\nFORTRAN MANUAL\n      SUBROUTINE DLASRT( ID, N, D, INFO )\n\n*  Purpose\n*  =======\n*\n*  Sort the numbers in D in increasing order (if ID = 'I') or\n*  in decreasing order (if ID = 'D' ).\n*\n*  Use Quick Sort, reverting to Insertion sort on arrays of\n*  size <= 20. Dimension of STACK limits N to about 2**32.\n*\n\n*  Arguments\n*  =========\n*\n*  ID      (input) CHARACTER*1\n*          = 'I': sort D in increasing order;\n*          = 'D': sort D in decreasing order.\n*\n*  N       (input) INTEGER\n*          The length of the array D.\n*\n*  D       (input/output) DOUBLE PRECISION array, dimension (N)\n*          On entry, the array to be sorted.\n*          On exit, D has been sorted into increasing order\n*          (D(1) <= ... <= D(N) ) or into decreasing order\n*          (D(1) >= ... >= D(N) ), depending on ID.\n*\n*  INFO    (output) INTEGER\n*          = 0:  successful exit\n*          < 0:  if INFO = -i, the i-th argument had an illegal value\n*\n\n*  =====================================================================\n*\n\n");
      rb_exit(0);
    }
    if (rb_hash_aref(rb_options, sUsage) == Qtrue) {
      printf("%s\n", "USAGE:\n  info, d = NumRu::Lapack.dlasrt( id, d, [:usage => usage, :help => help])\n");
      rb_exit(0);
    } 
  } else
    rb_options = Qnil;
  if (argc != 2)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 2)", argc);
  rblapack_id = argv[0];
  rblapack_d = argv[1];
  if (rb_options != Qnil) {
  }

  if (!NA_IsNArray(rblapack_d))
    rb_raise(rb_eArgError, "d (2th argument) must be NArray");
  if (NA_RANK(rblapack_d) != 1)
    rb_raise(rb_eArgError, "rank of d (2th argument) must be %d", 1);
  n = NA_SHAPE0(rblapack_d);
  if (NA_TYPE(rblapack_d) != NA_DFLOAT)
    rblapack_d = na_change_type(rblapack_d, NA_DFLOAT);
  d = NA_PTR_TYPE(rblapack_d, doublereal*);
  id = StringValueCStr(rblapack_id)[0];
  {
    int shape[1];
    shape[0] = n;
    rblapack_d_out__ = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  d_out__ = NA_PTR_TYPE(rblapack_d_out__, doublereal*);
  MEMCPY(d_out__, d, doublereal, NA_TOTAL(rblapack_d));
  rblapack_d = rblapack_d_out__;
  d = d_out__;

  dlasrt_(&id, &n, d, &info);

  rblapack_info = INT2NUM(info);
  return rb_ary_new3(2, rblapack_info, rblapack_d);
}

void
init_lapack_dlasrt(VALUE mLapack){
  rb_define_module_function(mLapack, "dlasrt", rblapack_dlasrt, -1);
  sHelp = ID2SYM(rb_intern("help"));
  sUsage = ID2SYM(rb_intern("usage"));
}
