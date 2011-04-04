#include "rb_lapack.h"

extern VOID dppsvx_(char *fact, char *uplo, integer *n, integer *nrhs, doublereal *ap, doublereal *afp, char *equed, doublereal *s, doublereal *b, integer *ldb, doublereal *x, integer *ldx, doublereal *rcond, doublereal *ferr, doublereal *berr, doublereal *work, integer *iwork, integer *info);

static VALUE
rb_dppsvx(int argc, VALUE *argv, VALUE self){
  VALUE rb_fact;
  char fact; 
  VALUE rb_uplo;
  char uplo; 
  VALUE rb_ap;
  doublereal *ap; 
  VALUE rb_afp;
  doublereal *afp; 
  VALUE rb_equed;
  char equed; 
  VALUE rb_s;
  doublereal *s; 
  VALUE rb_b;
  doublereal *b; 
  VALUE rb_x;
  doublereal *x; 
  VALUE rb_rcond;
  doublereal rcond; 
  VALUE rb_ferr;
  doublereal *ferr; 
  VALUE rb_berr;
  doublereal *berr; 
  VALUE rb_info;
  integer info; 
  VALUE rb_ap_out__;
  doublereal *ap_out__;
  VALUE rb_afp_out__;
  doublereal *afp_out__;
  VALUE rb_s_out__;
  doublereal *s_out__;
  VALUE rb_b_out__;
  doublereal *b_out__;
  doublereal *work;
  integer *iwork;

  integer n;
  integer ldb;
  integer nrhs;
  integer ldx;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  x, rcond, ferr, berr, info, ap, afp, equed, s, b = NumRu::Lapack.dppsvx( fact, uplo, ap, afp, equed, s, b)\n    or\n  NumRu::Lapack.dppsvx  # print help\n\n\nFORTRAN MANUAL\n\n");
    return Qnil;
  }
  if (argc != 7)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 7)", argc);
  rb_fact = argv[0];
  rb_uplo = argv[1];
  rb_ap = argv[2];
  rb_afp = argv[3];
  rb_equed = argv[4];
  rb_s = argv[5];
  rb_b = argv[6];

  uplo = StringValueCStr(rb_uplo)[0];
  equed = StringValueCStr(rb_equed)[0];
  if (!NA_IsNArray(rb_b))
    rb_raise(rb_eArgError, "b (7th argument) must be NArray");
  if (NA_RANK(rb_b) != 2)
    rb_raise(rb_eArgError, "rank of b (7th argument) must be %d", 2);
  nrhs = NA_SHAPE1(rb_b);
  ldb = NA_SHAPE0(rb_b);
  if (NA_TYPE(rb_b) != NA_DFLOAT)
    rb_b = na_change_type(rb_b, NA_DFLOAT);
  b = NA_PTR_TYPE(rb_b, doublereal*);
  fact = StringValueCStr(rb_fact)[0];
  if (!NA_IsNArray(rb_s))
    rb_raise(rb_eArgError, "s (6th argument) must be NArray");
  if (NA_RANK(rb_s) != 1)
    rb_raise(rb_eArgError, "rank of s (6th argument) must be %d", 1);
  n = NA_SHAPE0(rb_s);
  if (NA_TYPE(rb_s) != NA_DFLOAT)
    rb_s = na_change_type(rb_s, NA_DFLOAT);
  s = NA_PTR_TYPE(rb_s, doublereal*);
  if (!NA_IsNArray(rb_ap))
    rb_raise(rb_eArgError, "ap (3th argument) must be NArray");
  if (NA_RANK(rb_ap) != 1)
    rb_raise(rb_eArgError, "rank of ap (3th argument) must be %d", 1);
  if (NA_SHAPE0(rb_ap) != (n*(n+1)/2))
    rb_raise(rb_eRuntimeError, "shape 0 of ap must be %d", n*(n+1)/2);
  if (NA_TYPE(rb_ap) != NA_DFLOAT)
    rb_ap = na_change_type(rb_ap, NA_DFLOAT);
  ap = NA_PTR_TYPE(rb_ap, doublereal*);
  if (!NA_IsNArray(rb_afp))
    rb_raise(rb_eArgError, "afp (4th argument) must be NArray");
  if (NA_RANK(rb_afp) != 1)
    rb_raise(rb_eArgError, "rank of afp (4th argument) must be %d", 1);
  if (NA_SHAPE0(rb_afp) != (n*(n+1)/2))
    rb_raise(rb_eRuntimeError, "shape 0 of afp must be %d", n*(n+1)/2);
  if (NA_TYPE(rb_afp) != NA_DFLOAT)
    rb_afp = na_change_type(rb_afp, NA_DFLOAT);
  afp = NA_PTR_TYPE(rb_afp, doublereal*);
  ldx = MAX(1,n);
  {
    int shape[2];
    shape[0] = ldx;
    shape[1] = nrhs;
    rb_x = na_make_object(NA_DFLOAT, 2, shape, cNArray);
  }
  x = NA_PTR_TYPE(rb_x, doublereal*);
  {
    int shape[1];
    shape[0] = nrhs;
    rb_ferr = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  ferr = NA_PTR_TYPE(rb_ferr, doublereal*);
  {
    int shape[1];
    shape[0] = nrhs;
    rb_berr = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  berr = NA_PTR_TYPE(rb_berr, doublereal*);
  {
    int shape[1];
    shape[0] = n*(n+1)/2;
    rb_ap_out__ = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  ap_out__ = NA_PTR_TYPE(rb_ap_out__, doublereal*);
  MEMCPY(ap_out__, ap, doublereal, NA_TOTAL(rb_ap));
  rb_ap = rb_ap_out__;
  ap = ap_out__;
  {
    int shape[1];
    shape[0] = n*(n+1)/2;
    rb_afp_out__ = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  afp_out__ = NA_PTR_TYPE(rb_afp_out__, doublereal*);
  MEMCPY(afp_out__, afp, doublereal, NA_TOTAL(rb_afp));
  rb_afp = rb_afp_out__;
  afp = afp_out__;
  {
    int shape[1];
    shape[0] = n;
    rb_s_out__ = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  s_out__ = NA_PTR_TYPE(rb_s_out__, doublereal*);
  MEMCPY(s_out__, s, doublereal, NA_TOTAL(rb_s));
  rb_s = rb_s_out__;
  s = s_out__;
  {
    int shape[2];
    shape[0] = ldb;
    shape[1] = nrhs;
    rb_b_out__ = na_make_object(NA_DFLOAT, 2, shape, cNArray);
  }
  b_out__ = NA_PTR_TYPE(rb_b_out__, doublereal*);
  MEMCPY(b_out__, b, doublereal, NA_TOTAL(rb_b));
  rb_b = rb_b_out__;
  b = b_out__;
  work = ALLOC_N(doublereal, (3*n));
  iwork = ALLOC_N(integer, (n));

  dppsvx_(&fact, &uplo, &n, &nrhs, ap, afp, &equed, s, b, &ldb, x, &ldx, &rcond, ferr, berr, work, iwork, &info);

  free(work);
  free(iwork);
  rb_rcond = rb_float_new((double)rcond);
  rb_info = INT2NUM(info);
  rb_equed = rb_str_new(&equed,1);
  return rb_ary_new3(10, rb_x, rb_rcond, rb_ferr, rb_berr, rb_info, rb_ap, rb_afp, rb_equed, rb_s, rb_b);
}

void
init_lapack_dppsvx(VALUE mLapack){
  rb_define_module_function(mLapack, "dppsvx", rb_dppsvx, -1);
}
