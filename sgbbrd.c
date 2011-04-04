#include "rb_lapack.h"

extern VOID sgbbrd_(char *vect, integer *m, integer *n, integer *ncc, integer *kl, integer *ku, real *ab, integer *ldab, real *d, real *e, real *q, integer *ldq, real *pt, integer *ldpt, real *c, integer *ldc, real *work, integer *info);

static VALUE
rb_sgbbrd(int argc, VALUE *argv, VALUE self){
  VALUE rb_vect;
  char vect; 
  VALUE rb_kl;
  integer kl; 
  VALUE rb_ku;
  integer ku; 
  VALUE rb_ab;
  real *ab; 
  VALUE rb_c;
  real *c; 
  VALUE rb_d;
  real *d; 
  VALUE rb_e;
  real *e; 
  VALUE rb_q;
  real *q; 
  VALUE rb_pt;
  real *pt; 
  VALUE rb_info;
  integer info; 
  VALUE rb_ab_out__;
  real *ab_out__;
  VALUE rb_c_out__;
  real *c_out__;
  real *work;

  integer ldab;
  integer n;
  integer ldc;
  integer ncc;
  integer ldq;
  integer m;
  integer ldpt;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  d, e, q, pt, info, ab, c = NumRu::Lapack.sgbbrd( vect, kl, ku, ab, c)\n    or\n  NumRu::Lapack.sgbbrd  # print help\n\n\nFORTRAN MANUAL\n\n");
    return Qnil;
  }
  if (argc != 5)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 5)", argc);
  rb_vect = argv[0];
  rb_kl = argv[1];
  rb_ku = argv[2];
  rb_ab = argv[3];
  rb_c = argv[4];

  if (!NA_IsNArray(rb_ab))
    rb_raise(rb_eArgError, "ab (4th argument) must be NArray");
  if (NA_RANK(rb_ab) != 2)
    rb_raise(rb_eArgError, "rank of ab (4th argument) must be %d", 2);
  n = NA_SHAPE1(rb_ab);
  ldab = NA_SHAPE0(rb_ab);
  if (NA_TYPE(rb_ab) != NA_SFLOAT)
    rb_ab = na_change_type(rb_ab, NA_SFLOAT);
  ab = NA_PTR_TYPE(rb_ab, real*);
  kl = NUM2INT(rb_kl);
  if (!NA_IsNArray(rb_c))
    rb_raise(rb_eArgError, "c (5th argument) must be NArray");
  if (NA_RANK(rb_c) != 2)
    rb_raise(rb_eArgError, "rank of c (5th argument) must be %d", 2);
  ncc = NA_SHAPE1(rb_c);
  ldc = NA_SHAPE0(rb_c);
  if (NA_TYPE(rb_c) != NA_SFLOAT)
    rb_c = na_change_type(rb_c, NA_SFLOAT);
  c = NA_PTR_TYPE(rb_c, real*);
  vect = StringValueCStr(rb_vect)[0];
  ku = NUM2INT(rb_ku);
  m = ldab;
  ldpt = ((lsame_(&vect,"P")) || (lsame_(&vect,"B"))) ? MAX(1,n) : 1;
  ldq = ((lsame_(&vect,"Q")) || (lsame_(&vect,"B"))) ? MAX(1,m) : 1;
  {
    int shape[1];
    shape[0] = MIN(m,n);
    rb_d = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  d = NA_PTR_TYPE(rb_d, real*);
  {
    int shape[1];
    shape[0] = MIN(m,n)-1;
    rb_e = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  e = NA_PTR_TYPE(rb_e, real*);
  {
    int shape[2];
    shape[0] = ldq;
    shape[1] = m;
    rb_q = na_make_object(NA_SFLOAT, 2, shape, cNArray);
  }
  q = NA_PTR_TYPE(rb_q, real*);
  {
    int shape[2];
    shape[0] = ldpt;
    shape[1] = n;
    rb_pt = na_make_object(NA_SFLOAT, 2, shape, cNArray);
  }
  pt = NA_PTR_TYPE(rb_pt, real*);
  {
    int shape[2];
    shape[0] = ldab;
    shape[1] = n;
    rb_ab_out__ = na_make_object(NA_SFLOAT, 2, shape, cNArray);
  }
  ab_out__ = NA_PTR_TYPE(rb_ab_out__, real*);
  MEMCPY(ab_out__, ab, real, NA_TOTAL(rb_ab));
  rb_ab = rb_ab_out__;
  ab = ab_out__;
  {
    int shape[2];
    shape[0] = ldc;
    shape[1] = ncc;
    rb_c_out__ = na_make_object(NA_SFLOAT, 2, shape, cNArray);
  }
  c_out__ = NA_PTR_TYPE(rb_c_out__, real*);
  MEMCPY(c_out__, c, real, NA_TOTAL(rb_c));
  rb_c = rb_c_out__;
  c = c_out__;
  work = ALLOC_N(real, (2*MAX(m,n)));

  sgbbrd_(&vect, &m, &n, &ncc, &kl, &ku, ab, &ldab, d, e, q, &ldq, pt, &ldpt, c, &ldc, work, &info);

  free(work);
  rb_info = INT2NUM(info);
  return rb_ary_new3(7, rb_d, rb_e, rb_q, rb_pt, rb_info, rb_ab, rb_c);
}

void
init_lapack_sgbbrd(VALUE mLapack){
  rb_define_module_function(mLapack, "sgbbrd", rb_sgbbrd, -1);
}
