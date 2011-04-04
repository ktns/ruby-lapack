#include "rb_lapack.h"

extern VOID dgbbrd_(char *vect, integer *m, integer *n, integer *ncc, integer *kl, integer *ku, doublereal *ab, integer *ldab, doublereal *d, doublereal *e, doublereal *q, integer *ldq, doublereal *pt, integer *ldpt, doublereal *c, integer *ldc, doublereal *work, integer *info);

static VALUE
rb_dgbbrd(int argc, VALUE *argv, VALUE self){
  VALUE rb_vect;
  char vect; 
  VALUE rb_kl;
  integer kl; 
  VALUE rb_ku;
  integer ku; 
  VALUE rb_ab;
  doublereal *ab; 
  VALUE rb_c;
  doublereal *c; 
  VALUE rb_d;
  doublereal *d; 
  VALUE rb_e;
  doublereal *e; 
  VALUE rb_q;
  doublereal *q; 
  VALUE rb_pt;
  doublereal *pt; 
  VALUE rb_info;
  integer info; 
  VALUE rb_ab_out__;
  doublereal *ab_out__;
  VALUE rb_c_out__;
  doublereal *c_out__;
  doublereal *work;

  integer ldab;
  integer n;
  integer ldc;
  integer ncc;
  integer ldq;
  integer m;
  integer ldpt;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  d, e, q, pt, info, ab, c = NumRu::Lapack.dgbbrd( vect, kl, ku, ab, c)\n    or\n  NumRu::Lapack.dgbbrd  # print help\n\n\nFORTRAN MANUAL\n\n");
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
  if (NA_TYPE(rb_ab) != NA_DFLOAT)
    rb_ab = na_change_type(rb_ab, NA_DFLOAT);
  ab = NA_PTR_TYPE(rb_ab, doublereal*);
  kl = NUM2INT(rb_kl);
  if (!NA_IsNArray(rb_c))
    rb_raise(rb_eArgError, "c (5th argument) must be NArray");
  if (NA_RANK(rb_c) != 2)
    rb_raise(rb_eArgError, "rank of c (5th argument) must be %d", 2);
  ncc = NA_SHAPE1(rb_c);
  ldc = NA_SHAPE0(rb_c);
  if (NA_TYPE(rb_c) != NA_DFLOAT)
    rb_c = na_change_type(rb_c, NA_DFLOAT);
  c = NA_PTR_TYPE(rb_c, doublereal*);
  vect = StringValueCStr(rb_vect)[0];
  ku = NUM2INT(rb_ku);
  m = ldab;
  ldpt = ((lsame_(&vect,"P")) || (lsame_(&vect,"B"))) ? MAX(1,n) : 1;
  ldq = ((lsame_(&vect,"Q")) || (lsame_(&vect,"B"))) ? MAX(1,m) : 1;
  {
    int shape[1];
    shape[0] = MIN(m,n);
    rb_d = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  d = NA_PTR_TYPE(rb_d, doublereal*);
  {
    int shape[1];
    shape[0] = MIN(m,n)-1;
    rb_e = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  e = NA_PTR_TYPE(rb_e, doublereal*);
  {
    int shape[2];
    shape[0] = ldq;
    shape[1] = m;
    rb_q = na_make_object(NA_DFLOAT, 2, shape, cNArray);
  }
  q = NA_PTR_TYPE(rb_q, doublereal*);
  {
    int shape[2];
    shape[0] = ldpt;
    shape[1] = n;
    rb_pt = na_make_object(NA_DFLOAT, 2, shape, cNArray);
  }
  pt = NA_PTR_TYPE(rb_pt, doublereal*);
  {
    int shape[2];
    shape[0] = ldab;
    shape[1] = n;
    rb_ab_out__ = na_make_object(NA_DFLOAT, 2, shape, cNArray);
  }
  ab_out__ = NA_PTR_TYPE(rb_ab_out__, doublereal*);
  MEMCPY(ab_out__, ab, doublereal, NA_TOTAL(rb_ab));
  rb_ab = rb_ab_out__;
  ab = ab_out__;
  {
    int shape[2];
    shape[0] = ldc;
    shape[1] = ncc;
    rb_c_out__ = na_make_object(NA_DFLOAT, 2, shape, cNArray);
  }
  c_out__ = NA_PTR_TYPE(rb_c_out__, doublereal*);
  MEMCPY(c_out__, c, doublereal, NA_TOTAL(rb_c));
  rb_c = rb_c_out__;
  c = c_out__;
  work = ALLOC_N(doublereal, (2*MAX(m,n)));

  dgbbrd_(&vect, &m, &n, &ncc, &kl, &ku, ab, &ldab, d, e, q, &ldq, pt, &ldpt, c, &ldc, work, &info);

  free(work);
  rb_info = INT2NUM(info);
  return rb_ary_new3(7, rb_d, rb_e, rb_q, rb_pt, rb_info, rb_ab, rb_c);
}

void
init_lapack_dgbbrd(VALUE mLapack){
  rb_define_module_function(mLapack, "dgbbrd", rb_dgbbrd, -1);
}
