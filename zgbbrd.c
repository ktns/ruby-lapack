#include "rb_lapack.h"

extern VOID zgbbrd_(char *vect, integer *m, integer *n, integer *ncc, integer *kl, integer *ku, doublecomplex *ab, integer *ldab, doublereal *d, doublereal *e, doublecomplex *q, integer *ldq, doublecomplex *pt, integer *ldpt, doublecomplex *c, integer *ldc, doublecomplex *work, doublereal *rwork, integer *info);

static VALUE
rb_zgbbrd(int argc, VALUE *argv, VALUE self){
  VALUE rb_vect;
  char vect; 
  VALUE rb_kl;
  integer kl; 
  VALUE rb_ku;
  integer ku; 
  VALUE rb_ab;
  doublecomplex *ab; 
  VALUE rb_c;
  doublecomplex *c; 
  VALUE rb_d;
  doublereal *d; 
  VALUE rb_e;
  doublereal *e; 
  VALUE rb_q;
  doublecomplex *q; 
  VALUE rb_pt;
  doublecomplex *pt; 
  VALUE rb_info;
  integer info; 
  VALUE rb_ab_out__;
  doublecomplex *ab_out__;
  VALUE rb_c_out__;
  doublecomplex *c_out__;
  doublecomplex *work;
  doublereal *rwork;

  integer ldab;
  integer n;
  integer ldc;
  integer ncc;
  integer ldq;
  integer m;
  integer ldpt;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  d, e, q, pt, info, ab, c = NumRu::Lapack.zgbbrd( vect, kl, ku, ab, c)\n    or\n  NumRu::Lapack.zgbbrd  # print help\n\n\nFORTRAN MANUAL\n\n");
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
  if (NA_TYPE(rb_ab) != NA_DCOMPLEX)
    rb_ab = na_change_type(rb_ab, NA_DCOMPLEX);
  ab = NA_PTR_TYPE(rb_ab, doublecomplex*);
  kl = NUM2INT(rb_kl);
  if (!NA_IsNArray(rb_c))
    rb_raise(rb_eArgError, "c (5th argument) must be NArray");
  if (NA_RANK(rb_c) != 2)
    rb_raise(rb_eArgError, "rank of c (5th argument) must be %d", 2);
  ncc = NA_SHAPE1(rb_c);
  ldc = NA_SHAPE0(rb_c);
  if (NA_TYPE(rb_c) != NA_DCOMPLEX)
    rb_c = na_change_type(rb_c, NA_DCOMPLEX);
  c = NA_PTR_TYPE(rb_c, doublecomplex*);
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
    rb_q = na_make_object(NA_DCOMPLEX, 2, shape, cNArray);
  }
  q = NA_PTR_TYPE(rb_q, doublecomplex*);
  {
    int shape[2];
    shape[0] = ldpt;
    shape[1] = n;
    rb_pt = na_make_object(NA_DCOMPLEX, 2, shape, cNArray);
  }
  pt = NA_PTR_TYPE(rb_pt, doublecomplex*);
  {
    int shape[2];
    shape[0] = ldab;
    shape[1] = n;
    rb_ab_out__ = na_make_object(NA_DCOMPLEX, 2, shape, cNArray);
  }
  ab_out__ = NA_PTR_TYPE(rb_ab_out__, doublecomplex*);
  MEMCPY(ab_out__, ab, doublecomplex, NA_TOTAL(rb_ab));
  rb_ab = rb_ab_out__;
  ab = ab_out__;
  {
    int shape[2];
    shape[0] = ldc;
    shape[1] = ncc;
    rb_c_out__ = na_make_object(NA_DCOMPLEX, 2, shape, cNArray);
  }
  c_out__ = NA_PTR_TYPE(rb_c_out__, doublecomplex*);
  MEMCPY(c_out__, c, doublecomplex, NA_TOTAL(rb_c));
  rb_c = rb_c_out__;
  c = c_out__;
  work = ALLOC_N(doublecomplex, (MAX(m,n)));
  rwork = ALLOC_N(doublereal, (MAX(m,n)));

  zgbbrd_(&vect, &m, &n, &ncc, &kl, &ku, ab, &ldab, d, e, q, &ldq, pt, &ldpt, c, &ldc, work, rwork, &info);

  free(work);
  free(rwork);
  rb_info = INT2NUM(info);
  return rb_ary_new3(7, rb_d, rb_e, rb_q, rb_pt, rb_info, rb_ab, rb_c);
}

void
init_lapack_zgbbrd(VALUE mLapack){
  rb_define_module_function(mLapack, "zgbbrd", rb_zgbbrd, -1);
}
