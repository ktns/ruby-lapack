#include "rb_lapack.h"

extern VOID zlabrd_(integer *m, integer *n, integer *nb, doublecomplex *a, integer *lda, doublereal *d, doublereal *e, doublecomplex *tauq, doublecomplex *taup, doublecomplex *x, integer *ldx, doublecomplex *y, integer *ldy);

static VALUE
rb_zlabrd(int argc, VALUE *argv, VALUE self){
  VALUE rb_m;
  integer m; 
  VALUE rb_nb;
  integer nb; 
  VALUE rb_a;
  doublecomplex *a; 
  VALUE rb_d;
  doublereal *d; 
  VALUE rb_e;
  doublereal *e; 
  VALUE rb_tauq;
  doublecomplex *tauq; 
  VALUE rb_taup;
  doublecomplex *taup; 
  VALUE rb_x;
  doublecomplex *x; 
  VALUE rb_y;
  doublecomplex *y; 
  VALUE rb_a_out__;
  doublecomplex *a_out__;

  integer lda;
  integer n;
  integer ldx;
  integer ldy;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  d, e, tauq, taup, x, y, a = NumRu::Lapack.zlabrd( m, nb, a)\n    or\n  NumRu::Lapack.zlabrd  # print help\n\n\nFORTRAN MANUAL\n\n");
    return Qnil;
  }
  if (argc != 3)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 3)", argc);
  rb_m = argv[0];
  rb_nb = argv[1];
  rb_a = argv[2];

  nb = NUM2INT(rb_nb);
  m = NUM2INT(rb_m);
  if (!NA_IsNArray(rb_a))
    rb_raise(rb_eArgError, "a (3th argument) must be NArray");
  if (NA_RANK(rb_a) != 2)
    rb_raise(rb_eArgError, "rank of a (3th argument) must be %d", 2);
  n = NA_SHAPE1(rb_a);
  lda = NA_SHAPE0(rb_a);
  if (NA_TYPE(rb_a) != NA_DCOMPLEX)
    rb_a = na_change_type(rb_a, NA_DCOMPLEX);
  a = NA_PTR_TYPE(rb_a, doublecomplex*);
  ldy = MAX(1,n);
  ldx = MAX(1,m);
  {
    int shape[1];
    shape[0] = MAX(1,nb);
    rb_d = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  d = NA_PTR_TYPE(rb_d, doublereal*);
  {
    int shape[1];
    shape[0] = MAX(1,nb);
    rb_e = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  e = NA_PTR_TYPE(rb_e, doublereal*);
  {
    int shape[1];
    shape[0] = MAX(1,nb);
    rb_tauq = na_make_object(NA_DCOMPLEX, 1, shape, cNArray);
  }
  tauq = NA_PTR_TYPE(rb_tauq, doublecomplex*);
  {
    int shape[1];
    shape[0] = MAX(1,nb);
    rb_taup = na_make_object(NA_DCOMPLEX, 1, shape, cNArray);
  }
  taup = NA_PTR_TYPE(rb_taup, doublecomplex*);
  {
    int shape[2];
    shape[0] = ldx;
    shape[1] = MAX(1,nb);
    rb_x = na_make_object(NA_DCOMPLEX, 2, shape, cNArray);
  }
  x = NA_PTR_TYPE(rb_x, doublecomplex*);
  {
    int shape[2];
    shape[0] = ldy;
    shape[1] = MAX(1,nb);
    rb_y = na_make_object(NA_DCOMPLEX, 2, shape, cNArray);
  }
  y = NA_PTR_TYPE(rb_y, doublecomplex*);
  {
    int shape[2];
    shape[0] = lda;
    shape[1] = n;
    rb_a_out__ = na_make_object(NA_DCOMPLEX, 2, shape, cNArray);
  }
  a_out__ = NA_PTR_TYPE(rb_a_out__, doublecomplex*);
  MEMCPY(a_out__, a, doublecomplex, NA_TOTAL(rb_a));
  rb_a = rb_a_out__;
  a = a_out__;

  zlabrd_(&m, &n, &nb, a, &lda, d, e, tauq, taup, x, &ldx, y, &ldy);

  return rb_ary_new3(7, rb_d, rb_e, rb_tauq, rb_taup, rb_x, rb_y, rb_a);
}

void
init_lapack_zlabrd(VALUE mLapack){
  rb_define_module_function(mLapack, "zlabrd", rb_zlabrd, -1);
}
