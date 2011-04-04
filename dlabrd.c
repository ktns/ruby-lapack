#include "rb_lapack.h"

extern VOID dlabrd_(integer *m, integer *n, integer *nb, doublereal *a, integer *lda, doublereal *d, doublereal *e, doublereal *tauq, doublereal *taup, doublereal *x, integer *ldx, doublereal *y, integer *ldy);

static VALUE
rb_dlabrd(int argc, VALUE *argv, VALUE self){
  VALUE rb_m;
  integer m; 
  VALUE rb_nb;
  integer nb; 
  VALUE rb_a;
  doublereal *a; 
  VALUE rb_d;
  doublereal *d; 
  VALUE rb_e;
  doublereal *e; 
  VALUE rb_tauq;
  doublereal *tauq; 
  VALUE rb_taup;
  doublereal *taup; 
  VALUE rb_x;
  doublereal *x; 
  VALUE rb_y;
  doublereal *y; 
  VALUE rb_a_out__;
  doublereal *a_out__;

  integer lda;
  integer n;
  integer ldx;
  integer ldy;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  d, e, tauq, taup, x, y, a = NumRu::Lapack.dlabrd( m, nb, a)\n    or\n  NumRu::Lapack.dlabrd  # print help\n\n\nFORTRAN MANUAL\n\n");
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
  if (NA_TYPE(rb_a) != NA_DFLOAT)
    rb_a = na_change_type(rb_a, NA_DFLOAT);
  a = NA_PTR_TYPE(rb_a, doublereal*);
  ldy = n;
  ldx = m;
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
    rb_tauq = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  tauq = NA_PTR_TYPE(rb_tauq, doublereal*);
  {
    int shape[1];
    shape[0] = MAX(1,nb);
    rb_taup = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  taup = NA_PTR_TYPE(rb_taup, doublereal*);
  {
    int shape[2];
    shape[0] = ldx;
    shape[1] = MAX(1,nb);
    rb_x = na_make_object(NA_DFLOAT, 2, shape, cNArray);
  }
  x = NA_PTR_TYPE(rb_x, doublereal*);
  {
    int shape[2];
    shape[0] = ldy;
    shape[1] = MAX(1,nb);
    rb_y = na_make_object(NA_DFLOAT, 2, shape, cNArray);
  }
  y = NA_PTR_TYPE(rb_y, doublereal*);
  {
    int shape[2];
    shape[0] = lda;
    shape[1] = n;
    rb_a_out__ = na_make_object(NA_DFLOAT, 2, shape, cNArray);
  }
  a_out__ = NA_PTR_TYPE(rb_a_out__, doublereal*);
  MEMCPY(a_out__, a, doublereal, NA_TOTAL(rb_a));
  rb_a = rb_a_out__;
  a = a_out__;

  dlabrd_(&m, &n, &nb, a, &lda, d, e, tauq, taup, x, &ldx, y, &ldy);

  return rb_ary_new3(7, rb_d, rb_e, rb_tauq, rb_taup, rb_x, rb_y, rb_a);
}

void
init_lapack_dlabrd(VALUE mLapack){
  rb_define_module_function(mLapack, "dlabrd", rb_dlabrd, -1);
}
