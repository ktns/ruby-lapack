#include "rb_lapack.h"

extern VOID slaexc_(logical *wantq, integer *n, real *t, integer *ldt, real *q, integer *ldq, integer *j1, integer *n1, integer *n2, real *work, integer *info);

static VALUE
rb_slaexc(int argc, VALUE *argv, VALUE self){
  VALUE rb_wantq;
  logical wantq; 
  VALUE rb_t;
  real *t; 
  VALUE rb_q;
  real *q; 
  VALUE rb_j1;
  integer j1; 
  VALUE rb_n1;
  integer n1; 
  VALUE rb_n2;
  integer n2; 
  VALUE rb_info;
  integer info; 
  VALUE rb_t_out__;
  real *t_out__;
  VALUE rb_q_out__;
  real *q_out__;
  real *work;

  integer ldt;
  integer n;
  integer ldq;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  info, t, q = NumRu::Lapack.slaexc( wantq, t, q, j1, n1, n2)\n    or\n  NumRu::Lapack.slaexc  # print help\n\n\nFORTRAN MANUAL\n\n");
    return Qnil;
  }
  if (argc != 6)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 6)", argc);
  rb_wantq = argv[0];
  rb_t = argv[1];
  rb_q = argv[2];
  rb_j1 = argv[3];
  rb_n1 = argv[4];
  rb_n2 = argv[5];

  n1 = NUM2INT(rb_n1);
  wantq = (rb_wantq == Qtrue);
  n2 = NUM2INT(rb_n2);
  if (!NA_IsNArray(rb_q))
    rb_raise(rb_eArgError, "q (3th argument) must be NArray");
  if (NA_RANK(rb_q) != 2)
    rb_raise(rb_eArgError, "rank of q (3th argument) must be %d", 2);
  n = NA_SHAPE1(rb_q);
  ldq = NA_SHAPE0(rb_q);
  if (NA_TYPE(rb_q) != NA_SFLOAT)
    rb_q = na_change_type(rb_q, NA_SFLOAT);
  q = NA_PTR_TYPE(rb_q, real*);
  j1 = NUM2INT(rb_j1);
  if (!NA_IsNArray(rb_t))
    rb_raise(rb_eArgError, "t (2th argument) must be NArray");
  if (NA_RANK(rb_t) != 2)
    rb_raise(rb_eArgError, "rank of t (2th argument) must be %d", 2);
  if (NA_SHAPE1(rb_t) != n)
    rb_raise(rb_eRuntimeError, "shape 1 of t must be the same as shape 1 of q");
  ldt = NA_SHAPE0(rb_t);
  if (NA_TYPE(rb_t) != NA_SFLOAT)
    rb_t = na_change_type(rb_t, NA_SFLOAT);
  t = NA_PTR_TYPE(rb_t, real*);
  {
    int shape[2];
    shape[0] = ldt;
    shape[1] = n;
    rb_t_out__ = na_make_object(NA_SFLOAT, 2, shape, cNArray);
  }
  t_out__ = NA_PTR_TYPE(rb_t_out__, real*);
  MEMCPY(t_out__, t, real, NA_TOTAL(rb_t));
  rb_t = rb_t_out__;
  t = t_out__;
  {
    int shape[2];
    shape[0] = ldq;
    shape[1] = n;
    rb_q_out__ = na_make_object(NA_SFLOAT, 2, shape, cNArray);
  }
  q_out__ = NA_PTR_TYPE(rb_q_out__, real*);
  MEMCPY(q_out__, q, real, NA_TOTAL(rb_q));
  rb_q = rb_q_out__;
  q = q_out__;
  work = ALLOC_N(real, (n));

  slaexc_(&wantq, &n, t, &ldt, q, &ldq, &j1, &n1, &n2, work, &info);

  free(work);
  rb_info = INT2NUM(info);
  return rb_ary_new3(3, rb_info, rb_t, rb_q);
}

void
init_lapack_slaexc(VALUE mLapack){
  rb_define_module_function(mLapack, "slaexc", rb_slaexc, -1);
}
