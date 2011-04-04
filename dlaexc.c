#include "rb_lapack.h"

extern VOID dlaexc_(logical *wantq, integer *n, doublereal *t, integer *ldt, doublereal *q, integer *ldq, integer *j1, integer *n1, integer *n2, doublereal *work, integer *info);

static VALUE
rb_dlaexc(int argc, VALUE *argv, VALUE self){
  VALUE rb_wantq;
  logical wantq; 
  VALUE rb_t;
  doublereal *t; 
  VALUE rb_q;
  doublereal *q; 
  VALUE rb_j1;
  integer j1; 
  VALUE rb_n1;
  integer n1; 
  VALUE rb_n2;
  integer n2; 
  VALUE rb_info;
  integer info; 
  VALUE rb_t_out__;
  doublereal *t_out__;
  VALUE rb_q_out__;
  doublereal *q_out__;
  doublereal *work;

  integer ldt;
  integer n;
  integer ldq;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  info, t, q = NumRu::Lapack.dlaexc( wantq, t, q, j1, n1, n2)\n    or\n  NumRu::Lapack.dlaexc  # print help\n\n\nFORTRAN MANUAL\n      SUBROUTINE DLAEXC( WANTQ, N, T, LDT, Q, LDQ, J1, N1, N2, WORK, INFO )\n\n*  Purpose\n*  =======\n*\n*  DLAEXC swaps adjacent diagonal blocks T11 and T22 of order 1 or 2 in\n*  an upper quasi-triangular matrix T by an orthogonal similarity\n*  transformation.\n*\n*  T must be in Schur canonical form, that is, block upper triangular\n*  with 1-by-1 and 2-by-2 diagonal blocks; each 2-by-2 diagonal block\n*  has its diagonal elemnts equal and its off-diagonal elements of\n*  opposite sign.\n*\n\n*  Arguments\n*  =========\n*\n*  WANTQ   (input) LOGICAL\n*          = .TRUE. : accumulate the transformation in the matrix Q;\n*          = .FALSE.: do not accumulate the transformation.\n*\n*  N       (input) INTEGER\n*          The order of the matrix T. N >= 0.\n*\n*  T       (input/output) DOUBLE PRECISION array, dimension (LDT,N)\n*          On entry, the upper quasi-triangular matrix T, in Schur\n*          canonical form.\n*          On exit, the updated matrix T, again in Schur canonical form.\n*\n*  LDT     (input) INTEGER\n*          The leading dimension of the array T. LDT >= max(1,N).\n*\n*  Q       (input/output) DOUBLE PRECISION array, dimension (LDQ,N)\n*          On entry, if WANTQ is .TRUE., the orthogonal matrix Q.\n*          On exit, if WANTQ is .TRUE., the updated matrix Q.\n*          If WANTQ is .FALSE., Q is not referenced.\n*\n*  LDQ     (input) INTEGER\n*          The leading dimension of the array Q.\n*          LDQ >= 1; and if WANTQ is .TRUE., LDQ >= N.\n*\n*  J1      (input) INTEGER\n*          The index of the first row of the first block T11.\n*\n*  N1      (input) INTEGER\n*          The order of the first block T11. N1 = 0, 1 or 2.\n*\n*  N2      (input) INTEGER\n*          The order of the second block T22. N2 = 0, 1 or 2.\n*\n*  WORK    (workspace) DOUBLE PRECISION array, dimension (N)\n*\n*  INFO    (output) INTEGER\n*          = 0: successful exit\n*          = 1: the transformed matrix T would be too far from Schur\n*               form; the blocks are not swapped and T and Q are\n*               unchanged.\n*\n\n*  =====================================================================\n*\n\n");
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
  if (NA_TYPE(rb_q) != NA_DFLOAT)
    rb_q = na_change_type(rb_q, NA_DFLOAT);
  q = NA_PTR_TYPE(rb_q, doublereal*);
  j1 = NUM2INT(rb_j1);
  if (!NA_IsNArray(rb_t))
    rb_raise(rb_eArgError, "t (2th argument) must be NArray");
  if (NA_RANK(rb_t) != 2)
    rb_raise(rb_eArgError, "rank of t (2th argument) must be %d", 2);
  if (NA_SHAPE1(rb_t) != n)
    rb_raise(rb_eRuntimeError, "shape 1 of t must be the same as shape 1 of q");
  ldt = NA_SHAPE0(rb_t);
  if (NA_TYPE(rb_t) != NA_DFLOAT)
    rb_t = na_change_type(rb_t, NA_DFLOAT);
  t = NA_PTR_TYPE(rb_t, doublereal*);
  {
    int shape[2];
    shape[0] = ldt;
    shape[1] = n;
    rb_t_out__ = na_make_object(NA_DFLOAT, 2, shape, cNArray);
  }
  t_out__ = NA_PTR_TYPE(rb_t_out__, doublereal*);
  MEMCPY(t_out__, t, doublereal, NA_TOTAL(rb_t));
  rb_t = rb_t_out__;
  t = t_out__;
  {
    int shape[2];
    shape[0] = ldq;
    shape[1] = n;
    rb_q_out__ = na_make_object(NA_DFLOAT, 2, shape, cNArray);
  }
  q_out__ = NA_PTR_TYPE(rb_q_out__, doublereal*);
  MEMCPY(q_out__, q, doublereal, NA_TOTAL(rb_q));
  rb_q = rb_q_out__;
  q = q_out__;
  work = ALLOC_N(doublereal, (n));

  dlaexc_(&wantq, &n, t, &ldt, q, &ldq, &j1, &n1, &n2, work, &info);

  free(work);
  rb_info = INT2NUM(info);
  return rb_ary_new3(3, rb_info, rb_t, rb_q);
}

void
init_lapack_dlaexc(VALUE mLapack){
  rb_define_module_function(mLapack, "dlaexc", rb_dlaexc, -1);
}
