#include "rb_lapack.h"

extern VOID slaexc_(logical *wantq, integer *n, real *t, integer *ldt, real *q, integer *ldq, integer *j1, integer *n1, integer *n2, real *work, integer *info);

static VALUE sHelp, sUsage;

static VALUE
rblapack_slaexc(int argc, VALUE *argv, VALUE self){
  VALUE rblapack_wantq;
  logical wantq; 
  VALUE rblapack_t;
  real *t; 
  VALUE rblapack_q;
  real *q; 
  VALUE rblapack_j1;
  integer j1; 
  VALUE rblapack_n1;
  integer n1; 
  VALUE rblapack_n2;
  integer n2; 
  VALUE rblapack_info;
  integer info; 
  VALUE rblapack_t_out__;
  real *t_out__;
  VALUE rblapack_q_out__;
  real *q_out__;
  real *work;

  integer ldt;
  integer n;
  integer ldq;

  VALUE rb_options;
  if (argc > 0 && TYPE(argv[argc-1]) == T_HASH) {
    argc--;
    rb_options = argv[argc];
    if (rb_hash_aref(rb_options, sHelp) == Qtrue) {
      printf("%s\n", "USAGE:\n  info, t, q = NumRu::Lapack.slaexc( wantq, t, q, j1, n1, n2, [:usage => usage, :help => help])\n\n\nFORTRAN MANUAL\n      SUBROUTINE SLAEXC( WANTQ, N, T, LDT, Q, LDQ, J1, N1, N2, WORK, INFO )\n\n*  Purpose\n*  =======\n*\n*  SLAEXC swaps adjacent diagonal blocks T11 and T22 of order 1 or 2 in\n*  an upper quasi-triangular matrix T by an orthogonal similarity\n*  transformation.\n*\n*  T must be in Schur canonical form, that is, block upper triangular\n*  with 1-by-1 and 2-by-2 diagonal blocks; each 2-by-2 diagonal block\n*  has its diagonal elemnts equal and its off-diagonal elements of\n*  opposite sign.\n*\n\n*  Arguments\n*  =========\n*\n*  WANTQ   (input) LOGICAL\n*          = .TRUE. : accumulate the transformation in the matrix Q;\n*          = .FALSE.: do not accumulate the transformation.\n*\n*  N       (input) INTEGER\n*          The order of the matrix T. N >= 0.\n*\n*  T       (input/output) REAL array, dimension (LDT,N)\n*          On entry, the upper quasi-triangular matrix T, in Schur\n*          canonical form.\n*          On exit, the updated matrix T, again in Schur canonical form.\n*\n*  LDT     (input)  INTEGER\n*          The leading dimension of the array T. LDT >= max(1,N).\n*\n*  Q       (input/output) REAL array, dimension (LDQ,N)\n*          On entry, if WANTQ is .TRUE., the orthogonal matrix Q.\n*          On exit, if WANTQ is .TRUE., the updated matrix Q.\n*          If WANTQ is .FALSE., Q is not referenced.\n*\n*  LDQ     (input) INTEGER\n*          The leading dimension of the array Q.\n*          LDQ >= 1; and if WANTQ is .TRUE., LDQ >= N.\n*\n*  J1      (input) INTEGER\n*          The index of the first row of the first block T11.\n*\n*  N1      (input) INTEGER\n*          The order of the first block T11. N1 = 0, 1 or 2.\n*\n*  N2      (input) INTEGER\n*          The order of the second block T22. N2 = 0, 1 or 2.\n*\n*  WORK    (workspace) REAL array, dimension (N)\n*\n*  INFO    (output) INTEGER\n*          = 0: successful exit\n*          = 1: the transformed matrix T would be too far from Schur\n*               form; the blocks are not swapped and T and Q are\n*               unchanged.\n*\n\n*  =====================================================================\n*\n\n");
      rb_exit(0);
    }
    if (rb_hash_aref(rb_options, sUsage) == Qtrue) {
      printf("%s\n", "USAGE:\n  info, t, q = NumRu::Lapack.slaexc( wantq, t, q, j1, n1, n2, [:usage => usage, :help => help])\n");
      rb_exit(0);
    } 
  } else
    rb_options = Qnil;
  if (argc != 6)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 6)", argc);
  rblapack_wantq = argv[0];
  rblapack_t = argv[1];
  rblapack_q = argv[2];
  rblapack_j1 = argv[3];
  rblapack_n1 = argv[4];
  rblapack_n2 = argv[5];
  if (rb_options != Qnil) {
  }

  n1 = NUM2INT(rblapack_n1);
  wantq = (rblapack_wantq == Qtrue);
  n2 = NUM2INT(rblapack_n2);
  if (!NA_IsNArray(rblapack_q))
    rb_raise(rb_eArgError, "q (3th argument) must be NArray");
  if (NA_RANK(rblapack_q) != 2)
    rb_raise(rb_eArgError, "rank of q (3th argument) must be %d", 2);
  n = NA_SHAPE1(rblapack_q);
  ldq = NA_SHAPE0(rblapack_q);
  if (NA_TYPE(rblapack_q) != NA_SFLOAT)
    rblapack_q = na_change_type(rblapack_q, NA_SFLOAT);
  q = NA_PTR_TYPE(rblapack_q, real*);
  j1 = NUM2INT(rblapack_j1);
  if (!NA_IsNArray(rblapack_t))
    rb_raise(rb_eArgError, "t (2th argument) must be NArray");
  if (NA_RANK(rblapack_t) != 2)
    rb_raise(rb_eArgError, "rank of t (2th argument) must be %d", 2);
  if (NA_SHAPE1(rblapack_t) != n)
    rb_raise(rb_eRuntimeError, "shape 1 of t must be the same as shape 1 of q");
  ldt = NA_SHAPE0(rblapack_t);
  if (NA_TYPE(rblapack_t) != NA_SFLOAT)
    rblapack_t = na_change_type(rblapack_t, NA_SFLOAT);
  t = NA_PTR_TYPE(rblapack_t, real*);
  {
    int shape[2];
    shape[0] = ldt;
    shape[1] = n;
    rblapack_t_out__ = na_make_object(NA_SFLOAT, 2, shape, cNArray);
  }
  t_out__ = NA_PTR_TYPE(rblapack_t_out__, real*);
  MEMCPY(t_out__, t, real, NA_TOTAL(rblapack_t));
  rblapack_t = rblapack_t_out__;
  t = t_out__;
  {
    int shape[2];
    shape[0] = ldq;
    shape[1] = n;
    rblapack_q_out__ = na_make_object(NA_SFLOAT, 2, shape, cNArray);
  }
  q_out__ = NA_PTR_TYPE(rblapack_q_out__, real*);
  MEMCPY(q_out__, q, real, NA_TOTAL(rblapack_q));
  rblapack_q = rblapack_q_out__;
  q = q_out__;
  work = ALLOC_N(real, (n));

  slaexc_(&wantq, &n, t, &ldt, q, &ldq, &j1, &n1, &n2, work, &info);

  free(work);
  rblapack_info = INT2NUM(info);
  return rb_ary_new3(3, rblapack_info, rblapack_t, rblapack_q);
}

void
init_lapack_slaexc(VALUE mLapack){
  rb_define_module_function(mLapack, "slaexc", rblapack_slaexc, -1);
  sHelp = ID2SYM(rb_intern("help"));
  sUsage = ID2SYM(rb_intern("usage"));
}
