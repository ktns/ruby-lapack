#include "rb_lapack.h"

extern VOID zlaqps_(integer *m, integer *n, integer *offset, integer *nb, integer *kb, doublecomplex *a, integer *lda, integer *jpvt, doublecomplex *tau, doublereal *vn1, doublereal *vn2, doublecomplex *auxv, doublecomplex *f, integer *ldf);

static VALUE sHelp, sUsage;

static VALUE
rblapack_zlaqps(int argc, VALUE *argv, VALUE self){
  VALUE rblapack_m;
  integer m; 
  VALUE rblapack_offset;
  integer offset; 
  VALUE rblapack_a;
  doublecomplex *a; 
  VALUE rblapack_jpvt;
  integer *jpvt; 
  VALUE rblapack_vn1;
  doublereal *vn1; 
  VALUE rblapack_vn2;
  doublereal *vn2; 
  VALUE rblapack_auxv;
  doublecomplex *auxv; 
  VALUE rblapack_f;
  doublecomplex *f; 
  VALUE rblapack_kb;
  integer kb; 
  VALUE rblapack_tau;
  doublecomplex *tau; 
  VALUE rblapack_a_out__;
  doublecomplex *a_out__;
  VALUE rblapack_jpvt_out__;
  integer *jpvt_out__;
  VALUE rblapack_vn1_out__;
  doublereal *vn1_out__;
  VALUE rblapack_vn2_out__;
  doublereal *vn2_out__;
  VALUE rblapack_auxv_out__;
  doublecomplex *auxv_out__;
  VALUE rblapack_f_out__;
  doublecomplex *f_out__;

  integer lda;
  integer n;
  integer nb;
  integer ldf;

  VALUE rb_options;
  if (argc > 0 && TYPE(argv[argc-1]) == T_HASH) {
    argc--;
    rb_options = argv[argc];
    if (rb_hash_aref(rb_options, sHelp) == Qtrue) {
      printf("%s\n", "USAGE:\n  kb, tau, a, jpvt, vn1, vn2, auxv, f = NumRu::Lapack.zlaqps( m, offset, a, jpvt, vn1, vn2, auxv, f, [:usage => usage, :help => help])\n\n\nFORTRAN MANUAL\n      SUBROUTINE ZLAQPS( M, N, OFFSET, NB, KB, A, LDA, JPVT, TAU, VN1, VN2, AUXV, F, LDF )\n\n*  Purpose\n*  =======\n*\n*  ZLAQPS computes a step of QR factorization with column pivoting\n*  of a complex M-by-N matrix A by using Blas-3.  It tries to factorize\n*  NB columns from A starting from the row OFFSET+1, and updates all\n*  of the matrix with Blas-3 xGEMM.\n*\n*  In some cases, due to catastrophic cancellations, it cannot\n*  factorize NB columns.  Hence, the actual number of factorized\n*  columns is returned in KB.\n*\n*  Block A(1:OFFSET,1:N) is accordingly pivoted, but not factorized.\n*\n\n*  Arguments\n*  =========\n*\n*  M       (input) INTEGER\n*          The number of rows of the matrix A. M >= 0.\n*\n*  N       (input) INTEGER\n*          The number of columns of the matrix A. N >= 0\n*\n*  OFFSET  (input) INTEGER\n*          The number of rows of A that have been factorized in\n*          previous steps.\n*\n*  NB      (input) INTEGER\n*          The number of columns to factorize.\n*\n*  KB      (output) INTEGER\n*          The number of columns actually factorized.\n*\n*  A       (input/output) COMPLEX*16 array, dimension (LDA,N)\n*          On entry, the M-by-N matrix A.\n*          On exit, block A(OFFSET+1:M,1:KB) is the triangular\n*          factor obtained and block A(1:OFFSET,1:N) has been\n*          accordingly pivoted, but no factorized.\n*          The rest of the matrix, block A(OFFSET+1:M,KB+1:N) has\n*          been updated.\n*\n*  LDA     (input) INTEGER\n*          The leading dimension of the array A. LDA >= max(1,M).\n*\n*  JPVT    (input/output) INTEGER array, dimension (N)\n*          JPVT(I) = K <==> Column K of the full matrix A has been\n*          permuted into position I in AP.\n*\n*  TAU     (output) COMPLEX*16 array, dimension (KB)\n*          The scalar factors of the elementary reflectors.\n*\n*  VN1     (input/output) DOUBLE PRECISION array, dimension (N)\n*          The vector with the partial column norms.\n*\n*  VN2     (input/output) DOUBLE PRECISION array, dimension (N)\n*          The vector with the exact column norms.\n*\n*  AUXV    (input/output) COMPLEX*16 array, dimension (NB)\n*          Auxiliar vector.\n*\n*  F       (input/output) COMPLEX*16 array, dimension (LDF,NB)\n*          Matrix F' = L*Y'*A.\n*\n*  LDF     (input) INTEGER\n*          The leading dimension of the array F. LDF >= max(1,N).\n*\n\n*  Further Details\n*  ===============\n*\n*  Based on contributions by\n*    G. Quintana-Orti, Depto. de Informatica, Universidad Jaime I, Spain\n*    X. Sun, Computer Science Dept., Duke University, USA\n*\n*  =====================================================================\n*\n\n");
      rb_exit(0);
    }
    if (rb_hash_aref(rb_options, sUsage) == Qtrue) {
      printf("%s\n", "USAGE:\n  kb, tau, a, jpvt, vn1, vn2, auxv, f = NumRu::Lapack.zlaqps( m, offset, a, jpvt, vn1, vn2, auxv, f, [:usage => usage, :help => help])\n");
      rb_exit(0);
    } 
  } else
    rb_options = Qnil;
  if (argc != 8)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 8)", argc);
  rblapack_m = argv[0];
  rblapack_offset = argv[1];
  rblapack_a = argv[2];
  rblapack_jpvt = argv[3];
  rblapack_vn1 = argv[4];
  rblapack_vn2 = argv[5];
  rblapack_auxv = argv[6];
  rblapack_f = argv[7];
  if (rb_options != Qnil) {
  }

  if (!NA_IsNArray(rblapack_auxv))
    rb_raise(rb_eArgError, "auxv (7th argument) must be NArray");
  if (NA_RANK(rblapack_auxv) != 1)
    rb_raise(rb_eArgError, "rank of auxv (7th argument) must be %d", 1);
  nb = NA_SHAPE0(rblapack_auxv);
  if (NA_TYPE(rblapack_auxv) != NA_DCOMPLEX)
    rblapack_auxv = na_change_type(rblapack_auxv, NA_DCOMPLEX);
  auxv = NA_PTR_TYPE(rblapack_auxv, doublecomplex*);
  offset = NUM2INT(rblapack_offset);
  if (!NA_IsNArray(rblapack_a))
    rb_raise(rb_eArgError, "a (3th argument) must be NArray");
  if (NA_RANK(rblapack_a) != 2)
    rb_raise(rb_eArgError, "rank of a (3th argument) must be %d", 2);
  n = NA_SHAPE1(rblapack_a);
  lda = NA_SHAPE0(rblapack_a);
  if (NA_TYPE(rblapack_a) != NA_DCOMPLEX)
    rblapack_a = na_change_type(rblapack_a, NA_DCOMPLEX);
  a = NA_PTR_TYPE(rblapack_a, doublecomplex*);
  m = NUM2INT(rblapack_m);
  if (!NA_IsNArray(rblapack_jpvt))
    rb_raise(rb_eArgError, "jpvt (4th argument) must be NArray");
  if (NA_RANK(rblapack_jpvt) != 1)
    rb_raise(rb_eArgError, "rank of jpvt (4th argument) must be %d", 1);
  if (NA_SHAPE0(rblapack_jpvt) != n)
    rb_raise(rb_eRuntimeError, "shape 0 of jpvt must be the same as shape 1 of a");
  if (NA_TYPE(rblapack_jpvt) != NA_LINT)
    rblapack_jpvt = na_change_type(rblapack_jpvt, NA_LINT);
  jpvt = NA_PTR_TYPE(rblapack_jpvt, integer*);
  if (!NA_IsNArray(rblapack_vn2))
    rb_raise(rb_eArgError, "vn2 (6th argument) must be NArray");
  if (NA_RANK(rblapack_vn2) != 1)
    rb_raise(rb_eArgError, "rank of vn2 (6th argument) must be %d", 1);
  if (NA_SHAPE0(rblapack_vn2) != n)
    rb_raise(rb_eRuntimeError, "shape 0 of vn2 must be the same as shape 1 of a");
  if (NA_TYPE(rblapack_vn2) != NA_DFLOAT)
    rblapack_vn2 = na_change_type(rblapack_vn2, NA_DFLOAT);
  vn2 = NA_PTR_TYPE(rblapack_vn2, doublereal*);
  if (!NA_IsNArray(rblapack_f))
    rb_raise(rb_eArgError, "f (8th argument) must be NArray");
  if (NA_RANK(rblapack_f) != 2)
    rb_raise(rb_eArgError, "rank of f (8th argument) must be %d", 2);
  if (NA_SHAPE1(rblapack_f) != nb)
    rb_raise(rb_eRuntimeError, "shape 1 of f must be the same as shape 0 of auxv");
  ldf = NA_SHAPE0(rblapack_f);
  if (NA_TYPE(rblapack_f) != NA_DCOMPLEX)
    rblapack_f = na_change_type(rblapack_f, NA_DCOMPLEX);
  f = NA_PTR_TYPE(rblapack_f, doublecomplex*);
  if (!NA_IsNArray(rblapack_vn1))
    rb_raise(rb_eArgError, "vn1 (5th argument) must be NArray");
  if (NA_RANK(rblapack_vn1) != 1)
    rb_raise(rb_eArgError, "rank of vn1 (5th argument) must be %d", 1);
  if (NA_SHAPE0(rblapack_vn1) != n)
    rb_raise(rb_eRuntimeError, "shape 0 of vn1 must be the same as shape 1 of a");
  if (NA_TYPE(rblapack_vn1) != NA_DFLOAT)
    rblapack_vn1 = na_change_type(rblapack_vn1, NA_DFLOAT);
  vn1 = NA_PTR_TYPE(rblapack_vn1, doublereal*);
  kb = nb;
  {
    int shape[1];
    shape[0] = kb;
    rblapack_tau = na_make_object(NA_DCOMPLEX, 1, shape, cNArray);
  }
  tau = NA_PTR_TYPE(rblapack_tau, doublecomplex*);
  {
    int shape[2];
    shape[0] = lda;
    shape[1] = n;
    rblapack_a_out__ = na_make_object(NA_DCOMPLEX, 2, shape, cNArray);
  }
  a_out__ = NA_PTR_TYPE(rblapack_a_out__, doublecomplex*);
  MEMCPY(a_out__, a, doublecomplex, NA_TOTAL(rblapack_a));
  rblapack_a = rblapack_a_out__;
  a = a_out__;
  {
    int shape[1];
    shape[0] = n;
    rblapack_jpvt_out__ = na_make_object(NA_LINT, 1, shape, cNArray);
  }
  jpvt_out__ = NA_PTR_TYPE(rblapack_jpvt_out__, integer*);
  MEMCPY(jpvt_out__, jpvt, integer, NA_TOTAL(rblapack_jpvt));
  rblapack_jpvt = rblapack_jpvt_out__;
  jpvt = jpvt_out__;
  {
    int shape[1];
    shape[0] = n;
    rblapack_vn1_out__ = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  vn1_out__ = NA_PTR_TYPE(rblapack_vn1_out__, doublereal*);
  MEMCPY(vn1_out__, vn1, doublereal, NA_TOTAL(rblapack_vn1));
  rblapack_vn1 = rblapack_vn1_out__;
  vn1 = vn1_out__;
  {
    int shape[1];
    shape[0] = n;
    rblapack_vn2_out__ = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  vn2_out__ = NA_PTR_TYPE(rblapack_vn2_out__, doublereal*);
  MEMCPY(vn2_out__, vn2, doublereal, NA_TOTAL(rblapack_vn2));
  rblapack_vn2 = rblapack_vn2_out__;
  vn2 = vn2_out__;
  {
    int shape[1];
    shape[0] = nb;
    rblapack_auxv_out__ = na_make_object(NA_DCOMPLEX, 1, shape, cNArray);
  }
  auxv_out__ = NA_PTR_TYPE(rblapack_auxv_out__, doublecomplex*);
  MEMCPY(auxv_out__, auxv, doublecomplex, NA_TOTAL(rblapack_auxv));
  rblapack_auxv = rblapack_auxv_out__;
  auxv = auxv_out__;
  {
    int shape[2];
    shape[0] = ldf;
    shape[1] = nb;
    rblapack_f_out__ = na_make_object(NA_DCOMPLEX, 2, shape, cNArray);
  }
  f_out__ = NA_PTR_TYPE(rblapack_f_out__, doublecomplex*);
  MEMCPY(f_out__, f, doublecomplex, NA_TOTAL(rblapack_f));
  rblapack_f = rblapack_f_out__;
  f = f_out__;

  zlaqps_(&m, &n, &offset, &nb, &kb, a, &lda, jpvt, tau, vn1, vn2, auxv, f, &ldf);

  rblapack_kb = INT2NUM(kb);
  return rb_ary_new3(8, rblapack_kb, rblapack_tau, rblapack_a, rblapack_jpvt, rblapack_vn1, rblapack_vn2, rblapack_auxv, rblapack_f);
}

void
init_lapack_zlaqps(VALUE mLapack){
  rb_define_module_function(mLapack, "zlaqps", rblapack_zlaqps, -1);
  sHelp = ID2SYM(rb_intern("help"));
  sUsage = ID2SYM(rb_intern("usage"));
}
