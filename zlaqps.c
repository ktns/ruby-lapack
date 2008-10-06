#include "rb_lapack.h"

static VALUE
rb_zlaqps(int argc, VALUE *argv, VALUE self){
  VALUE rb_m;
  integer m; 
  VALUE rb_offset;
  integer offset; 
  VALUE rb_a;
  doublecomplex *a; 
  VALUE rb_jpvt;
  integer *jpvt; 
  VALUE rb_vn1;
  doublereal *vn1; 
  VALUE rb_vn2;
  doublereal *vn2; 
  VALUE rb_auxv;
  doublecomplex *auxv; 
  VALUE rb_f;
  doublecomplex *f; 
  VALUE rb_kb;
  integer kb; 
  VALUE rb_tau;
  doublecomplex *tau; 
  VALUE rb_a_out__;
  doublecomplex *a_out__;
  VALUE rb_jpvt_out__;
  integer *jpvt_out__;
  VALUE rb_vn1_out__;
  doublereal *vn1_out__;
  VALUE rb_vn2_out__;
  doublereal *vn2_out__;
  VALUE rb_auxv_out__;
  doublecomplex *auxv_out__;
  VALUE rb_f_out__;
  doublecomplex *f_out__;

  integer lda;
  integer n;
  integer nb;
  integer ldf;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  kb, tau, a, jpvt, vn1, vn2, auxv, f = NumRu::Lapack.zlaqps( m, offset, a, jpvt, vn1, vn2, auxv, f)\n    or\n  NumRu::Lapack.zlaqps  # print help\n\n\nFORTRAN MANUAL\n      SUBROUTINE ZLAQPS( M, N, OFFSET, NB, KB, A, LDA, JPVT, TAU, VN1, VN2, AUXV, F, LDF )\n\n*  Purpose\n*  =======\n*\n*  ZLAQPS computes a step of QR factorization with column pivoting\n*  of a complex M-by-N matrix A by using Blas-3.  It tries to factorize\n*  NB columns from A starting from the row OFFSET+1, and updates all\n*  of the matrix with Blas-3 xGEMM.\n*\n*  In some cases, due to catastrophic cancellations, it cannot\n*  factorize NB columns.  Hence, the actual number of factorized\n*  columns is returned in KB.\n*\n*  Block A(1:OFFSET,1:N) is accordingly pivoted, but not factorized.\n*\n\n*  Arguments\n*  =========\n*\n*  M       (input) INTEGER\n*          The number of rows of the matrix A. M >= 0.\n*\n*  N       (input) INTEGER\n*          The number of columns of the matrix A. N >= 0\n*\n*  OFFSET  (input) INTEGER\n*          The number of rows of A that have been factorized in\n*          previous steps.\n*\n*  NB      (input) INTEGER\n*          The number of columns to factorize.\n*\n*  KB      (output) INTEGER\n*          The number of columns actually factorized.\n*\n*  A       (input/output) COMPLEX*16 array, dimension (LDA,N)\n*          On entry, the M-by-N matrix A.\n*          On exit, block A(OFFSET+1:M,1:KB) is the triangular\n*          factor obtained and block A(1:OFFSET,1:N) has been\n*          accordingly pivoted, but no factorized.\n*          The rest of the matrix, block A(OFFSET+1:M,KB+1:N) has\n*          been updated.\n*\n*  LDA     (input) INTEGER\n*          The leading dimension of the array A. LDA >= max(1,M).\n*\n*  JPVT    (input/output) INTEGER array, dimension (N)\n*          JPVT(I) = K <==> Column K of the full matrix A has been\n*          permuted into position I in AP.\n*\n*  TAU     (output) COMPLEX*16 array, dimension (KB)\n*          The scalar factors of the elementary reflectors.\n*\n*  VN1     (input/output) DOUBLE PRECISION array, dimension (N)\n*          The vector with the partial column norms.\n*\n*  VN2     (input/output) DOUBLE PRECISION array, dimension (N)\n*          The vector with the exact column norms.\n*\n*  AUXV    (input/output) COMPLEX*16 array, dimension (NB)\n*          Auxiliar vector.\n*\n*  F       (input/output) COMPLEX*16 array, dimension (LDF,NB)\n*          Matrix F' = L*Y'*A.\n*\n*  LDF     (input) INTEGER\n*          The leading dimension of the array F. LDF >= max(1,N).\n*\n\n*  Further Details\n*  ===============\n*\n*  Based on contributions by\n*    G. Quintana-Orti, Depto. de Informatica, Universidad Jaime I, Spain\n*    X. Sun, Computer Science Dept., Duke University, USA\n*\n*  =====================================================================\n*\n\n");
    return Qnil;
  }
  if (argc != 8)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 8)", argc);
  rb_m = argv[0];
  rb_offset = argv[1];
  rb_a = argv[2];
  rb_jpvt = argv[3];
  rb_vn1 = argv[4];
  rb_vn2 = argv[5];
  rb_auxv = argv[6];
  rb_f = argv[7];

  m = NUM2INT(rb_m);
  offset = NUM2INT(rb_offset);
  if (!NA_IsNArray(rb_a))
    rb_raise(rb_eArgError, "a (3th argument) must be NArray");
  if (NA_RANK(rb_a) != 2)
    rb_raise(rb_eArgError, "rank of a (3th argument) must be %d", 2);
  lda = NA_SHAPE0(rb_a);
  n = NA_SHAPE1(rb_a);
  if (NA_TYPE(rb_a) != NA_DCOMPLEX)
    rb_a = na_change_type(rb_a, NA_DCOMPLEX);
  a = NA_PTR_TYPE(rb_a, doublecomplex*);
  if (!NA_IsNArray(rb_jpvt))
    rb_raise(rb_eArgError, "jpvt (4th argument) must be NArray");
  if (NA_RANK(rb_jpvt) != 1)
    rb_raise(rb_eArgError, "rank of jpvt (4th argument) must be %d", 1);
  if (NA_SHAPE0(rb_jpvt) != n)
    rb_raise(rb_eRuntimeError, "shape 0 of jpvt must be the same as shape 1 of a");
  if (NA_TYPE(rb_jpvt) != NA_LINT)
    rb_jpvt = na_change_type(rb_jpvt, NA_LINT);
  jpvt = NA_PTR_TYPE(rb_jpvt, integer*);
  if (!NA_IsNArray(rb_vn1))
    rb_raise(rb_eArgError, "vn1 (5th argument) must be NArray");
  if (NA_RANK(rb_vn1) != 1)
    rb_raise(rb_eArgError, "rank of vn1 (5th argument) must be %d", 1);
  if (NA_SHAPE0(rb_vn1) != n)
    rb_raise(rb_eRuntimeError, "shape 0 of vn1 must be the same as shape 1 of a");
  if (NA_TYPE(rb_vn1) != NA_DFLOAT)
    rb_vn1 = na_change_type(rb_vn1, NA_DFLOAT);
  vn1 = NA_PTR_TYPE(rb_vn1, doublereal*);
  if (!NA_IsNArray(rb_vn2))
    rb_raise(rb_eArgError, "vn2 (6th argument) must be NArray");
  if (NA_RANK(rb_vn2) != 1)
    rb_raise(rb_eArgError, "rank of vn2 (6th argument) must be %d", 1);
  if (NA_SHAPE0(rb_vn2) != n)
    rb_raise(rb_eRuntimeError, "shape 0 of vn2 must be the same as shape 1 of a");
  if (NA_TYPE(rb_vn2) != NA_DFLOAT)
    rb_vn2 = na_change_type(rb_vn2, NA_DFLOAT);
  vn2 = NA_PTR_TYPE(rb_vn2, doublereal*);
  if (!NA_IsNArray(rb_auxv))
    rb_raise(rb_eArgError, "auxv (7th argument) must be NArray");
  if (NA_RANK(rb_auxv) != 1)
    rb_raise(rb_eArgError, "rank of auxv (7th argument) must be %d", 1);
  nb = NA_SHAPE0(rb_auxv);
  if (NA_TYPE(rb_auxv) != NA_DCOMPLEX)
    rb_auxv = na_change_type(rb_auxv, NA_DCOMPLEX);
  auxv = NA_PTR_TYPE(rb_auxv, doublecomplex*);
  if (!NA_IsNArray(rb_f))
    rb_raise(rb_eArgError, "f (8th argument) must be NArray");
  if (NA_RANK(rb_f) != 2)
    rb_raise(rb_eArgError, "rank of f (8th argument) must be %d", 2);
  ldf = NA_SHAPE0(rb_f);
  if (NA_SHAPE1(rb_f) != nb)
    rb_raise(rb_eRuntimeError, "shape 1 of f must be the same as shape 0 of auxv");
  if (NA_TYPE(rb_f) != NA_DCOMPLEX)
    rb_f = na_change_type(rb_f, NA_DCOMPLEX);
  f = NA_PTR_TYPE(rb_f, doublecomplex*);
  kb = nb;
  {
    int shape[1];
    shape[0] = kb;
    rb_tau = na_make_object(NA_DCOMPLEX, 1, shape, cNArray);
  }
  tau = NA_PTR_TYPE(rb_tau, doublecomplex*);
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
  {
    int shape[1];
    shape[0] = n;
    rb_jpvt_out__ = na_make_object(NA_LINT, 1, shape, cNArray);
  }
  jpvt_out__ = NA_PTR_TYPE(rb_jpvt_out__, integer*);
  MEMCPY(jpvt_out__, jpvt, integer, NA_TOTAL(rb_jpvt));
  rb_jpvt = rb_jpvt_out__;
  jpvt = jpvt_out__;
  {
    int shape[1];
    shape[0] = n;
    rb_vn1_out__ = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  vn1_out__ = NA_PTR_TYPE(rb_vn1_out__, doublereal*);
  MEMCPY(vn1_out__, vn1, doublereal, NA_TOTAL(rb_vn1));
  rb_vn1 = rb_vn1_out__;
  vn1 = vn1_out__;
  {
    int shape[1];
    shape[0] = n;
    rb_vn2_out__ = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  vn2_out__ = NA_PTR_TYPE(rb_vn2_out__, doublereal*);
  MEMCPY(vn2_out__, vn2, doublereal, NA_TOTAL(rb_vn2));
  rb_vn2 = rb_vn2_out__;
  vn2 = vn2_out__;
  {
    int shape[1];
    shape[0] = nb;
    rb_auxv_out__ = na_make_object(NA_DCOMPLEX, 1, shape, cNArray);
  }
  auxv_out__ = NA_PTR_TYPE(rb_auxv_out__, doublecomplex*);
  MEMCPY(auxv_out__, auxv, doublecomplex, NA_TOTAL(rb_auxv));
  rb_auxv = rb_auxv_out__;
  auxv = auxv_out__;
  {
    int shape[2];
    shape[0] = ldf;
    shape[1] = nb;
    rb_f_out__ = na_make_object(NA_DCOMPLEX, 2, shape, cNArray);
  }
  f_out__ = NA_PTR_TYPE(rb_f_out__, doublecomplex*);
  MEMCPY(f_out__, f, doublecomplex, NA_TOTAL(rb_f));
  rb_f = rb_f_out__;
  f = f_out__;

  zlaqps_(&m, &n, &offset, &nb, &kb, a, &lda, jpvt, tau, vn1, vn2, auxv, f, &ldf);

  rb_kb = INT2NUM(kb);
  return rb_ary_new3(8, rb_kb, rb_tau, rb_a, rb_jpvt, rb_vn1, rb_vn2, rb_auxv, rb_f);
}

void
init_lapack_zlaqps(VALUE mLapack){
  rb_define_module_function(mLapack, "zlaqps", rb_zlaqps, -1);
}
