#include "rb_lapack.h"

extern VOID ctgsen_(integer *ijob, logical *wantq, logical *wantz, logical *select, integer *n, complex *a, integer *lda, complex *b, integer *ldb, complex *alpha, complex *beta, complex *q, integer *ldq, complex *z, integer *ldz, integer *m, real *pl, real *pr, real *dif, complex *work, integer *lwork, integer *iwork, integer *liwork, integer *info);

static VALUE sHelp, sUsage;

static VALUE
rblapack_ctgsen(int argc, VALUE *argv, VALUE self){
  VALUE rblapack_ijob;
  integer ijob; 
  VALUE rblapack_wantq;
  logical wantq; 
  VALUE rblapack_wantz;
  logical wantz; 
  VALUE rblapack_select;
  logical *select; 
  VALUE rblapack_a;
  complex *a; 
  VALUE rblapack_b;
  complex *b; 
  VALUE rblapack_q;
  complex *q; 
  VALUE rblapack_z;
  complex *z; 
  VALUE rblapack_lwork;
  integer lwork; 
  VALUE rblapack_liwork;
  integer liwork; 
  VALUE rblapack_alpha;
  complex *alpha; 
  VALUE rblapack_beta;
  complex *beta; 
  VALUE rblapack_m;
  integer m; 
  VALUE rblapack_pl;
  real pl; 
  VALUE rblapack_pr;
  real pr; 
  VALUE rblapack_dif;
  real *dif; 
  VALUE rblapack_work;
  complex *work; 
  VALUE rblapack_iwork;
  integer *iwork; 
  VALUE rblapack_info;
  integer info; 
  VALUE rblapack_a_out__;
  complex *a_out__;
  VALUE rblapack_b_out__;
  complex *b_out__;
  VALUE rblapack_q_out__;
  complex *q_out__;
  VALUE rblapack_z_out__;
  complex *z_out__;

  integer n;
  integer lda;
  integer ldb;
  integer ldq;
  integer ldz;

  VALUE rb_options;
  if (argc > 0 && TYPE(argv[argc-1]) == T_HASH) {
    argc--;
    rb_options = argv[argc];
    if (rb_hash_aref(rb_options, sHelp) == Qtrue) {
      printf("%s\n", "USAGE:\n  alpha, beta, m, pl, pr, dif, work, iwork, info, a, b, q, z = NumRu::Lapack.ctgsen( ijob, wantq, wantz, select, a, b, q, z, lwork, liwork, [:usage => usage, :help => help])\n\n\nFORTRAN MANUAL\n      SUBROUTINE CTGSEN( IJOB, WANTQ, WANTZ, SELECT, N, A, LDA, B, LDB, ALPHA, BETA, Q, LDQ, Z, LDZ, M, PL, PR, DIF, WORK, LWORK, IWORK, LIWORK, INFO )\n\n*  Purpose\n*  =======\n*\n*  CTGSEN reorders the generalized Schur decomposition of a complex\n*  matrix pair (A, B) (in terms of an unitary equivalence trans-\n*  formation Q' * (A, B) * Z), so that a selected cluster of eigenvalues\n*  appears in the leading diagonal blocks of the pair (A,B). The leading\n*  columns of Q and Z form unitary bases of the corresponding left and\n*  right eigenspaces (deflating subspaces). (A, B) must be in\n*  generalized Schur canonical form, that is, A and B are both upper\n*  triangular.\n*\n*  CTGSEN also computes the generalized eigenvalues\n*\n*           w(j)= ALPHA(j) / BETA(j)\n*\n*  of the reordered matrix pair (A, B).\n*\n*  Optionally, the routine computes estimates of reciprocal condition\n*  numbers for eigenvalues and eigenspaces. These are Difu[(A11,B11),\n*  (A22,B22)] and Difl[(A11,B11), (A22,B22)], i.e. the separation(s)\n*  between the matrix pairs (A11, B11) and (A22,B22) that correspond to\n*  the selected cluster and the eigenvalues outside the cluster, resp.,\n*  and norms of \"projections\" onto left and right eigenspaces w.r.t.\n*  the selected cluster in the (1,1)-block.\n*\n*\n\n*  Arguments\n*  =========\n*\n*  IJOB    (input) integer\n*          Specifies whether condition numbers are required for the\n*          cluster of eigenvalues (PL and PR) or the deflating subspaces\n*          (Difu and Difl):\n*           =0: Only reorder w.r.t. SELECT. No extras.\n*           =1: Reciprocal of norms of \"projections\" onto left and right\n*               eigenspaces w.r.t. the selected cluster (PL and PR).\n*           =2: Upper bounds on Difu and Difl. F-norm-based estimate\n*               (DIF(1:2)).\n*           =3: Estimate of Difu and Difl. 1-norm-based estimate\n*               (DIF(1:2)).\n*               About 5 times as expensive as IJOB = 2.\n*           =4: Compute PL, PR and DIF (i.e. 0, 1 and 2 above): Economic\n*               version to get it all.\n*           =5: Compute PL, PR and DIF (i.e. 0, 1 and 3 above)\n*\n*  WANTQ   (input) LOGICAL\n*          .TRUE. : update the left transformation matrix Q;\n*          .FALSE.: do not update Q.\n*\n*  WANTZ   (input) LOGICAL\n*          .TRUE. : update the right transformation matrix Z;\n*          .FALSE.: do not update Z.\n*\n*  SELECT  (input) LOGICAL array, dimension (N)\n*          SELECT specifies the eigenvalues in the selected cluster. To\n*          select an eigenvalue w(j), SELECT(j) must be set to\n*          .TRUE..\n*\n*  N       (input) INTEGER\n*          The order of the matrices A and B. N >= 0.\n*\n*  A       (input/output) COMPLEX array, dimension(LDA,N)\n*          On entry, the upper triangular matrix A, in generalized\n*          Schur canonical form.\n*          On exit, A is overwritten by the reordered matrix A.\n*\n*  LDA     (input) INTEGER\n*          The leading dimension of the array A. LDA >= max(1,N).\n*\n*  B       (input/output) COMPLEX array, dimension(LDB,N)\n*          On entry, the upper triangular matrix B, in generalized\n*          Schur canonical form.\n*          On exit, B is overwritten by the reordered matrix B.\n*\n*  LDB     (input) INTEGER\n*          The leading dimension of the array B. LDB >= max(1,N).\n*\n*  ALPHA   (output) COMPLEX array, dimension (N)\n*  BETA    (output) COMPLEX array, dimension (N)\n*          The diagonal elements of A and B, respectively,\n*          when the pair (A,B) has been reduced to generalized Schur\n*          form.  ALPHA(i)/BETA(i) i=1,...,N are the generalized\n*          eigenvalues.\n*\n*  Q       (input/output) COMPLEX array, dimension (LDQ,N)\n*          On entry, if WANTQ = .TRUE., Q is an N-by-N matrix.\n*          On exit, Q has been postmultiplied by the left unitary\n*          transformation matrix which reorder (A, B); The leading M\n*          columns of Q form orthonormal bases for the specified pair of\n*          left eigenspaces (deflating subspaces).\n*          If WANTQ = .FALSE., Q is not referenced.\n*\n*  LDQ     (input) INTEGER\n*          The leading dimension of the array Q. LDQ >= 1.\n*          If WANTQ = .TRUE., LDQ >= N.\n*\n*  Z       (input/output) COMPLEX array, dimension (LDZ,N)\n*          On entry, if WANTZ = .TRUE., Z is an N-by-N matrix.\n*          On exit, Z has been postmultiplied by the left unitary\n*          transformation matrix which reorder (A, B); The leading M\n*          columns of Z form orthonormal bases for the specified pair of\n*          left eigenspaces (deflating subspaces).\n*          If WANTZ = .FALSE., Z is not referenced.\n*\n*  LDZ     (input) INTEGER\n*          The leading dimension of the array Z. LDZ >= 1.\n*          If WANTZ = .TRUE., LDZ >= N.\n*\n*  M       (output) INTEGER\n*          The dimension of the specified pair of left and right\n*          eigenspaces, (deflating subspaces) 0 <= M <= N.\n*\n*  PL	   (output) REAL\n*  PR	   (output) REAL\n*          If IJOB = 1, 4 or 5, PL, PR are lower bounds on the\n*          reciprocal  of the norm of \"projections\" onto left and right\n*          eigenspace with respect to the selected cluster.\n*          0 < PL, PR <= 1.\n*          If M = 0 or M = N, PL = PR  = 1.\n*          If IJOB = 0, 2 or 3 PL, PR are not referenced.\n*\n*  DIF     (output) REAL array, dimension (2).\n*          If IJOB >= 2, DIF(1:2) store the estimates of Difu and Difl.\n*          If IJOB = 2 or 4, DIF(1:2) are F-norm-based upper bounds on\n*          Difu and Difl. If IJOB = 3 or 5, DIF(1:2) are 1-norm-based\n*          estimates of Difu and Difl, computed using reversed\n*          communication with CLACN2.\n*          If M = 0 or N, DIF(1:2) = F-norm([A, B]).\n*          If IJOB = 0 or 1, DIF is not referenced.\n*\n*  WORK    (workspace/output) COMPLEX array, dimension (MAX(1,LWORK))\n*          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.\n*\n*  LWORK   (input) INTEGER\n*          The dimension of the array WORK. LWORK >=  1\n*          If IJOB = 1, 2 or 4, LWORK >=  2*M*(N-M)\n*          If IJOB = 3 or 5, LWORK >=  4*M*(N-M)\n*\n*          If LWORK = -1, then a workspace query is assumed; the routine\n*          only calculates the optimal size of the WORK array, returns\n*          this value as the first entry of the WORK array, and no error\n*          message related to LWORK is issued by XERBLA.\n*\n*  IWORK   (workspace/output) INTEGER array, dimension (MAX(1,LIWORK))\n*          On exit, if INFO = 0, IWORK(1) returns the optimal LIWORK.\n*\n*  LIWORK  (input) INTEGER\n*          The dimension of the array IWORK. LIWORK >= 1.\n*          If IJOB = 1, 2 or 4, LIWORK >=  N+2;\n*          If IJOB = 3 or 5, LIWORK >= MAX(N+2, 2*M*(N-M));\n*\n*          If LIWORK = -1, then a workspace query is assumed; the\n*          routine only calculates the optimal size of the IWORK array,\n*          returns this value as the first entry of the IWORK array, and\n*          no error message related to LIWORK is issued by XERBLA.\n*\n*  INFO    (output) INTEGER\n*            =0: Successful exit.\n*            <0: If INFO = -i, the i-th argument had an illegal value.\n*            =1: Reordering of (A, B) failed because the transformed\n*                matrix pair (A, B) would be too far from generalized\n*                Schur form; the problem is very ill-conditioned.\n*                (A, B) may have been partially reordered.\n*                If requested, 0 is returned in DIF(*), PL and PR.\n*\n*\n\n*  Further Details\n*  ===============\n*\n*  CTGSEN first collects the selected eigenvalues by computing unitary\n*  U and W that move them to the top left corner of (A, B). In other\n*  words, the selected eigenvalues are the eigenvalues of (A11, B11) in\n*\n*                U'*(A, B)*W = (A11 A12) (B11 B12) n1\n*                              ( 0  A22),( 0  B22) n2\n*                                n1  n2    n1  n2\n*\n*  where N = n1+n2 and U' means the conjugate transpose of U. The first\n*  n1 columns of U and W span the specified pair of left and right\n*  eigenspaces (deflating subspaces) of (A, B).\n*\n*  If (A, B) has been obtained from the generalized real Schur\n*  decomposition of a matrix pair (C, D) = Q*(A, B)*Z', then the\n*  reordered generalized Schur form of (C, D) is given by\n*\n*           (C, D) = (Q*U)*(U'*(A, B)*W)*(Z*W)',\n*\n*  and the first n1 columns of Q*U and Z*W span the corresponding\n*  deflating subspaces of (C, D) (Q and Z store Q*U and Z*W, resp.).\n*\n*  Note that if the selected eigenvalue is sufficiently ill-conditioned,\n*  then its value may differ significantly from its value before\n*  reordering.\n*\n*  The reciprocal condition numbers of the left and right eigenspaces\n*  spanned by the first n1 columns of U and W (or Q*U and Z*W) may\n*  be returned in DIF(1:2), corresponding to Difu and Difl, resp.\n*\n*  The Difu and Difl are defined as:\n*\n*       Difu[(A11, B11), (A22, B22)] = sigma-min( Zu )\n*  and\n*       Difl[(A11, B11), (A22, B22)] = Difu[(A22, B22), (A11, B11)],\n*\n*  where sigma-min(Zu) is the smallest singular value of the\n*  (2*n1*n2)-by-(2*n1*n2) matrix\n*\n*       Zu = [ kron(In2, A11)  -kron(A22', In1) ]\n*            [ kron(In2, B11)  -kron(B22', In1) ].\n*\n*  Here, Inx is the identity matrix of size nx and A22' is the\n*  transpose of A22. kron(X, Y) is the Kronecker product between\n*  the matrices X and Y.\n*\n*  When DIF(2) is small, small changes in (A, B) can cause large changes\n*  in the deflating subspace. An approximate (asymptotic) bound on the\n*  maximum angular error in the computed deflating subspaces is\n*\n*       EPS * norm((A, B)) / DIF(2),\n*\n*  where EPS is the machine precision.\n*\n*  The reciprocal norm of the projectors on the left and right\n*  eigenspaces associated with (A11, B11) may be returned in PL and PR.\n*  They are computed as follows. First we compute L and R so that\n*  P*(A, B)*Q is block diagonal, where\n*\n*       P = ( I -L ) n1           Q = ( I R ) n1\n*           ( 0  I ) n2    and        ( 0 I ) n2\n*             n1 n2                    n1 n2\n*\n*  and (L, R) is the solution to the generalized Sylvester equation\n*\n*       A11*R - L*A22 = -A12\n*       B11*R - L*B22 = -B12\n*\n*  Then PL = (F-norm(L)**2+1)**(-1/2) and PR = (F-norm(R)**2+1)**(-1/2).\n*  An approximate (asymptotic) bound on the average absolute error of\n*  the selected eigenvalues is\n*\n*       EPS * norm((A, B)) / PL.\n*\n*  There are also global error bounds which valid for perturbations up\n*  to a certain restriction:  A lower bound (x) on the smallest\n*  F-norm(E,F) for which an eigenvalue of (A11, B11) may move and\n*  coalesce with an eigenvalue of (A22, B22) under perturbation (E,F),\n*  (i.e. (A + E, B + F), is\n*\n*   x = min(Difu,Difl)/((1/(PL*PL)+1/(PR*PR))**(1/2)+2*max(1/PL,1/PR)).\n*\n*  An approximate bound on x can be computed from DIF(1:2), PL and PR.\n*\n*  If y = ( F-norm(E,F) / x) <= 1, the angles between the perturbed\n*  (L', R') and unperturbed (L, R) left and right deflating subspaces\n*  associated with the selected cluster in the (1,1)-blocks can be\n*  bounded as\n*\n*   max-angle(L, L') <= arctan( y * PL / (1 - y * (1 - PL * PL)**(1/2))\n*   max-angle(R, R') <= arctan( y * PR / (1 - y * (1 - PR * PR)**(1/2))\n*\n*  See LAPACK User's Guide section 4.11 or the following references\n*  for more information.\n*\n*  Note that if the default method for computing the Frobenius-norm-\n*  based estimate DIF is not wanted (see CLATDF), then the parameter\n*  IDIFJB (see below) should be changed from 3 to 4 (routine CLATDF\n*  (IJOB = 2 will be used)). See CTGSYL for more details.\n*\n*  Based on contributions by\n*     Bo Kagstrom and Peter Poromaa, Department of Computing Science,\n*     Umea University, S-901 87 Umea, Sweden.\n*\n*  References\n*  ==========\n*\n*  [1] B. Kagstrom; A Direct Method for Reordering Eigenvalues in the\n*      Generalized Real Schur Form of a Regular Matrix Pair (A, B), in\n*      M.S. Moonen et al (eds), Linear Algebra for Large Scale and\n*      Real-Time Applications, Kluwer Academic Publ. 1993, pp 195-218.\n*\n*  [2] B. Kagstrom and P. Poromaa; Computing Eigenspaces with Specified\n*      Eigenvalues of a Regular Matrix Pair (A, B) and Condition\n*      Estimation: Theory, Algorithms and Software, Report\n*      UMINF - 94.04, Department of Computing Science, Umea University,\n*      S-901 87 Umea, Sweden, 1994. Also as LAPACK Working Note 87.\n*      To appear in Numerical Algorithms, 1996.\n*\n*  [3] B. Kagstrom and P. Poromaa, LAPACK-Style Algorithms and Software\n*      for Solving the Generalized Sylvester Equation and Estimating the\n*      Separation between Regular Matrix Pairs, Report UMINF - 93.23,\n*      Department of Computing Science, Umea University, S-901 87 Umea,\n*      Sweden, December 1993, Revised April 1994, Also as LAPACK working\n*      Note 75. To appear in ACM Trans. on Math. Software, Vol 22, No 1,\n*      1996.\n*\n*  =====================================================================\n*\n\n");
      rb_exit(0);
    }
    if (rb_hash_aref(rb_options, sUsage) == Qtrue) {
      printf("%s\n", "USAGE:\n  alpha, beta, m, pl, pr, dif, work, iwork, info, a, b, q, z = NumRu::Lapack.ctgsen( ijob, wantq, wantz, select, a, b, q, z, lwork, liwork, [:usage => usage, :help => help])\n");
      rb_exit(0);
    } 
  } else
    rb_options = Qnil;
  if (argc != 10)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 10)", argc);
  rblapack_ijob = argv[0];
  rblapack_wantq = argv[1];
  rblapack_wantz = argv[2];
  rblapack_select = argv[3];
  rblapack_a = argv[4];
  rblapack_b = argv[5];
  rblapack_q = argv[6];
  rblapack_z = argv[7];
  rblapack_lwork = argv[8];
  rblapack_liwork = argv[9];
  if (rb_options != Qnil) {
  }

  ijob = NUM2INT(rblapack_ijob);
  if (!NA_IsNArray(rblapack_a))
    rb_raise(rb_eArgError, "a (5th argument) must be NArray");
  if (NA_RANK(rblapack_a) != 2)
    rb_raise(rb_eArgError, "rank of a (5th argument) must be %d", 2);
  n = NA_SHAPE1(rblapack_a);
  lda = NA_SHAPE0(rblapack_a);
  if (NA_TYPE(rblapack_a) != NA_SCOMPLEX)
    rblapack_a = na_change_type(rblapack_a, NA_SCOMPLEX);
  a = NA_PTR_TYPE(rblapack_a, complex*);
  wantz = (rblapack_wantz == Qtrue);
  if (!NA_IsNArray(rblapack_b))
    rb_raise(rb_eArgError, "b (6th argument) must be NArray");
  if (NA_RANK(rblapack_b) != 2)
    rb_raise(rb_eArgError, "rank of b (6th argument) must be %d", 2);
  if (NA_SHAPE1(rblapack_b) != n)
    rb_raise(rb_eRuntimeError, "shape 1 of b must be the same as shape 1 of a");
  ldb = NA_SHAPE0(rblapack_b);
  if (NA_TYPE(rblapack_b) != NA_SCOMPLEX)
    rblapack_b = na_change_type(rblapack_b, NA_SCOMPLEX);
  b = NA_PTR_TYPE(rblapack_b, complex*);
  liwork = NUM2INT(rblapack_liwork);
  wantq = (rblapack_wantq == Qtrue);
  if (!NA_IsNArray(rblapack_z))
    rb_raise(rb_eArgError, "z (8th argument) must be NArray");
  if (NA_RANK(rblapack_z) != 2)
    rb_raise(rb_eArgError, "rank of z (8th argument) must be %d", 2);
  if (NA_SHAPE1(rblapack_z) != n)
    rb_raise(rb_eRuntimeError, "shape 1 of z must be the same as shape 1 of a");
  ldz = NA_SHAPE0(rblapack_z);
  if (NA_TYPE(rblapack_z) != NA_SCOMPLEX)
    rblapack_z = na_change_type(rblapack_z, NA_SCOMPLEX);
  z = NA_PTR_TYPE(rblapack_z, complex*);
  lwork = NUM2INT(rblapack_lwork);
  if (!NA_IsNArray(rblapack_q))
    rb_raise(rb_eArgError, "q (7th argument) must be NArray");
  if (NA_RANK(rblapack_q) != 2)
    rb_raise(rb_eArgError, "rank of q (7th argument) must be %d", 2);
  if (NA_SHAPE1(rblapack_q) != n)
    rb_raise(rb_eRuntimeError, "shape 1 of q must be the same as shape 1 of a");
  ldq = NA_SHAPE0(rblapack_q);
  if (NA_TYPE(rblapack_q) != NA_SCOMPLEX)
    rblapack_q = na_change_type(rblapack_q, NA_SCOMPLEX);
  q = NA_PTR_TYPE(rblapack_q, complex*);
  if (!NA_IsNArray(rblapack_select))
    rb_raise(rb_eArgError, "select (4th argument) must be NArray");
  if (NA_RANK(rblapack_select) != 1)
    rb_raise(rb_eArgError, "rank of select (4th argument) must be %d", 1);
  if (NA_SHAPE0(rblapack_select) != n)
    rb_raise(rb_eRuntimeError, "shape 0 of select must be the same as shape 1 of a");
  if (NA_TYPE(rblapack_select) != NA_LINT)
    rblapack_select = na_change_type(rblapack_select, NA_LINT);
  select = NA_PTR_TYPE(rblapack_select, logical*);
  {
    int shape[1];
    shape[0] = n;
    rblapack_alpha = na_make_object(NA_SCOMPLEX, 1, shape, cNArray);
  }
  alpha = NA_PTR_TYPE(rblapack_alpha, complex*);
  {
    int shape[1];
    shape[0] = n;
    rblapack_beta = na_make_object(NA_SCOMPLEX, 1, shape, cNArray);
  }
  beta = NA_PTR_TYPE(rblapack_beta, complex*);
  {
    int shape[1];
    shape[0] = 2;
    rblapack_dif = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  dif = NA_PTR_TYPE(rblapack_dif, real*);
  {
    int shape[1];
    shape[0] = ijob==0 ? 0 : MAX(1,lwork);
    rblapack_work = na_make_object(NA_SCOMPLEX, 1, shape, cNArray);
  }
  work = NA_PTR_TYPE(rblapack_work, complex*);
  {
    int shape[1];
    shape[0] = ijob==0 ? 0 : MAX(1,liwork);
    rblapack_iwork = na_make_object(NA_LINT, 1, shape, cNArray);
  }
  iwork = NA_PTR_TYPE(rblapack_iwork, integer*);
  {
    int shape[2];
    shape[0] = lda;
    shape[1] = n;
    rblapack_a_out__ = na_make_object(NA_SCOMPLEX, 2, shape, cNArray);
  }
  a_out__ = NA_PTR_TYPE(rblapack_a_out__, complex*);
  MEMCPY(a_out__, a, complex, NA_TOTAL(rblapack_a));
  rblapack_a = rblapack_a_out__;
  a = a_out__;
  {
    int shape[2];
    shape[0] = ldb;
    shape[1] = n;
    rblapack_b_out__ = na_make_object(NA_SCOMPLEX, 2, shape, cNArray);
  }
  b_out__ = NA_PTR_TYPE(rblapack_b_out__, complex*);
  MEMCPY(b_out__, b, complex, NA_TOTAL(rblapack_b));
  rblapack_b = rblapack_b_out__;
  b = b_out__;
  {
    int shape[2];
    shape[0] = ldq;
    shape[1] = n;
    rblapack_q_out__ = na_make_object(NA_SCOMPLEX, 2, shape, cNArray);
  }
  q_out__ = NA_PTR_TYPE(rblapack_q_out__, complex*);
  MEMCPY(q_out__, q, complex, NA_TOTAL(rblapack_q));
  rblapack_q = rblapack_q_out__;
  q = q_out__;
  {
    int shape[2];
    shape[0] = ldz;
    shape[1] = n;
    rblapack_z_out__ = na_make_object(NA_SCOMPLEX, 2, shape, cNArray);
  }
  z_out__ = NA_PTR_TYPE(rblapack_z_out__, complex*);
  MEMCPY(z_out__, z, complex, NA_TOTAL(rblapack_z));
  rblapack_z = rblapack_z_out__;
  z = z_out__;

  ctgsen_(&ijob, &wantq, &wantz, select, &n, a, &lda, b, &ldb, alpha, beta, q, &ldq, z, &ldz, &m, &pl, &pr, dif, work, &lwork, iwork, &liwork, &info);

  rblapack_m = INT2NUM(m);
  rblapack_pl = rb_float_new((double)pl);
  rblapack_pr = rb_float_new((double)pr);
  rblapack_info = INT2NUM(info);
  return rb_ary_new3(13, rblapack_alpha, rblapack_beta, rblapack_m, rblapack_pl, rblapack_pr, rblapack_dif, rblapack_work, rblapack_iwork, rblapack_info, rblapack_a, rblapack_b, rblapack_q, rblapack_z);
}

void
init_lapack_ctgsen(VALUE mLapack){
  rb_define_module_function(mLapack, "ctgsen", rblapack_ctgsen, -1);
  sHelp = ID2SYM(rb_intern("help"));
  sUsage = ID2SYM(rb_intern("usage"));
}
