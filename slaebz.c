#include "rb_lapack.h"

extern VOID slaebz_(integer *ijob, integer *nitmax, integer *n, integer *mmax, integer *minp, integer *nbmin, real *abstol, real *reltol, real *pivmin, real *d, real *e, real *e2, integer *nval, real *ab, real *c, integer *mout, integer *nab, real *work, integer *iwork, integer *info);

static VALUE
rb_slaebz(int argc, VALUE *argv, VALUE self){
  VALUE rb_ijob;
  integer ijob; 
  VALUE rb_nitmax;
  integer nitmax; 
  VALUE rb_minp;
  integer minp; 
  VALUE rb_nbmin;
  integer nbmin; 
  VALUE rb_abstol;
  real abstol; 
  VALUE rb_reltol;
  real reltol; 
  VALUE rb_pivmin;
  real pivmin; 
  VALUE rb_d;
  real *d; 
  VALUE rb_e;
  real *e; 
  VALUE rb_e2;
  real *e2; 
  VALUE rb_nval;
  integer *nval; 
  VALUE rb_ab;
  real *ab; 
  VALUE rb_c;
  real *c; 
  VALUE rb_nab;
  integer *nab; 
  VALUE rb_mout;
  integer mout; 
  VALUE rb_info;
  integer info; 
  VALUE rb_nval_out__;
  integer *nval_out__;
  VALUE rb_ab_out__;
  real *ab_out__;
  VALUE rb_c_out__;
  real *c_out__;
  VALUE rb_nab_out__;
  integer *nab_out__;
  real *work;
  integer *iwork;

  integer n;
  integer mmax;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  mout, info, nval, ab, c, nab = NumRu::Lapack.slaebz( ijob, nitmax, minp, nbmin, abstol, reltol, pivmin, d, e, e2, nval, ab, c, nab)\n    or\n  NumRu::Lapack.slaebz  # print help\n\n\nFORTRAN MANUAL\n      SUBROUTINE SLAEBZ( IJOB, NITMAX, N, MMAX, MINP, NBMIN, ABSTOL, RELTOL, PIVMIN, D, E, E2, NVAL, AB, C, MOUT, NAB, WORK, IWORK, INFO )\n\n*  Purpose\n*  =======\n*\n*  SLAEBZ contains the iteration loops which compute and use the\n*  function N(w), which is the count of eigenvalues of a symmetric\n*  tridiagonal matrix T less than or equal to its argument  w.  It\n*  performs a choice of two types of loops:\n*\n*  IJOB=1, followed by\n*  IJOB=2: It takes as input a list of intervals and returns a list of\n*          sufficiently small intervals whose union contains the same\n*          eigenvalues as the union of the original intervals.\n*          The input intervals are (AB(j,1),AB(j,2)], j=1,...,MINP.\n*          The output interval (AB(j,1),AB(j,2)] will contain\n*          eigenvalues NAB(j,1)+1,...,NAB(j,2), where 1 <= j <= MOUT.\n*\n*  IJOB=3: It performs a binary search in each input interval\n*          (AB(j,1),AB(j,2)] for a point  w(j)  such that\n*          N(w(j))=NVAL(j), and uses  C(j)  as the starting point of\n*          the search.  If such a w(j) is found, then on output\n*          AB(j,1)=AB(j,2)=w.  If no such w(j) is found, then on output\n*          (AB(j,1),AB(j,2)] will be a small interval containing the\n*          point where N(w) jumps through NVAL(j), unless that point\n*          lies outside the initial interval.\n*\n*  Note that the intervals are in all cases half-open intervals,\n*  i.e., of the form  (a,b] , which includes  b  but not  a .\n*\n*  To avoid underflow, the matrix should be scaled so that its largest\n*  element is no greater than  overflow**(1/2) * underflow**(1/4)\n*  in absolute value.  To assure the most accurate computation\n*  of small eigenvalues, the matrix should be scaled to be\n*  not much smaller than that, either.\n*\n*  See W. Kahan \"Accurate Eigenvalues of a Symmetric Tridiagonal\n*  Matrix\", Report CS41, Computer Science Dept., Stanford\n*  University, July 21, 1966\n*\n*  Note: the arguments are, in general, *not* checked for unreasonable\n*  values.\n*\n\n*  Arguments\n*  =========\n*\n*  IJOB    (input) INTEGER\n*          Specifies what is to be done:\n*          = 1:  Compute NAB for the initial intervals.\n*          = 2:  Perform bisection iteration to find eigenvalues of T.\n*          = 3:  Perform bisection iteration to invert N(w), i.e.,\n*                to find a point which has a specified number of\n*                eigenvalues of T to its left.\n*          Other values will cause SLAEBZ to return with INFO=-1.\n*\n*  NITMAX  (input) INTEGER\n*          The maximum number of \"levels\" of bisection to be\n*          performed, i.e., an interval of width W will not be made\n*          smaller than 2^(-NITMAX) * W.  If not all intervals\n*          have converged after NITMAX iterations, then INFO is set\n*          to the number of non-converged intervals.\n*\n*  N       (input) INTEGER\n*          The dimension n of the tridiagonal matrix T.  It must be at\n*          least 1.\n*\n*  MMAX    (input) INTEGER\n*          The maximum number of intervals.  If more than MMAX intervals\n*          are generated, then SLAEBZ will quit with INFO=MMAX+1.\n*\n*  MINP    (input) INTEGER\n*          The initial number of intervals.  It may not be greater than\n*          MMAX.\n*\n*  NBMIN   (input) INTEGER\n*          The smallest number of intervals that should be processed\n*          using a vector loop.  If zero, then only the scalar loop\n*          will be used.\n*\n*  ABSTOL  (input) REAL\n*          The minimum (absolute) width of an interval.  When an\n*          interval is narrower than ABSTOL, or than RELTOL times the\n*          larger (in magnitude) endpoint, then it is considered to be\n*          sufficiently small, i.e., converged.  This must be at least\n*          zero.\n*\n*  RELTOL  (input) REAL\n*          The minimum relative width of an interval.  When an interval\n*          is narrower than ABSTOL, or than RELTOL times the larger (in\n*          magnitude) endpoint, then it is considered to be\n*          sufficiently small, i.e., converged.  Note: this should\n*          always be at least radix*machine epsilon.\n*\n*  PIVMIN  (input) REAL\n*          The minimum absolute value of a \"pivot\" in the Sturm\n*          sequence loop.  This *must* be at least  max |e(j)**2| *\n*          safe_min  and at least safe_min, where safe_min is at least\n*          the smallest number that can divide one without overflow.\n*\n*  D       (input) REAL array, dimension (N)\n*          The diagonal elements of the tridiagonal matrix T.\n*\n*  E       (input) REAL array, dimension (N)\n*          The offdiagonal elements of the tridiagonal matrix T in\n*          positions 1 through N-1.  E(N) is arbitrary.\n*\n*  E2      (input) REAL array, dimension (N)\n*          The squares of the offdiagonal elements of the tridiagonal\n*          matrix T.  E2(N) is ignored.\n*\n*  NVAL    (input/output) INTEGER array, dimension (MINP)\n*          If IJOB=1 or 2, not referenced.\n*          If IJOB=3, the desired values of N(w).  The elements of NVAL\n*          will be reordered to correspond with the intervals in AB.\n*          Thus, NVAL(j) on output will not, in general be the same as\n*          NVAL(j) on input, but it will correspond with the interval\n*          (AB(j,1),AB(j,2)] on output.\n*\n*  AB      (input/output) REAL array, dimension (MMAX,2)\n*          The endpoints of the intervals.  AB(j,1) is  a(j), the left\n*          endpoint of the j-th interval, and AB(j,2) is b(j), the\n*          right endpoint of the j-th interval.  The input intervals\n*          will, in general, be modified, split, and reordered by the\n*          calculation.\n*\n*  C       (input/output) REAL array, dimension (MMAX)\n*          If IJOB=1, ignored.\n*          If IJOB=2, workspace.\n*          If IJOB=3, then on input C(j) should be initialized to the\n*          first search point in the binary search.\n*\n*  MOUT    (output) INTEGER\n*          If IJOB=1, the number of eigenvalues in the intervals.\n*          If IJOB=2 or 3, the number of intervals output.\n*          If IJOB=3, MOUT will equal MINP.\n*\n*  NAB     (input/output) INTEGER array, dimension (MMAX,2)\n*          If IJOB=1, then on output NAB(i,j) will be set to N(AB(i,j)).\n*          If IJOB=2, then on input, NAB(i,j) should be set.  It must\n*             satisfy the condition:\n*             N(AB(i,1)) <= NAB(i,1) <= NAB(i,2) <= N(AB(i,2)),\n*             which means that in interval i only eigenvalues\n*             NAB(i,1)+1,...,NAB(i,2) will be considered.  Usually,\n*             NAB(i,j)=N(AB(i,j)), from a previous call to SLAEBZ with\n*             IJOB=1.\n*             On output, NAB(i,j) will contain\n*             max(na(k),min(nb(k),N(AB(i,j)))), where k is the index of\n*             the input interval that the output interval\n*             (AB(j,1),AB(j,2)] came from, and na(k) and nb(k) are the\n*             the input values of NAB(k,1) and NAB(k,2).\n*          If IJOB=3, then on output, NAB(i,j) contains N(AB(i,j)),\n*             unless N(w) > NVAL(i) for all search points  w , in which\n*             case NAB(i,1) will not be modified, i.e., the output\n*             value will be the same as the input value (modulo\n*             reorderings -- see NVAL and AB), or unless N(w) < NVAL(i)\n*             for all search points  w , in which case NAB(i,2) will\n*             not be modified.  Normally, NAB should be set to some\n*             distinctive value(s) before SLAEBZ is called.\n*\n*  WORK    (workspace) REAL array, dimension (MMAX)\n*          Workspace.\n*\n*  IWORK   (workspace) INTEGER array, dimension (MMAX)\n*          Workspace.\n*\n*  INFO    (output) INTEGER\n*          = 0:       All intervals converged.\n*          = 1--MMAX: The last INFO intervals did not converge.\n*          = MMAX+1:  More than MMAX intervals were generated.\n*\n\n*  Further Details\n*  ===============\n*\n*      This routine is intended to be called only by other LAPACK\n*  routines, thus the interface is less user-friendly.  It is intended\n*  for two purposes:\n*\n*  (a) finding eigenvalues.  In this case, SLAEBZ should have one or\n*      more initial intervals set up in AB, and SLAEBZ should be called\n*      with IJOB=1.  This sets up NAB, and also counts the eigenvalues.\n*      Intervals with no eigenvalues would usually be thrown out at\n*      this point.  Also, if not all the eigenvalues in an interval i\n*      are desired, NAB(i,1) can be increased or NAB(i,2) decreased.\n*      For example, set NAB(i,1)=NAB(i,2)-1 to get the largest\n*      eigenvalue.  SLAEBZ is then called with IJOB=2 and MMAX\n*      no smaller than the value of MOUT returned by the call with\n*      IJOB=1.  After this (IJOB=2) call, eigenvalues NAB(i,1)+1\n*      through NAB(i,2) are approximately AB(i,1) (or AB(i,2)) to the\n*      tolerance specified by ABSTOL and RELTOL.\n*\n*  (b) finding an interval (a',b'] containing eigenvalues w(f),...,w(l).\n*      In this case, start with a Gershgorin interval  (a,b).  Set up\n*      AB to contain 2 search intervals, both initially (a,b).  One\n*      NVAL element should contain  f-1  and the other should contain  l\n*      , while C should contain a and b, resp.  NAB(i,1) should be -1\n*      and NAB(i,2) should be N+1, to flag an error if the desired\n*      interval does not lie in (a,b).  SLAEBZ is then called with\n*      IJOB=3.  On exit, if w(f-1) < w(f), then one of the intervals --\n*      j -- will have AB(j,1)=AB(j,2) and NAB(j,1)=NAB(j,2)=f-1, while\n*      if, to the specified tolerance, w(f-k)=...=w(f+r), k > 0 and r\n*      >= 0, then the interval will have  N(AB(j,1))=NAB(j,1)=f-k and\n*      N(AB(j,2))=NAB(j,2)=f+r.  The cases w(l) < w(l+1) and\n*      w(l-r)=...=w(l+k) are handled similarly.\n*\n*  =====================================================================\n*\n\n");
    return Qnil;
  }
  if (argc != 14)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 14)", argc);
  rb_ijob = argv[0];
  rb_nitmax = argv[1];
  rb_minp = argv[2];
  rb_nbmin = argv[3];
  rb_abstol = argv[4];
  rb_reltol = argv[5];
  rb_pivmin = argv[6];
  rb_d = argv[7];
  rb_e = argv[8];
  rb_e2 = argv[9];
  rb_nval = argv[10];
  rb_ab = argv[11];
  rb_c = argv[12];
  rb_nab = argv[13];

  abstol = (real)NUM2DBL(rb_abstol);
  ijob = NUM2INT(rb_ijob);
  if (!NA_IsNArray(rb_ab))
    rb_raise(rb_eArgError, "ab (12th argument) must be NArray");
  if (NA_RANK(rb_ab) != 2)
    rb_raise(rb_eArgError, "rank of ab (12th argument) must be %d", 2);
  if (NA_SHAPE1(rb_ab) != (2))
    rb_raise(rb_eRuntimeError, "shape 1 of ab must be %d", 2);
  mmax = NA_SHAPE0(rb_ab);
  if (NA_TYPE(rb_ab) != NA_SFLOAT)
    rb_ab = na_change_type(rb_ab, NA_SFLOAT);
  ab = NA_PTR_TYPE(rb_ab, real*);
  if (!NA_IsNArray(rb_e2))
    rb_raise(rb_eArgError, "e2 (10th argument) must be NArray");
  if (NA_RANK(rb_e2) != 1)
    rb_raise(rb_eArgError, "rank of e2 (10th argument) must be %d", 1);
  n = NA_SHAPE0(rb_e2);
  if (NA_TYPE(rb_e2) != NA_SFLOAT)
    rb_e2 = na_change_type(rb_e2, NA_SFLOAT);
  e2 = NA_PTR_TYPE(rb_e2, real*);
  nitmax = NUM2INT(rb_nitmax);
  pivmin = (real)NUM2DBL(rb_pivmin);
  if (!NA_IsNArray(rb_nab))
    rb_raise(rb_eArgError, "nab (14th argument) must be NArray");
  if (NA_RANK(rb_nab) != 2)
    rb_raise(rb_eArgError, "rank of nab (14th argument) must be %d", 2);
  if (NA_SHAPE1(rb_nab) != (2))
    rb_raise(rb_eRuntimeError, "shape 1 of nab must be %d", 2);
  if (NA_SHAPE0(rb_nab) != mmax)
    rb_raise(rb_eRuntimeError, "shape 0 of nab must be the same as shape 0 of ab");
  if (NA_TYPE(rb_nab) != NA_LINT)
    rb_nab = na_change_type(rb_nab, NA_LINT);
  nab = NA_PTR_TYPE(rb_nab, integer*);
  nbmin = NUM2INT(rb_nbmin);
  reltol = (real)NUM2DBL(rb_reltol);
  if (!NA_IsNArray(rb_e))
    rb_raise(rb_eArgError, "e (9th argument) must be NArray");
  if (NA_RANK(rb_e) != 1)
    rb_raise(rb_eArgError, "rank of e (9th argument) must be %d", 1);
  if (NA_SHAPE0(rb_e) != n)
    rb_raise(rb_eRuntimeError, "shape 0 of e must be the same as shape 0 of e2");
  if (NA_TYPE(rb_e) != NA_SFLOAT)
    rb_e = na_change_type(rb_e, NA_SFLOAT);
  e = NA_PTR_TYPE(rb_e, real*);
  minp = NUM2INT(rb_minp);
  if (!NA_IsNArray(rb_d))
    rb_raise(rb_eArgError, "d (8th argument) must be NArray");
  if (NA_RANK(rb_d) != 1)
    rb_raise(rb_eArgError, "rank of d (8th argument) must be %d", 1);
  if (NA_SHAPE0(rb_d) != n)
    rb_raise(rb_eRuntimeError, "shape 0 of d must be the same as shape 0 of e2");
  if (NA_TYPE(rb_d) != NA_SFLOAT)
    rb_d = na_change_type(rb_d, NA_SFLOAT);
  d = NA_PTR_TYPE(rb_d, real*);
  if (!NA_IsNArray(rb_nval))
    rb_raise(rb_eArgError, "nval (11th argument) must be NArray");
  if (NA_RANK(rb_nval) != 1)
    rb_raise(rb_eArgError, "rank of nval (11th argument) must be %d", 1);
  if (NA_SHAPE0(rb_nval) != ((ijob==1||ijob==2) ? 0 : ijob==3 ? minp : 0))
    rb_raise(rb_eRuntimeError, "shape 0 of nval must be %d", (ijob==1||ijob==2) ? 0 : ijob==3 ? minp : 0);
  if (NA_TYPE(rb_nval) != NA_LINT)
    rb_nval = na_change_type(rb_nval, NA_LINT);
  nval = NA_PTR_TYPE(rb_nval, integer*);
  if (!NA_IsNArray(rb_c))
    rb_raise(rb_eArgError, "c (13th argument) must be NArray");
  if (NA_RANK(rb_c) != 1)
    rb_raise(rb_eArgError, "rank of c (13th argument) must be %d", 1);
  if (NA_SHAPE0(rb_c) != (ijob==1 ? 0 : (ijob==2||ijob==3) ? mmax : 0))
    rb_raise(rb_eRuntimeError, "shape 0 of c must be %d", ijob==1 ? 0 : (ijob==2||ijob==3) ? mmax : 0);
  if (NA_TYPE(rb_c) != NA_SFLOAT)
    rb_c = na_change_type(rb_c, NA_SFLOAT);
  c = NA_PTR_TYPE(rb_c, real*);
  {
    int shape[1];
    shape[0] = (ijob==1||ijob==2) ? 0 : ijob==3 ? minp : 0;
    rb_nval_out__ = na_make_object(NA_LINT, 1, shape, cNArray);
  }
  nval_out__ = NA_PTR_TYPE(rb_nval_out__, integer*);
  MEMCPY(nval_out__, nval, integer, NA_TOTAL(rb_nval));
  rb_nval = rb_nval_out__;
  nval = nval_out__;
  {
    int shape[2];
    shape[0] = mmax;
    shape[1] = 2;
    rb_ab_out__ = na_make_object(NA_SFLOAT, 2, shape, cNArray);
  }
  ab_out__ = NA_PTR_TYPE(rb_ab_out__, real*);
  MEMCPY(ab_out__, ab, real, NA_TOTAL(rb_ab));
  rb_ab = rb_ab_out__;
  ab = ab_out__;
  {
    int shape[1];
    shape[0] = ijob==1 ? 0 : (ijob==2||ijob==3) ? mmax : 0;
    rb_c_out__ = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  c_out__ = NA_PTR_TYPE(rb_c_out__, real*);
  MEMCPY(c_out__, c, real, NA_TOTAL(rb_c));
  rb_c = rb_c_out__;
  c = c_out__;
  {
    int shape[2];
    shape[0] = mmax;
    shape[1] = 2;
    rb_nab_out__ = na_make_object(NA_LINT, 2, shape, cNArray);
  }
  nab_out__ = NA_PTR_TYPE(rb_nab_out__, integer*);
  MEMCPY(nab_out__, nab, integer, NA_TOTAL(rb_nab));
  rb_nab = rb_nab_out__;
  nab = nab_out__;
  work = ALLOC_N(real, (mmax));
  iwork = ALLOC_N(integer, (mmax));

  slaebz_(&ijob, &nitmax, &n, &mmax, &minp, &nbmin, &abstol, &reltol, &pivmin, d, e, e2, nval, ab, c, &mout, nab, work, iwork, &info);

  free(work);
  free(iwork);
  rb_mout = INT2NUM(mout);
  rb_info = INT2NUM(info);
  return rb_ary_new3(6, rb_mout, rb_info, rb_nval, rb_ab, rb_c, rb_nab);
}

void
init_lapack_slaebz(VALUE mLapack){
  rb_define_module_function(mLapack, "slaebz", rb_slaebz, -1);
}
