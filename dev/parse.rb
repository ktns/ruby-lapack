require "pp"

RBPREFIX = "rb_"

CTYPES = {
  "INTEGER" => "integer",
  "CHARACTER" => "char",
  "REAL" => "real",
  "DOUBLE PRECISION" => "doublereal",
  "COMPLEX" => "complex",
  "COMPLEX*16" => "doublecomplex",
  "LOGICAL" => "logical"
}
NATYPES = {
  "integer" => "NA_LINT",
  "real" => "NA_SFLOAT",
  "doublereal" => "NA_DFLOAT",
  "complex" => "NA_SCOMPLEX",
  "doublecomplex" => "NA_DCOMPLEX",
  "logical" => "NA_LINT",
}


ARGS = {
  "slaqr1" => {
    "si1" => {:intent => "input", :type => "real"},
    "sr2" => {:intent => "input", :type => "real"},
    "si2" => {:intent => "input", :type => "real"},
  },
  "dlaqr1" => {
    "si1" => {:intent => "input", :type => "doublereal"},
    "sr2" => {:intent => "input", :type => "doublereal"},
    "si2" => {:intent => "input", :type => "doublereal"},
  },
  "claqr1" => {
    "s2" => {:intent => "input", :type => "complex"},
  },
  "zlaqr1" => {
    "s2" => {:intent => "input", :type => "doublecomplex"},
  },
  "slarrf" => {
    "clgapl" => {:intent => "input", :type => "real"},
    "clgapr" => {:intent => "input", :type => "real"},
    "info" => {:intent => "output", :type => "integer"},
  },
  "dlarrf" => {
    "clgapl" => {:intent => "input", :type => "doublereal"},
    "clgapr" => {:intent => "input", :type => "doublereal"},
    "info" => {:intent => "output", :type => "integer"},
  }
}

DIMS = {
  "stgex2" => {
    "work" => ["lwork"]
  },
  "dtgex2" => {
    "work" => ["lwork"]
  },
  "dlasd1" => {
    "d" => ["n"]
  },
  "slaruv" => {
    "x" => ["MAX(1,n)"]
  },
  "dlaruv" => {
    "x" => ["MAX(1,n)"]
  },
  "slasyf" => {
    "w" => ["ldw","MAX(1,nb)"]
  },
  "dlasyf" => {
    "w" => ["ldw","MAX(1,nb)"]
  },
  "clasyf" => {
    "w" => ["ldw","MAX(1,nb)"]
  },
  "zlasyf" => {
    "w" => ["ldw","MAX(1,nb)"]
  },
  "slaeda" => {
    "qptr" => ["ldqptr"]
  },
  "dlaeda" => {
    "qptr" => ["ldqptr"]
  },
  "slasdt" => {
    "inode" => ["MAX(1,n)"],
    "ndiml" => ["MAX(1,n)"],
    "ndimr" => ["MAX(1,n)"],
  },
  "dlasdt" => {
    "inode" => ["MAX(1,n)"],
    "ndiml" => ["MAX(1,n)"],
    "ndimr" => ["MAX(1,n)"],
  },
  "sgbequ" => {
    "r" => ["MAX(1,m)"]
  },
  "dgbequ" => {
    "r" => ["MAX(1,m)"]
  },
  "cgbequ" => {
    "r" => ["MAX(1,m)"]
  },
  "zgbequ" => {
    "r" => ["MAX(1,m)"]
  },
  "slaed9" => {
    "d" => ["MAX(1,n)"],
    "q" => ["ldq", "MAX(1,n)"]
  },
  "dlaed9" => {
    "d" => ["MAX(1,n)"],
    "q" => ["ldq", "MAX(1,n)"]
  },
  "slarnv" => {
    "x" => ["MAX(1,n)"],
  },
  "dlarnv" => {
    "x" => ["MAX(1,n)"],
  },
  "clarnv" => {
    "x" => ["MAX(1,n)"],
  },
  "zlarnv" => {
    "x" => ["MAX(1,n)"]
  },
  "slabrd" => {
    "d" => ["MAX(1,nb)"],
    "e" => ["MAX(1,nb)"],
    "tauq" => ["MAX(1,nb)"],
    "taup" => ["MAX(1,nb)"],
    "x" => ["ldx", "MAX(1,nb)"],
    "y" => ["ldy", "MAX(1,nb)"],
  },
  "dlabrd" => {
    "d" => ["MAX(1,nb)"],
    "e" => ["MAX(1,nb)"],
    "tauq" => ["MAX(1,nb)"],
    "taup" => ["MAX(1,nb)"],
    "x" => ["ldx", "MAX(1,nb)"],
    "y" => ["ldy", "MAX(1,nb)"],
  },
  "clabrd" => {
    "d" => ["MAX(1,nb)"],
    "e" => ["MAX(1,nb)"],
    "tauq" => ["MAX(1,nb)"],
    "taup" => ["MAX(1,nb)"],
    "x" => ["ldx", "MAX(1,nb)"],
    "y" => ["ldy", "MAX(1,nb)"],
  },
  "zlabrd" => {
    "d" => ["MAX(1,nb)"],
    "e" => ["MAX(1,nb)"],
    "tauq" => ["MAX(1,nb)"],
    "taup" => ["MAX(1,nb)"],
    "x" => ["ldx", "MAX(1,nb)"],
    "y" => ["ldy", "MAX(1,nb)"],
  },
  "slahr2" => {
    "tau" => ["MAX(1,nb)"],
    "t" => ["ldt","MAX(1,nb)"],
    "y" => ["ldy","MAX(1,nb)"],
  },
  "dlahr2" => {
    "tau" => ["MAX(1,nb)"],
    "t" => ["ldt","MAX(1,nb)"],
    "y" => ["ldy","MAX(1,nb)"],
  },
  "clahr2" => {
    "tau" => ["MAX(1,nb)"],
    "t" => ["ldt","MAX(1,nb)"],
    "y" => ["ldy","MAX(1,nb)"],
  },
  "zlahr2" => {
    "tau" => ["MAX(1,nb)"],
    "t" => ["ldt","MAX(1,nb)"],
    "y" => ["ldy","MAX(1,nb)"],
  },
  "slahrd" => {
    "tau" => ["MAX(1,nb)"],
    "t" => ["ldt","MAX(1,nb)"],
    "y" => ["ldy","MAX(1,nb)"],
  },
  "dlahrd" => {
    "tau" => ["MAX(1,nb)"],
    "t" => ["ldt","MAX(1,nb)"],
    "y" => ["ldy","MAX(1,nb)"],
  },
  "clahrd" => {
    "tau" => ["MAX(1,nb)"],
    "t" => ["ldt","MAX(1,nb)"],
    "y" => ["ldy","MAX(1,nb)"],
  },
  "zlahrd" => {
    "tau" => ["MAX(1,nb)"],
    "t" => ["ldt","MAX(1,nb)"],
    "y" => ["ldy","MAX(1,nb)"],
  },
  "slaqr2" => {
    "sr" => ["MAX(1,kbot)"],
    "si" => ["MAX(1,kbot)"],
    "v" => ["ldv", "MAX(1,nw)"],
    "t" => ["ldt", "MAX(1,nw)"],
    "wv" => ["ldwv", "MAX(1,nw)"],
  },
  "dlaqr2" => {
    "sr" => ["MAX(1,kbot)"],
    "si" => ["MAX(1,kbot)"],
    "v" => ["ldv", "MAX(1,nw)"],
    "t" => ["ldt", "MAX(1,nw)"],
    "wv" => ["ldwv", "MAX(1,nw)"],
  },
  "claqr2" => {
    "sh" => ["MAX(1,kbot)"],
    "v" => ["ldv", "MAX(1,nw)"],
    "t" => ["ldv", "MAX(1,nw)"],
    "wv" => ["ldv", "MAX(1,nw)"],
  },
  "zlaqr2" => {
    "sh" => ["MAX(1,kbot)"],
    "v" => ["ldv", "MAX(1,nw)"],
    "t" => ["ldv", "MAX(1,nw)"],
    "wv" => ["ldv", "MAX(1,nw)"],
  },
  "slaqr3" => {
    "sr" => ["MAX(1,kbot)"],
    "si" => ["MAX(1,kbot)"],
    "v" => ["ldv", "MAX(1,nw)"],
    "t" => ["ldt", "MAX(1,nw)"],
    "wv" => ["ldwv", "MAX(1,nw)"],
  },
  "dlaqr3" => {
    "sr" => ["MAX(1,kbot)"],
    "si" => ["MAX(1,kbot)"],
    "v" => ["ldv", "MAX(1,nw)"],
    "t" => ["ldt", "MAX(1,nw)"],
    "wv" => ["ldwv", "MAX(1,nw)"],
  },
  "claqr3" => {
    "sh" => ["MAX(1,kbot)"],
    "v" => ["ldv", "MAX(1,nw)"],
    "t" => ["ldv", "MAX(1,nw)"],
    "wv" => ["ldv", "MAX(1,nw)"],
  },
  "zlaqr3" => {
    "sh" => ["MAX(1,kbot)"],
    "v" => ["ldv", "MAX(1,nw)"],
    "t" => ["ldv", "MAX(1,nw)"],
    "wv" => ["ldv", "MAX(1,nw)"],
  },
  "slatrd" => {
    "w" => ["ldw", "MAX(n,nb)"]
  },
  "dlatrd" => {
    "w" => ["ldw", "MAX(n,nb)"]
  },
  "clatrd" => {
    "w" => ["ldw", "MAX(n,nb)"]
  },
  "zlatrd" => {
    "w" => ["ldw", "MAX(n,nb)"]
  },
  "clahef" => {
    "w" => ["ldw", "MAX(n,nb)"]
  },
  "zlahef" => {
    "w" => ["ldw", "MAX(n,nb)"]
  },
  "sopgtr" => {
    "ap" => ["ldap"],
    "tau" => ["ldtau"],
  },
  "dopgtr" => {
    "ap" => ["ldap"],
    "tau" => ["ldtau"],
  },
  "cupgtr" => {
    "ap" => ["ldap"],
    "tau" => ["ldtau"],
  },
  "zupgtr" => {
    "ap" => ["ldap"],
    "tau" => ["ldtau"],
  },
  "ssptrd" => {
    "ap" => ["ldap"]
  },
  "dsptrd" => {
    "ap" => ["ldap"]
  },
  "chptrd" => {
    "ap" => ["ldap"]
  },
  "zhptrd" => {
    "ap" => ["ldap"]
  },
  "sspgv" => {
    "ap" => ["ldap"]
  },
  "chpev" => {
    "ap" => ["ldap"],
  },
  "zhpev" => {
    "ap" => ["ldap"],
  },
  "chpevx" => {
    "ap" => ["ldap"],
  },
  "zhpevx" => {
    "ap" => ["ldap"],
  },
  "dspgv" => {
    "ap" => ["ldap"]
  },
  "sppequ" => {
    "ap" => ["ldap"]
  },
  "dppequ" => {
    "ap" => ["ldap"]
  },
  "cppequ" => {
    "ap" => ["ldap"]
  },
  "zppequ" => {
    "ap" => ["ldap"]
  },
  "sspgvd" => {
    "ap" => ["ldap"]
  },
  "dspgvd" => {
    "ap" => ["ldap"]
  },
  "stpcon" => {
    "ap" => ["ldap"]
  },
  "dtpcon" => {
    "ap" => ["ldap"]
  },
  "ctpcon" => {
    "ap" => ["ldap"]
  },
  "ztpcon" => {
    "ap" => ["ldap"]
  },
  "sspgvx" => {
    "ap" => ["ldap"],
  },
  "dspgvx" => {
    "ap" => ["ldap"],
  },
  "chpevd" => {
    "ap" => ["ldap"],
  },
  "zhpevd" => {
    "ap" => ["ldap"],
  },
  "sspev" => {
    "ap" => ["ldap"],
  },
  "dspev" => {
    "ap" => ["ldap"],
  },
  "sspevd" => {
    "ap" => ["ldap"],
  },
  "dspevd" => {
    "ap" => ["ldap"],
  },
  "sspevx" => {
    "ap" => ["ldap"],
  },
  "dspevx" => {
    "ap" => ["ldap"],
  },
  "sppcon" => {
    "ap" => ["ldap"]
  },
  "dppcon" => {
    "ap" => ["ldap"]
  },
  "cppcon" => {
    "ap" => ["ldap"]
  },
  "zppcon" => {
    "ap" => ["ldap"]
  },
  "chpgv" => {
    "ap" => ["ldap"]
  },
  "zhpgv" => {
    "ap" => ["ldap"]
  },
  "chpgvd" => {
    "ap" => ["ldap"]
  },
  "zhpgvd" => {
    "ap" => ["ldap"]
  },
  "chpgvx" => {
    "ap" => ["ldap"]
  },
  "zhpgvx" => {
    "ap" => ["ldap"]
  },
  "chptrf" => {
    "ap" => ["ldap"]
  },
  "zhptrf" => {
    "ap" => ["ldap"]
  },
  "ssptrf" => {
    "ap" => ["ldap"]
  },
  "dsptrf" => {
    "ap" => ["ldap"]
  },
  "csptrf" => {
    "ap" => ["ldap"]
  },
  "zsptrf" => {
    "ap" => ["ldap"]
  },
  "slaqr0" => {
    "z" => ["ldz","ihi"]
  },
  "dlaqr0" => {
    "z" => ["ldz","ihi"]
  },
  "claqr0" => {
    "z" => ["ldz","ihi"],
    "work" => ["MAX(1,lwork)"]
  },
  "zlaqr0" => {
    "z" => ["wantz ? ldz : 0","wantz ? ihi : 0"],
    "work" => ["MAX(1,lwork)"]
  },
  "slaqr4" => {
    "z" => ["ldz","ihi"]
  },
  "dlaqr4" => {
    "z" => ["ldz","ihi"]
  },
  "claqr4" => {
    "z" => ["ldz","ihi"]
  },
  "zlaqr4" => {
    "z" => ["ldz","ihi"]
  },
  "slaqr5" => {
    "z" => ["wantz ? ldz : 0","wantz ? ihiz : 0"],
    "wh" => ["ldwh", "MAX(1,nh)"]
  },
  "dlaqr5" => {
    "z" => ["wantz ? ldz : 0","wantz ? ihiz : 0"],
    "wh" => ["ldwh", "MAX(1,nh)"]
  },
  "claqr5" => {
    "z" => ["wantz ? ldz : 0","wantz ? ihiz : 0"],
    "wh" => ["ldwh", "MAX(1,nh)"]
  },
  "zlaqr5" => {
    "z" => ["wantz ? ldz : 0","wantz ? ihiz : 0"],
    "wh" => ["ldwh", "MAX(1,nh)"]
  },
  "slasd0" => {
    "u" => ["ldu", "n"]
  },
  "dlasd0" => {
    "u" => ["ldu", "n"]
  },
  "slasd4" => {
    "delta" => ["n"],
    "work" => ["n"]
  },
  "dlasd4" => {
    "delta" => ["n"],
    "work" => ["n"]
  },
  "slasdq" => {
    "e" => ['sqre==0 ? n-1 : sqre==1 ? n : 0'],
    "work" => ["4*n"]
  },
  "dlasdq" => {
    "e" => ['sqre==0 ? n-1 : sqre==1 ? n : 0'],
    "work" => ["4*n"]
  },
  "slaed3" => {
    "q2" => ["n","n"]
  },
  "dlaed3" => {
    "q2" => ["n","n"]
  },
  "slaed4" => {
    "delta" => ["n"]
  },
  "dlaed4" => {
    "delta" => ["n"]
  },
  "slaed8" => {
    "q" => ['icompq==0 ? 0 : ldq', 'icompq==0 ? 0 : n'],
    "q2" => ['icompq==0 ? 0 : ldq2', 'icompq==0 ? 0 : n']
  },
  "dlaed8" => {
    "q" => ['icompq==0 ? 0 : ldq', 'icompq==0 ? 0 : n'],
    "q2" => ['icompq==0 ? 0 : ldq2', 'icompq==0 ? 0 : n']
  },
  "claed8" => {
    "q2" => ['ldq2', 'n']
  },
  "zlaed8" => {
    "q2" => ['ldq2', 'n']
  },
  "slasq3" => {
    "z" => ['4*n0']
  },
  "dlasq3" => {
    "z" => ['4*n0']
  },
  "slasq4" => {
    "z" => ['4*n0']
  },
  "dlasq4" => {
    "z" => ['4*n0']
  },
  "slasq5" => {
    "z" => ['4*n0']
  },
  "dlasq5" => {
    "z" => ['4*n0']
  },
  "slasq6" => {
    "z" => ['4*n0']
  },
  "dlasq6" => {
    "z" => ['4*n0']
  },
  "slazq3" => {
    "z" => ['4*n0']
  },
  "dlazq3" => {
    "z" => ['4*n0']
  },
  "slazq4" => {
    "z" => ['4*n0']
  },
  "dlazq4" => {
    "z" => ['4*n0']
  },
  "slahqr" => {
    "z" => ['wantz ? ldz : 0', 'wantz ? n : 0']
  },
  "dlahqr" => {
    "z" => ['wantz ? ldz : 0', 'wantz ? n : 0']
  },
  "clahqr" => {
    "z" => ['wantz ? ldz : 0', 'wantz ? n : 0']
  },
  "zlahqr" => {
    "z" => ['wantz ? ldz : 0', 'wantz ? n : 0']
  },
  "cpbsvx" => {
    "afb" => ["ldafb","n"]
  },
  "dgesdd" => {
    "vt" => ["ldvt","n"]
  },
  "slasda" => {
    "u" => ["ldu", "MAX(1,smlsiz)"],
    "s" => ['icompq==1 ? n : icompq==0 ? 1 : 0']
  },
  "dlasda" => {
    "u" => ["ldu", "MAX(1,smlsiz)"],
    "s" => ['icompq==1 ? n : icompq==0 ? 1 : 0']
  },
  "clalsa" => {
    "rwork" => ['MAX(n,(smlsiz+1)*nrhs*3)']
  },
  "zlalsa" => {
    "rwork" => ['MAX(n,(smlsiz+1)*nrhs*3)']
  },
  "stgsen" => {
    "iwork" => ['ijob==0 ? 0 : MAX(1,liwork)']
  },
  "dtgsen" => {
    "iwork" => ['ijob==0 ? 0 : MAX(1,liwork)']
  },
  "ztgsen" => {
    "work" => ['ijob==0 ? 0 : MAX(1,lwork)'],
    "iwork" => ['ijob==0 ? 0 : MAX(1,liwork)']
  },
  "ctgsen" => {
    "work" => ['ijob==0 ? 0 : MAX(1,lwork)'],
    "iwork" => ['ijob==0 ? 0 : MAX(1,liwork)']
  },
  "slarfx" => {
    "work" => ['lsame_(&side,"L") ? n : lsame_(&side,"R") ? m : 0']
  },
  "dlarfx" => {
    "work" => ['lsame_(&side,"L") ? n : lsame_(&side,"R") ? m : 0']
  },
  "clarfx" => {
    "work" => ['lsame_(&side,"L") ? n : lsame_(&side,"R") ? m : 0']
  },
  "zlarfx" => {
    "work" => ['lsame_(&side,"L") ? n : lsame_(&side,"R") ? m : 0']
  },
  "slaebz" => {
    "nval" => ['(ijob==1||ijob==2) ? 0 : ijob==3 ? minp : 0'],
    "c" => ['ijob==1 ? 0 : (ijob==2||ijob==3) ? mmax : 0'],
    "nab" => ['mmax', '2']
  },
  "dlaebz" => {
    "nval" => ['(ijob==1||ijob==2) ? 0 : ijob==3 ? minp : 0'],
    "c" => ['ijob==1 ? 0 : (ijob==2||ijob==3) ? mmax : 0'],
    "nab" => ['mmax', '2']
  },
  "sbdsdc" => {
    "u" => ['lsame_(&compq,"I") ? ldu : 0','lsame_(&compq,"I") ? n : 0'],
    "vt" => ['lsame_(&compq,"I") ? ldvt : 0','lsame_(&compq,"I") ? n : 0'],
    "q" => ['lsame_(&compq,"I") ? ldq : 0'],
    "iq" => ['lsame_(&compq,"I") ? ldiq : 0'],
    "work" => ['MAX(1,lwork)']
  },
  "dbdsdc" => {
    "u" => ['lsame_(&compq,"I") ? ldu : 0','lsame_(&compq,"I") ? n : 0'],
    "vt" => ['lsame_(&compq,"I") ? ldvt : 0','lsame_(&compq,"I") ? n : 0'],
    "q" => ['lsame_(&compq,"I") ? ldq : 0'],
    "iq" => ['lsame_(&compq,"I") ? ldiq : 0'],
    "work" => ['MAX(1,lwork)']
  },
  "sgelsx" => {
    "work" => ['MAX((MIN(m,n))+3*n,2*(MIN(m,n))*nrhs)']
  },
  "dgelsx" => {
    "work" => ['MAX((MIN(m,n))+3*n,2*(MIN(m,n))*nrhs)']
  },
  "cgelsx" => {
    "work" => ['MIN(m,n) + MAX(n,2*(MIN(m,n))+nrhs)']
  },
  "zgelsx" => {
    "work" => ['MIN(m,n) + MAX(n,2*(MIN(m,n))+nrhs)']
  },
  "sgeevx" => {
    "iwork" => ['(lsame_(&sense,"N")||lsame_(&sense,"E")) ? 0 : 2*n-2']
  },
  "dgeevx" => {
    "iwork" => ['(lsame_(&sense,"N")||lsame_(&sense,"E")) ? 0 : 2*n-2']
  },
  "ctgex2" => {
    "q" => ['wantq ? ldq : 0', 'wantq ? n : 0'],
    "z" => ['wantq ? ldz : 0', 'wantq ? n : 0']
  },
  "ztgex2" => {
    "q" => ['wantq ? ldq : 0', 'wantq ? n : 0'],
    "z" => ['wantq ? ldz : 0', 'wantq ? n : 0']
  },
  "clalsd" => {
    "rwork" => ['9*n+2*n*smlsiz+8*n*nlvl+3*smlsiz*nrhs+(smlsiz+1)*(smlsiz+1)']
  },
  "zlalsd" => {
    "rwork" => ['9*n+2*n*smlsiz+8*n*nlvl+3*smlsiz*nrhs+(smlsiz+1)*(smlsiz+1)']
  },
  "ssbgvx" => {
    "work" => ['7*n'],
    "iwork" => ['5*n'],
  }
}


TYPES = {
  "slarrc" => {
    "vl" => "real",
    "vu" => "real",
    "d" => "real",
    "e" => "real",
    "pivmin" => "real"
  },
  "slarrb" => {
    "pivmin" => "real",
    "spdiam" => "real"
  },
  "slarre" => {
    "pivmin" => "real"
  },
  "slarrf" => {
    "pivmin" => "real",
    "spdiam" => "real"
  },
  "dlarrf" => {
    "spdiam" => "doublereal"
  },
  "slarrj" => {
    "pivmin" => "real",
    "spdiam" => "real"
  },
  "slarrv" => {
    "pivmin" => "real"
  },
  "clarrv" =>{
    "pivmin" => "real"
  },
  "zggbal" => {
    "work" => "doublereal"
  },
  "clarcm" => {
    "b" => "complex"
  },
  "zlarcm" => {
    "b" => "doublecomplex"
  },
  "clatdf" => {
    "z" => "complex",
    "rhs" => "complex",
  },
  "zlatdf" => {
    "z" => "doublecomplex",
    "rhs" => "doublecomplex",
  },
  "dggbal" => {
    "work" => "doublereal"
  },
  "zggevx" => {
    "rwork" => "doublereal"
  },
  "dlazq3" => {
    "dmin1" => "doublereal",
    "dmin2" => "doublereal",
    "dn" => "doublereal",
    "dn1" => "doublereal",
    "dn2" => "doublereal",
    "tau" => "doublereal",
  },
  "zlag2c" => {
    "a" => "doublecomplex",
    "sa" => "complex"
  },
  "clag2z" => {
    "a" => "doublecomplex",
    "sa" => "complex"
  },
  "cptts2" => {
    "b" => "complex"
  },
  "zptts2" => {
    "b" => "doublecomplex"
  },
  "cpttrs" => {
    "b" => "complex"
  },
  "zpttrs" => {
    "b" => "doublecomplex"
  },
  "zpbrfs" => {
    "ab" => "doublecomplex"
  }
}



SUBSTS = {
  "dlasd1" => {
    "d" => {"n" => "nl+nr+1"}
  },
  "sgbbrd" => {
    "d" => {"m" => "ldab"}
  },
  "dgbbrd" => {
    "d" => {"m" => "ldab"}
  },
  "cgbbrd" => {
    "d" => {"m" => "ldab"}
  },
  "zgbbrd" => {
    "d" => {"m" => "ldab"}
  },
  "slaeda" => {
    "givptr" => {"n" => "ldqptr-2"}
  },
  "dlaeda" => {
    "givptr" => {"n" => "ldqptr-2"}
  },
  "sggsvp" => {
    "u" => {"m" => "lda"}
  },
  "dggsvp" => {
    "u" => {"m" => "lda"}
  },
  "cggsvp" => {
    "u" => {"m" => "lda"}
  },
  "zggsvp" => {
    "u" => {"m" => "lda"}
  },
  "sgeequ" => {
    "r" => {"m" => "lda"}
  },
  "dgeequ" => {
    "r" => {"m" => "lda"}
  },
  "cgeequ" => {
    "r" => {"m" => "lda"}
  },
  "zgeequ" => {
    "r" => {"m" => "lda"}
  },
  "cungr2" => {
    "work" => {"m" => "lda"}
  },
  "zungr2" => {
    "work" => {"m" => "lda"}
  },
  "cungl2" => {
    "work" => {"m" => "lda"}
  },
  "zungl2" => {
    "work" => {"m" => "lda"}
  },
  "sorgr2" => {
    "work" => {"m" => "lda"}
  },
  "dorgr2" => {
    "work" => {"m" => "lda"}
  },
  "sorgl2" => {
    "work" => {"m" => "lda"}
  },
  "dorgl2" => {
    "work" => {"m" => "lda"}
  },
  "stzrqf" => {
    "tau" => {"m" => "lda"}
  },
  "dtzrqf" => {
    "tau" => {"m" => "lda"}
  },
  "ctzrqf" => {
    "tau" => {"m" => "lda"}
  },
  "ztzrqf" => {
    "tau" => {"m" => "lda"}
  },
  "stzrzf" => {
    "tau" => {"m" => "lda"}
  },
  "dtzrzf" => {
    "tau" => {"m" => "lda"}
  },
  "ctzrzf" => {
    "tau" => {"m" => "lda"}
  },
  "ztzrzf" => {
    "tau" => {"m" => "lda"}
  },
  "sgerq2" => {
    "tau" => {"m" => "lda"}
  },
  "dgerq2" => {
    "tau" => {"m" => "lda"}
  },
  "cgerq2" => {
    "tau" => {"m" => "lda"}
  },
  "zgerq2" => {
    "tau" => {"m" => "lda"}
  },
  "sgelq2" => {
    "tau" => {"m" => "lda"}
  },
  "dgelq2" => {
    "tau" => {"m" => "lda"}
  },
  "cgelq2" => {
    "tau" => {"m" => "lda"}
  },
  "zgelq2" => {
    "tau" => {"m" => "lda"}
  },
  "slatrz" => {
    "tau" => {"m" => "lda"}
  },
  "dlatrz" => {
    "tau" => {"m" => "lda"}
  },
  "clatrz" => {
    "tau" => {"m" => "lda"}
  },
  "zlatrz" => {
    "tau" => {"m" => "lda"}
  },
  "slaqps" => {
    "tau" => {"kb" => "nb"}
  },
  "dlaqps" => {
    "tau" => {"kb" => "nb"}
  },
  "claqps" => {
    "tau" => {"kb" => "nb"}
  },
  "zlaqps" => {
    "tau" => {"kb" => "nb"}
  },
  "sstevx" => {
    "z" => {"m" => "n"}
  },
  "dstevx" => {
    "z" => {"m" => "n"}
  },
  "sstein" => {
    "z" => {"m" => "n"}
  },
  "dstein" => {
    "z" => {"m" => "n"}
  },
  "cstein" => {
    "z" => {"m" => "n"}
  },
  "zstein" => {
    "z" => {"m" => "n"}
  },
  "sopgtr" => {
    "ap" => {"n" => "ldtau+1"},
  },
  "dopgtr" => {
    "ap" => {"n" => "ldtau+1"},
  },
  "cupgtr" => {
    "ap" => {"n" => "ldtau+1"},
  },
  "zupgtr" => {
    "ap" => {"n" => "ldtau+1"},
  },
  "stgsna" => {
    "s" => {"mm" => "m"}
  },
  "dtgsna" => {
    "s" => {"mm" => "m"}
  },
  "ctgsna" => {
    "s" => {"mm" => "m"}
  },
  "ztgsna" => {
    "s" => {"mm" => "m"}
  },
  "strsna" => {
    "s" => {"mm" => "m"}
  },
  "dtrsna" => {
    "s" => {"mm" => "m"}
  },
  "ctrsna" => {
    "s" => {"mm" => "m"}
  },
  "ztrsna" => {
    "s" => {"mm" => "m"}
  },
  "sggsvd" => {
    "u" => {"m" => "lda"},
    "v" => {"p" => "ldb"},
  },
  "dggsvd" => {
    "u" => {"m" => "lda"},
    "v" => {"p" => "ldb"},
  },
  "cggsvd" => {
    "u" => {"m" => "lda"},
    "v" => {"p" => "ldb"},
  },
  "zggsvd" => {
    "u" => {"m" => "lda"},
    "v" => {"p" => "ldb"},
  },
  "dlansp" => {
    "work" => {"lwork" => '(lsame_(&norm,"I") || lsame_(&norm,"1") || lsame_(&norm,"0")) ? n : 0'}
  },
  "ssptrd" => {
    "d" => {"n" => '(int)(sqrt((double)8*ldap+1)-1)/2'},
  },
  "dsptrd" => {
    "d" => {"n" => '(int)(sqrt((double)8*ldap+1)-1)/2'},
  },
  "chptrd" => {
    "d" => {"n" => '(int)(sqrt((double)8*ldap+1)-1)/2'},
  },
  "zhptrd" => {
    "d" => {"n" => '(int)(sqrt((double)8*ldap+1)-1)/2'},
  },
  "sppequ" => {
    "s" => {"n" => '(int)(sqrt((double)8*ldap+1)-1)/2'},
  },
  "dppequ" => {
    "s" => {"n" => '(int)(sqrt((double)8*ldap+1)-1)/2'},
  },
  "cppequ" => {
    "s" => {"n" => '(int)(sqrt((double)8*ldap+1)-1)/2'},
  },
  "zppequ" => {
    "s" => {"n" => '(int)(sqrt((double)8*ldap+1)-1)/2'},
  },
  "chpevd" => {
    "w" => {"n" => '(int)(sqrt((double)8*ldap+1)-1)/2'},
  },
  "zhpevd" => {
    "w" => {"n" => '(int)(sqrt((double)8*ldap+1)-1)/2'},
  },
  "sspev" => {
    "w" => {"n" => '(int)(sqrt((double)8*ldap+1)-1)/2'},
  },
  "dspev" => {
    "w" => {"n" => '(int)(sqrt((double)8*ldap+1)-1)/2'},
  },
  "sspevd" => {
    "w" => {"n" => '(int)(sqrt((double)8*ldap+1)-1)/2'},
  },
  "dspevd" => {
    "w" => {"n" => '(int)(sqrt((double)8*ldap+1)-1)/2'},
  },
  "chpev" => {
    "w" => {"n" => '(int)(sqrt((double)8*ldap+1)-1)/2'}
  },
  "zhpev" => {
    "w" => {"n" => '(int)(sqrt((double)8*ldap+1)-1)/2'}
  },
  "sspgv" => {
    "bp" => {"n" => '(int)(sqrt((double)8*ldap+1)-1)/2'},
  },
  "chpgvx" => {
    "bp" => {"n" => '(int)(sqrt((double)8*ldap+1)-1)/2'},
  },
  "zhpgvx" => {
    "bp" => {"n" => '(int)(sqrt((double)8*ldap+1)-1)/2'},
  },
  "dspgv" => {
    "bp" => {"n" => '(int)(sqrt((double)8*ldap+1)-1)/2'},
  },
  "sspgvd" => {
    "bp" => {"n" => '(int)(sqrt((double)8*ldap+1)-1)/2'},
  },
  "dspgvd" => {
    "bp" => {"n" => '(int)(sqrt((double)8*ldap+1)-1)/2'},
  },
  "chpgv" => {
    "bp" => {"n" => '(int)(sqrt((double)8*ldap+1)-1)/2'},
  },
  "zhpgv" => {
    "bp" => {"n" => '(int)(sqrt((double)8*ldap+1)-1)/2'},
  },
  "chpgvd" => {
    "bp" => {"n" => '(int)(sqrt((double)8*ldap+1)-1)/2'},
  },
  "zhpgvd" => {
    "bp" => {"n" => '(int)(sqrt((double)8*ldap+1)-1)/2'},
  },
  "chptrf" => {
    "ipiv" => {"n" => '(int)(sqrt((double)8*ldap+1)-1)/2'}
  },
  "zhptrf" => {
    "ipiv" => {"n" => '(int)(sqrt((double)8*ldap+1)-1)/2'}
  },
  "ssptrf" => {
    "ipiv" => {"n" => '(int)(sqrt((double)8*ldap+1)-1)/2'}
  },
  "dsptrf" => {
    "ipiv" => {"n" => '(int)(sqrt((double)8*ldap+1)-1)/2'}
  },
  "csptrf" => {
    "ipiv" => {"n" => '(int)(sqrt((double)8*ldap+1)-1)/2'}
  },
  "zsptrf" => {
    "ipiv" => {"n" => '(int)(sqrt((double)8*ldap+1)-1)/2'}
  },
  "stpcon" => {
    "work" => {"n" => '(int)(sqrt((double)8*ldap+1)-1)/2'},
  },
  "dtpcon" => {
    "work" => {"n" => '(int)(sqrt((double)8*ldap+1)-1)/2'},
  },
  "ctpcon" => {
    "work" => {"n" => '(int)(sqrt((double)8*ldap+1)-1)/2'},
  },
  "ztpcon" => {
    "work" => {"n" => '(int)(sqrt((double)8*ldap+1)-1)/2'},
  },
  "sppcon" => {
    "work" => {"n" => '(int)(sqrt((double)8*ldap+1)-1)/2'}
  },
  "dppcon" => {
    "work" => {"n" => '(int)(sqrt((double)8*ldap+1)-1)/2'}
  },
  "cppcon" => {
    "work" => {"n" => '(int)(sqrt((double)8*ldap+1)-1)/2'}
  },
  "zppcon" => {
    "work" => {"n" => '(int)(sqrt((double)8*ldap+1)-1)/2'}
  },
  "sspgvx" => {
    "bp" => {"n" => '(int)(sqrt((double)8*ldap+1)-1)/2'},
    "z" => {
      "m" => 'lsame_(&range,"A") ? n : lsame_(&range,"I") ? iu-il+1 : 0'
    }
  },
  "dspgvx" => {
    "bp" => {"n" => '(int)(sqrt((double)8*ldap+1)-1)/2'},
    "z" => {
      "m" => 'lsame_(&range,"A") ? n : lsame_(&range,"I") ? iu-il+1 : 0'
    }
  },
  "sspevx" => {
    "w" => {"n" => '(int)(sqrt((double)8*ldap+1)-1)/2'},
    "z" => {
      "m" => 'lsame_(&range,"A") ? n : lsame_(&range,"I") ? iu-il+1 : 0'
    }
  },
  "dspevx" => {
    "w" => {"n" => '(int)(sqrt((double)8*ldap+1)-1)/2'},
    "z" => {
      "m" => 'lsame_(&range,"A") ? n : lsame_(&range,"I") ? iu-il+1 : 0'
    }
  },
  "chpevx" => {
    "w" => {"n" => '(int)(sqrt((double)8*ldap+1)-1)/2'},
    "z" => {
      "m" => 'lsame_(&range,"A") ? n : lsame_(&range,"I") ? iu-il+1 : 0'
    }
  },
  "zhpevx" => {
    "w" => {"n" => '(int)(sqrt((double)8*ldap+1)-1)/2'},
    "z" => {
      "m" => 'lsame_(&range,"A") ? n : lsame_(&range,"I") ? iu-il+1 : 0'
    }
  },
  "sstevr" => {
    "z" => {"m" => 'lsame_(&range,"I") ? iu-il+1 : n'}
  },
  "dstevr" => {
    "z" => {"m" => 'lsame_(&range,"I") ? iu-il+1 : n'}
  },
  "ssyevr" => {
    "z" => {"m" => 'lsame_(&range,"I") ? iu-il+1 : n'}
  },
  "dsyevr" => {
    "z" => {"m" => 'lsame_(&range,"I") ? iu-il+1 : n'}
  },
  "dsyevx" => {
    "z" => {"m" => 'lsame_(&range,"I") ? iu-il+1 : n'}
  },
  "ssyevx" => {
    "z" => {"m" => 'lsame_(&range,"I") ? iu-il+1 : n'}
  },
  "ssygvx" => {
    "z" => {
      "m" => 'lsame_(&range,"A") ? n : lsame_(&range,"I") ? iu-il+1 : 0'
    }
  },
  "dsygvx" => {
    "z" => {
      "m" => 'lsame_(&range,"A") ? n : lsame_(&range,"I") ? iu-il+1 : 0'
    }
  },
  "cheevr" => {
    "z" => {
      "m" => 'lsame_(&range,"A") ? n : lsame_(&range,"I") ? iu-il+1 : 0'
    }
  },
  "zheevr" => {
    "z" => {
      "m" => 'lsame_(&range,"A") ? n : lsame_(&range,"I") ? iu-il+1 : 0'
    }
  },
  "chbevx" => {
    "z" => {
      "m" => 'lsame_(&range,"A") ? n : lsame_(&range,"I") ? iu-il+1 : 0'
    }
  },
  "zhbevx" => {
    "z" => {
      "m" => 'lsame_(&range,"A") ? n : lsame_(&range,"I") ? iu-il+1 : 0'
    }
  },
  "sstemr" => {
    "z" => {
      "m" => 'lsame_(&range,"A") ? n : lsame_(&range,"I") ? iu-il+1 : 0'
    }
  },
  "dstemr" => {
    "z" => {
      "m" => 'lsame_(&range,"A") ? n : lsame_(&range,"I") ? iu-il+1 : 0'
    }
  },
  "cstemr" => {
    "z" => {
      "m" => 'lsame_(&range,"A") ? n : lsame_(&range,"I") ? iu-il+1 : 0'
    }
  },
  "zstemr" => {
    "z" => {
      "m" => 'lsame_(&range,"A") ? n : lsame_(&range,"I") ? iu-il+1 : 0'
    }
  },
  "sstegr" => {
    "z" => {
      "m" => 'lsame_(&range,"A") ? n : lsame_(&range,"I") ? iu-il+1 : 0'
    }
  },
  "dstegr" => {
    "z" => {
      "m" => 'lsame_(&range,"A") ? n : lsame_(&range,"I") ? iu-il+1 : 0'
    }
  },
  "cstegr" => {
    "z" => {
      "m" => 'lsame_(&range,"A") ? n : lsame_(&range,"I") ? iu-il+1 : 0'
    }
  },
  "zstegr" => {
    "z" => {
      "m" => 'lsame_(&range,"A") ? n : lsame_(&range,"I") ? iu-il+1 : 0'
    }
  },
  "ssbevx" => {
    "z" => {
      "m" => 'lsame_(&range,"A") ? n : lsame_(&range,"I") ? iu-il+1 : 0'
    }
  },
  "dsbevx" => {
    "z" => {
      "m" => 'lsame_(&range,"A") ? n : lsame_(&range,"I") ? iu-il+1 : 0'
    }
  },
  "cheevx" => {
    "z" => {
      "m" => 'lsame_(&range,"A") ? n : lsame_(&range,"I") ? iu-il+1 : 0'
    }
  },
  "zheevx" => {
    "z" => {
      "m" => 'lsame_(&range,"A") ? n : lsame_(&range,"I") ? iu-il+1 : 0'
    }
  },
  "chegvx" => {
    "w" => {
      "m" => 'lsame_(&range,"A") ? n : lsame_(&range,"I") ? iu-il+1 : 0'
    }
  },
  "zhegvx" => {
    "w" => {
      "m" => 'lsame_(&range,"A") ? n : lsame_(&range,"I") ? iu-il+1 : 0'
    }
  },
  "ssbgvx" => {
    "ifail" => {
      "m" => 'lsame_(&range,"A") ? n : lsame_(&range,"I") ? iu-il+1 : 0'
    }
  },
  "dsbgvx" => {
    "ifail" => {
      "m" => 'lsame_(&range,"A") ? n : lsame_(&range,"I") ? iu-il+1 : 0'
    }
  },
  "stgex2" => {
    "work" => {"lwork" => 'MAX(1,(MAX(n*(n2+n1),(n2+n1)*(n2+n1)*2)))'}
  },
  "dtgex2" => {
    "work" => {"lwork" => 'MAX(1,(MAX(n*(n2+n1),(n2+n1)*(n2+n1)*2)))'}
  },
  "zlange" => {
    "work" => {"lwork" => 'lsame_(&norm,"I") ? m : 0'}
  },
  "zlanhs" => {
    "work" => {"lwork" => 'lsame_(&norm,"I") ? n : 0'}
  },
  "spprfs" => {
    "ap" => {"n" => "ldb"}
  },
  "dpprfs" => {
    "ap" => {"n" => "ldb"}
  },
  "cpprfs" => {
    "ap" => {"n" => "ldb"}
  },
  "zpprfs" => {
    "ap" => {"n" => "ldb"}
  },
  "chpsv" => {
    "ap" => {"n" => "ldb"}
  },
  "zhpsv" => {
    "ap" => {"n" => "ldb"}
  },
  "sspsv" => {
    "ap" => {"n" => 'ldb'}
  },
  "dspsv" => {
    "ap" => {"n" => 'ldb'}
  },
  "cspsv" => {
    "ap" => {"n" => 'ldb'}
  },
  "zspsv" => {
    "ap" => {"n" => 'ldb'}
  },
  "stprfs" => {
    "ap" => {"n" => 'ldb'}
  },
  "dtprfs" => {
    "ap" => {"n" => 'ldb'}
  },
  "ctprfs" => {
    "ap" => {"n" => 'ldb'}
  },
  "ztprfs" => {
    "ap" => {"n" => 'ldb'}
  },
  "slasd0" => {
    "e" => {
      "m" => "sqre == 0 ? n : sqre == 1 ? n+1 : 0",
      "ldu" => "n"
    },
    "vt" => {"ldvt" => "m"}
  },
  "dlasd0" => {
    "e" => {"m" => "sqre == 0 ? n : sqre == 1 ? n+1 : 0", "ldu" => "n"},
    "u" => {"ldu" => "n"},
    "vt" => {"ldvt" => "n"}
  },
  "slasd3" => {
    "u2" => {"ldu2" => "n", "n" => "nl + nr + 1"},
    "vt2" => {"ldvt2" => "n"},
    "vt" => {"m" => "n+sqre"},
    :order => {
      "u2" => ["n", "ldu2"]
    }
  },
  "dlasd3" => {
    "u2" => {"ldu2" => "n", "n" => "nl + nr + 1"},
    "vt2" => {"ldvt2" => "n"},
    "vt" => {"m" => "n + sqre"},
    :order => {
      "u2" => ["n", "ldu2"]
    }
  },
  "slasda" => {
    "e" => {"m" => "sqre == 0 ? n : sqre == 1 ? n+1 : 0", "ldu" => "n"}
  },
  "dlasda" => {
    "e" => {"m" => "sqre == 0 ? n : sqre == 1 ? n+1 : 0", "ldu" => "n"}
  },
  "slals0" => {
    "bx" => {"ldbx" => "n"},
  },
  "dlals0" => {
    "bx" => {"ldbx" => "n"},
  },
  "clals0" => {
    "bx" => {"ldbx" => "n"}
  },
  "zlals0" => {
    "bx" => {"ldbx" => "n"}
  },
  "slalsa" => {
    "bx" => {"ldbx" => "n"},
  },
  "dlalsa" => {
    "bx" => {"ldbx" => "n"},
  },
  "clalsa" => {
    "bx" => {"ldbx" => "n"},
  },
  "zlalsa" => {
    "bx" => {"ldbx" => "n"},
  },
  "slasd2" => {
    "u2" => {"ldu2" => "n"},
    "vt2" => {"ldvt2" => "m"}
  },
  "dlasd2" => {
    "u2" => {"ldu2" => "n"},
    "vt2" => {"ldvt2" => "m"}
  },
  "slasd6" => {
    "vf" => {"m" => "n + sqre", "n" => "nl + nr + 1"},
    :order => {"vf" => ["n", "m"]}
  },
  "dlasd6" => {
    "vf" => {"m" => "n + sqre", "n" => "nl + nr + 1"},
    :order => {"vf" => ["n", "m"]}
  },
  "dgesdd" => {
    "vt" => {"ldvt" => '(lsame_(&jobz,"A")||(lsame_(&jobz,"O")&&(m>=n))) ? n : lsame_(&jobz,"S") ? MIN(m,n) : 0'}
  },
  "zlaqr4" => {
    "z" => {"ldz" => 'wantz ? MAX(1,ihiz) : 1'}
  },
  "slaqr5" => {
    "z" => {"ldz" => 'n'}
  },
  "dlaqr5" => {
    "z" => {"ldz" => 'n'}
  },
  "sgelsd" => {
    "iwork" => {
      "c__0" => "0",
      "c__9" => "9",
      "smlsiz" => 'ilaenv_(&c__9,"DGELSD"," ",&c__0,&c__0,&c__0,&c__0,(ftnlen)6,(ftnlen)1)',
      "nlvl" => 'MAX(0,((int)(log(((double)(MIN(m,n)))/(smlsiz+1))/log(2.0))+1))',
      "liwork" => '3*(MIN(m,n))*nlvl+11*(MIN(m,n))',
    },
    :order => {
      "iwork" => ["c__0", "c__9", "smlsiz", "nlvl", "liwork"]
    }
  },
  "dgelsd" => {
    "iwork" => {
      "c__0" => "0",
      "c__9" => "9",
      "smlsiz" => 'ilaenv_(&c__9,"DGELSD"," ",&c__0,&c__0,&c__0,&c__0,(ftnlen)6,(ftnlen)1)',
      "nlvl" => 'MAX(0,((int)(log(((double)(MIN(m,n)))/(smlsiz+1))/log(2.0))+1))',
      "liwork" => '3*(MIN(m,n))*nlvl+11*(MIN(m,n))',
    },
    :order => {
      "iwork" => ["c__0", "c__9", "smlsiz", "nlvl", "liwork"]
    }
  },
  "cgelsd" => {
    "rwork" => {
      "c__0" => "0",
      "c__9" => "9",
      "smlsiz" => 'ilaenv_(&c__9,"CGELSD"," ",&c__0,&c__0,&c__0,&c__0,(ftnlen)6,(ftnlen)1)',
      "nlvl" => 'MAX(0,(int)(log(1.0*MIN(m,n)/(smlsiz+1))/log(2.0)))',
      "lrwork" => 'm>=n ? 10*n+2*n*smlsiz+8*n*nlvl+3*smlsiz*nrhs+(smlsiz+1)*(smlsiz+1) : 10*m+2*m*smlsiz+8*m*nlvl+2*smlsiz*nrhs+(smlsiz+1)*(smlsiz+1)',
    },
    "iwork" => {
      "liwork" => 'MAX(1,3*(MIN(m,n))*nlvl+11*(MIN(m,n)))'
    },
    :order => {
      "rwork" => ["c__0","c__9","smlsiz","nlvl","lrwork"]
    }
  },
  "zgelsd" => {
    "rwork" => {
      "c__9" => "9",
      "c__0" => "0",
      "smlsiz" => 'ilaenv_(&c__9,"ZGELSD"," ",&c__0,&c__0,&c__0,&c__0,(ftnlen)6,(ftnlen)1)',
      "nlvl" => 'MAX(0,(int)(log(1.0*MIN(m,n)/(smlsiz+1))/log(2.0)))',
      "lrwork" => 'm>=n ? 10*n+2*n*smlsiz+8*n*nlvl+3*smlsiz*nrhs+(smlsiz+1)*(smlsiz+1) : 10*m+2*m*smlsiz+8*m*nlvl+2*smlsiz*nrhs+(smlsiz+1)*(smlsiz+1)',
    },
    "iwork" => {
      "liwork" => 'MAX(1,3*(MIN(m,n))*nlvl+11*(MIN(m,n)))'
    },
    :order => {
      "rwork" => ["c__0","c__9","smlsiz","nlvl","lrwork"]
    }
  },
  "sormbr" => {
    "a" => {"nq" => 'lsame_(&side,"L") ? m : lsame_(&side,"R") ? n : 0'}
  },
  "dormbr" => {
    "a" => {"nq" => 'lsame_(&side,"L") ? m : lsame_(&side,"R") ? n : 0'}
  },
  "cunmbr" => {
    "a" => {"nq" => 'lsame_(&side,"L") ? m : lsame_(&side,"R") ? n : 0'}
  },
  "zunmbr" => {
    "a" => {"nq" => 'lsame_(&side,"L") ? m : lsame_(&side,"R") ? n : 0'}
  },
  "sbdsdc" => {
    "q" => {
      "c__0" => "0",
      "c__9" => "9",
      "smlsiz" => 'ilaenv_(&c__9, "SBDSDC", " ", &c__0, &c__0, &c__0, &c__0, (ftnlen)6, (ftnlen)1)',
      "ldq" => 'lsame_(&compq,"P") ? n*(11+2*smlsiz+8*(int)(log(((double)n)/(smlsiz+1))/log(2.0))) : 0'
    },
    "iq" => {
      "ldiq" => 'lsame_(&compq,"P") ? n*(3+3*(int)(log(((double)n)/(smlsiz+1))/log(2.0))) : 0'
    },
    "work" => {"lwork" => 'lsame_(&compq,"N") ? 4*n : lsame_(&compq,"P") ? 6*n : lsame_(&compq,"I") ? 3*n*n+4*n : 0'},
    :order => {
      "q" => ["c__0","c__9","smlsiz","ldq"]
    }
  },
  "dbdsdc" => {
    "q" => {
      "c__0" => "0",
      "c__9" => "9",
      "smlsiz" => 'ilaenv_(&c__9, "DBDSDC", " ", &c__0, &c__0, &c__0, &c__0, (ftnlen)6, (ftnlen)1)',
      "ldq" => 'lsame_(&compq,"P") ? n*(11+2*smlsiz+8*(int)(log(((double)n)/(smlsiz+1))/log(2.0))) : 0'
    },
    "iq" => {
      "ldiq" => 'lsame_(&compq,"P") ? n*(3+3*(int)(log(((double)n)/(smlsiz+1))/log(2.0))) : 0'
    },
    "work" => {"lwork" => 'lsame_(&compq,"N") ? 4*n : lsame_(&compq,"P") ? 6*n : lsame_(&compq,"I") ? 3*n*n+4*n : 0'},
    :order => {
      "q" => ["c__0","c__9","smlsiz","ldq"]
    }
  },
  "clalsd" => {
    "rwork" => {
      "nlvl" => '( (int)( log(((double)n)/(smlsiz+1))/log(2.0) ) ) + 1'
    },
  },
  "zlalsd" => {
    "rwork" => {
      "nlvl" => '( (int)( log(((double)n)/(smlsiz+1))/log(2.0) ) ) + 1'
    },
  },
  "claed8" => {
    "q2" => {"ldq2" => "n"}
  },
  "zlaed8" => {
    "q2" => {"ldq2" => "n"}
  }
}



#IGNORE = %w(sgees dgees cgees zgees sgeesx dgeesx cgeesx zgeesx sgges dgges cgges zgges sggesx dggesx cggesx zggesx)


def get_cobj(name, type)
  case type
  when "integer"
    return "  #{name} = NUM2INT(#{RBPREFIX}#{name});\n"
    dimdefs[name] = true
  when "real"
    return "  #{name} = (real)NUM2DBL(#{RBPREFIX}#{name});\n"
  when "doublereal"
    return "  #{name} = NUM2DBL(#{RBPREFIX}#{name});\n"
  when "complex"
    code =<<"EOF"
  #{name}.r = (real)NUM2DBL(rb_funcall(#{RBPREFIX}#{name}, rb_intern("real"), 0));
  #{name}.i = (real)NUM2DBL(rb_funcall(#{RBPREFIX}#{name}, rb_intern("imag"), 0));
EOF
    return code
  when "doublecomplex"
    code =<<"EOF"
  #{name}.r = NUM2DBL(rb_funcall(#{RBPREFIX}#{name}, rb_intern("real"), 0));
  #{name}.i = NUM2DBL(rb_funcall(#{RBPREFIX}#{name}, rb_intern("imag"), 0));
EOF
    return code
  when "char"
    return "  #{name} = StringValueCStr(#{RBPREFIX}#{name})[0];\n"
  when "logical"
    return "  #{name} = (#{RBPREFIX}#{name} == Qtrue);\n"
  else
    raise "type (#{type}) is not defined in #{name}"
  end
end

def get_robj(name, type, flag=false)
  case type
  when "integer"
    cname =  flag ? "(*#{name})" : name
    return "  #{RBPREFIX}#{name} = INT2NUM(#{cname});\n"
  when "real", "doublereal"
    cname =  flag ? "(*#{name})" : name
    return "  #{RBPREFIX}#{name} = rb_float_new((double)#{cname});\n"
  when "complex", "doublecomplex"
    if flag
      r = "(#{name}->r)"
      i = "(#{name}->i)"
    else
      r = "(#{name}.r)"
      i = "(#{name}.i)"
    end
    return "  #{RBPREFIX}#{name} = rb_funcall(rb_gv_get(\"Complex\"), rb_intern(\"new\"), 2, rb_float_new((double)#{r}), rb_float_new((double)#{i}));\n"
  when "char"
    return "  #{RBPREFIX}#{name} = rb_str_new(&#{name},1);\n"
  when "logical"
    return "  #{RBPREFIX}#{name} = #{name} ? Qtrue : Qfalse;\n"
  else
    raise "type (#{type}) is not defined in #{name}"
  end
end

def pow(str)
  str = str.gsub(/([A-Z\d]+?)\*\*([A-Z\d]+)/, 'pow(\\1,\\2)')
  str = str.gsub(/\(([^\)]+)\)\*\*([A-Z\d]+)/, 'pow(\\1,\\2)')
  str.gsub!(/lg (\w+)/, 'LG(\\1)')
  str.gsub!(/(\w) (LG\()/, '\\1*\\2')
  return str
end


def get_vname(str)
  v = pow(str).strip.downcase.gsub(/lg\(/,"LG(").gsub(/([^a-z])max/,'\\1MAX').gsub(/^max/,'MAX').gsub(/min/,"MIN")
  if /if / =~ v || (/\,/ =~ v && /(MIN|MAX|pow)\(/ !~ v) || / the / =~ v
    raise "vname is invalid #{v}"
  end
  v
end


def get_dims(str)
  if /^\((.*)\)\.?$/ =~ str
    str = $1.strip
  end
  if /(max|min)/i =~ str
    dims = Array.new
    while (str && str != "")
      if /^((((MAX|MIN)\s*\([^\,]+\,)|[^\,])+)\,?(.*)$/i =~ str
        dns = $1.strip
        str = $5.strip
        dims.push dns
        next
      else
        /^([^\,]+)(\,.*)$/ =~ str
        dims.push $1
        str = $2.strip
      end
      if /^\,(.*)/ =~ str
        str = $1.strip
      end
    end
  else
    dims = str.split(",").collect{|dim| dim.sub(/\.$/,"")}
  end
  dims.collect{|dim| get_vname(dim)}
end

AO = {"and"=>"&&","or"=>"||"}
def get_cond(cond,v=nil)
  cond.sub!(/\.$/,"")
  if /^(.+?)\s*(and|or)\s*(.+)$/ =~ cond
    c0 = $1.strip
    ao = $2
    c1 = $3.strip
    if /^([A-Z\d]+)\s*=/ =~ c0
      v = get_vname($1)
    else
      v = nil
    end
    cond = "((#{get_cond(c0,v)}) #{AO[ao]} (#{get_cond(c1,v)}))"
  elsif /^(.+?)\s*=\s*\'(.+)'$/ =~ cond
    cond = "lsame_(&#{$1.downcase},\"#{$2}\")"
  else
    if /=.*=/ =~ cond
      conds = cond.split("=").collect{|c| c.strip.downcase}
      cond = Array.new
      (conds.length-1).times{|i|
        cond.push "(#{conds[i]}==#{conds[i+1]})"
      }
      cond = cond.join("&&")
    elsif v && (/^\'(.+)\'$/ =~ cond)
      cond = "lsame_(&#{v},\"#{$1}\")"
    elsif v && (/^\d+$/ =~ cond)
      cond = "#{v} == #{cond}"
    else
      cond = cond.gsub(/=/,"==").downcase
    end
  end
  cond
end


def get_vars(dim)
  ary = Array.new
  dim.gsub(/MAX\(/,",").gsub(/MIN\(/,",").gsub(/log\(/,",").gsub(/abs\(/,",").gsub(/sqrt\(/,",").gsub(/pow\(/,",").gsub(/LG\(/,",").gsub(/lsame_\(\&[^,]+/,",").gsub(/[\(\)\+\-\*\/:\?=\&\|]+/,",").split(",").each{|d|
    d.strip!
    next if (d == "") || (/^\d+$/ =~ d) || /^\"[^\"]+\"$/ =~ d
    ary.push d
  }
  ary
end




def parse_file(fname)
  flag_sub = false
  subr = nil
  sub_type = nil
  flag_pur = false
  purpose = nil
  flag_arg = false
  args = nil
  flag_fd = false
  fd = nil
  File.foreach(fname){|line|

    if /This routine is not for general use\./ =~ line
      return nil
    end
    if /^\*\s+\.\. Parameters \.\./ =~ line || /^\*\s+\.\. Executable Statements \.\./ =~ line
      break
    end

    if flag_sub
      if (/^     \$\s* (.+)$/ =~ line) || (/^     \+\s* (.+)$/ =~ line)
        subr << " " << $1.chomp
      else
        flag_sub = false
      end
      next
    elsif /^      SUBROUTINE/ =~ line
      subr = line.chomp
      flag_sub = true
      next
    elsif /^      ([A-Z\s\d\*]+)\s+FUNCTION/ =~ line
      subr = line.chomp
      flag_sub = true
    elsif /^      FUNCTION/ =~ line && /^[sd]laneg/ =~ File.basename(fname)
      subr = line.chomp
      flag_sub = true
    end

    if flag_pur
      if /^\*\s+Arguments$/ =~ line
        flag_pur = false
        args = line
        flag_arg = true
        next
      else
        case File.basename(fname)
        when /^[sdcz]laqr1/
          if /^\*       N      \(input\) integer$/ =~ line
            flag_pur = false
            args = line
            flag_arg = true
            next
          end
        when /^[sdcz]laqr2/, /^[sdcz]laqr3/
          if /^\*     WANTT   \(input\) LOGICAL$/ =~ line
            flag_pur = false
            args = line
            flag_arg = true
            next
          end
        when /^[sdcz]laqr5/
          if /^\*      WANTT  \(input\) logical scalar$/ =~ line
            flag_pur = false
            args = line
            flag_arg = true
            next
          end
        when /^[sd]lasq4/, /^[sd]lazq4/
          if /^\*  I0    \(input\) INTEGER$/ =~ line
            flag_pur = false
            args = line
            flag_arg = true
            next
          end
        end
        purpose << line
        next
      end
    elsif /^\*\s+Purpose$/ =~ line
      purpose = line
      flag_pur = true
      next
    else
      case File.basename(fname)
      when /^[sdcz]laqr1/
        if /^\*       Given a 2-by-2 or 3-by-3 matrix H, [SDCZ]LAQR1 sets v to a$/ =~ line
          purpose = line
          flag_pur = true
        end
      when /^[sdcz]laqr2/
        if /^\*     This subroutine is identical to [SDCZ]LAQR3 except that it avoids$/ =~ line
          purpose = line
          flag_pur = true
        end
      when /^[sdcz]laqr3/
        if /^\*     Aggressive early deflation:$/ =~ line
          purpose = line
          flag_pur = true
        end
      when /^[sdcz]laqr5/
        if /^\*     This auxiliary subroutine called by [SDCZ]LAQR0 performs a$/ =~ line
          purpose = line
          flag_pur = true
        end
      end
    end

    if flag_arg
      if /^\*\s+Further Details$/ =~ line || /^\*\s+={40,}$/ =~ line
        flag_arg = false
        fd = line
        flag_fd = true
        next
      else
        args << line
        next
      end
    end
    
    if flag_fd
      fd << line
      next
    end
  }

  unless subr
    raise "subr not found #{fname}"
  end
  unless purpose
    raise "purpose not found #{fname}"
  end
  unless args
    raise "args not found #{fname}"
  end

  help = subr + "\n\n" + purpose + "\n" + args
  help << "\n" + fd if fd

  if /^      SUBROUTINE\s+([A-Z\d]+)\(([^\)]+)\)/ =~ subr
    sub_name = $1.downcase
    arg_names = $2
    sub_type = :sub
  elsif /^      ([A-Z\s\*\d]+[A-Z\d])\s+FUNCTION\s+([A-Z\d]+)\(([^\)]+)\)/ =~ subr
    f_type = $1.strip
    sub_name = $2.downcase
    arg_names = $3
    sub_type = :func
    func_type = CTYPES[f_type]
    unless func_type
      raise "func_type #{f_type} is not defined"
    end
  elsif /^      FUNCTION\s+([A-Z\d]+)\(([^\)]+)\)/ =~ subr
    sub_name = $1.downcase
    arg_names = $2
    sub_type = :func
    case File.basename(fname)
    when /^[sd]laneg/
      func_type = "integer"
    else
      raise "function name is invalid #{subr}"
    end
  else
    raise "subroutine or function name is invalid #{subr}"
  end
  arg_names = arg_names.split(",").collect{|arg| arg.strip.downcase}
  unless arg_names
    raise "arg_names is nil #{fname}"
  end

  if @@debug
    p sub_name
    p arg_names
  end

  flag = false
  flag1 = false
  ary = Array.new
  args.each{|line|
    case line
    when /^\*\s*$/, /^\*  Arguments$/, /^\*  =========$/
      flag = false
      flag1 = false
      next
    end

    if /^\*\s+([A-Z_\d]+)\s+\((input|output|workspace|in|external procedure)[^\)]*\)\s+([A-Za-z]+)/ =~ line
      name = $1
      intent = $2
      type = $3
      ary.push line.chomp
      if /array/i =~ line
        flag = true
      elsif (/^LD/=~name || /L.*WORK/=~name) && /input/=~intent && /INTEGER/i=~type
        flag1 = true
      end
    elsif flag
      line.sub!(/^\*\s*/,"")
      if /[Ii]f / =~ line || /where / =~ line || /^[a-z]/ =~ line || /[a-z]$/ =~ ary[-1] || /is not referenced/ =~ line || /dimension (of [A-Z\d]+ must be )?at least/ =~ ary[-1] || /^Otherwise, / =~ line
        if /^[^\s]+ is the/ =~ line || /^If [^\.]+ must contain/ =~ line || /where in / =~ line || /At each / =~ line || /^[^\.]+\, [A-Z\d]+ contains / =~ line || /[A-Z]+\(\d\) and [A-Z]+\(\d\) contain the/ =~ line
          flag = false
          next
        end
        if  /^\( lg\(\s*\w?\s*\) = smallest integer \w$/ =~ line || /^such that 2\^\w \>= \w \)$/ =~ line
          next
        end
        while /^(.*)On exit/ =~ line || /^(.*?)[^\.]+ on exit/ =~ line || /^(.*)On entry/ =~ line || /^(.*)The/ =~ line || /^(.*?)[^\.]+ the [^d]/=~ line || /^(.*?)[^\.]+ is an /=~ line
          line = $1.strip
          flag = false
        end
        if line != ""
          if /^[A-Z]/ =~ line && /[\.\,;]$/ !~ ary[-1] && /if$/ !~ ary[-1]
            ary[-1] << ";"
          end
          if /\)$/ =~ ary[-1] && /^\(/ =~ line
            ary[-1] << ";"
          end
          ary[-1] << " " << line.chomp
        end
      elsif /^(?:The d|D)imension must be at least (.*)$/ =~ line
        dim = $1
        ary[-1].sub!(/\.$/," ")
        ary[-1] << "(#{dim})"
        flag = false
      else
        flag = false
      end
    elsif flag1
      line.sub!(/^\*\s*/,"")
      ary[-1] << " " << line.chomp
    else
      if /^\*\s+If L.?WORK = -1, .*a workspace query/ =~ line
        line.sub!(/^\*\s*/,"")
        ary[-1] << " " << line.chomp
        flag1 = true
      end
    end
  }

  if @debug
    pp ary
  end

  args = Hash.new
  ary.each{|line|
    line.strip!
    /^\*\s+([A-Z\d_]+)\s+\(([^\)]+)\)\s*(.*)$/ =~ line
    name = $1.downcase
    intent = $2
    type = $3.sub(/\.$/,"")

    hash =  Hash.new

    intent = "input" if intent == "in"
    intent = "input or input/output" if intent == "input or input/ouptut"
    case sub_name
    when /^[cz]larcm$/
      case name
      when "c"
        intent = "output"
      end
    when /^[cz]lacrm$/
      case name
      when "c"
        intent = "output"
      end
    end
    hash[:intent] =  intent

    subst = Hash.new
    if /^l.*work/ =~ name && (/The (dimension|length)? of (the )?(array|work|WORK)/ =~ type || /The amount of workspace available/ =~ type)
      if (/^(.*?) The (?:dimension|length)? of (?:the )?array\s+(?:WORK\.\s*)?(.+)/ =~ type) || (/^(.*?) The (?:dimension|length)? of (?:the )?(?:array|work|WORK)\.?\s+(.+)/ =~ type)
        type = $1.strip
        str = $2
      elsif /^(.+?) The amount of workspace available/ =~ type
        type = $1.strip
        str = nil
      else
        raise "invalid:  #{type}, #{name}, #{sub_name}"
      end
      aname = name.sub(/^l/,"")
      fff = false
      if (/If #{name.upcase} = -1, (then )?.*a workspace query/ =~ str) || (/If #{name.upcase} = -1 .+? returns the optimal/ =~ str)
        fff = true
      elsif /^([cz]tgsna|[cz]hetrf)$/=~sub_name && name=="lwork"
        fff = true
      elsif str
        unless /^[sd]tgex2$/ =~ sub_name && name == "lwork"
          raise "invalid #{str} #{name} #{sub_name}"
        end
      end
      if fff
        unless ar = args[aname]
          raise "array not found #{aname} #{name} #{sub_name}"
        end
        d = ar[:dims]
        unless d.length == 1
          raise "array length is not 1 #{aname} #{name} #{sub_name}"
        end
        if /^[a-z\d_]+$/ =~ d[0]
          d[0] = "MAX(1,#{d[0]})"
        end
      end
    elsif /^ld/ =~ name && (/(leading|first) dimension of/ =~ type || /(Leading|First) dimension of/ =~ type)
      if /^([^,]+) T[Hh]e (?:leading|first) dimension of\s+(.+)$/ =~ type || /^(.+?) LD[A-Z]+ is the leading dimension of\s+(.+)$/ =~ type || /^(.+?) On entry, (?:LD[A-Z]+\s+specifies the )?(?:leading|first) dimension of (.+)$/ =~ type || /^([^,]+) (?:[Ll]eading|First) dimension of\s+(.+)$/ =~ type
        type = $1.strip
        str = $2.strip
      elsif /^([^,]+)\,\s+(.+?)\s+The (?:leading|first) dimension of\s+(.+)$/ =~ type
        type = $1.strip
        str0 = $2.strip
        str1 = $3.strip
        str = str1 + ", " + str0
      else
        raise "invalid #{type}, #{name} #{sub_name}"
      end
      if /^(?:the )?arrays (.+?),? and ([A-Z\d_]+)[\,\.]?\s*(.*)$/ =~ str
        anames = $1.strip
        aname1 = $2
        str = $3.strip
        anames = anames.split(",")
        anames.push aname1
        anames.collect!{|an| an.downcase}
      elsif /^([A-Z\d_]+) and ([A-Z\d_]+)\, must be at least (.*)$/ =~ str
        anames = [$1.downcase, $2.downcase]
        str = $3.strip
      elsif /^(?:the )?(?:array |matrix )?([A-Z\d_]+)\.?\,?\s*(.*)$/ =~ str
        anames = [$1.strip.downcase]
        str = $2.strip
        case sub_name
        when /^[cz]laqr[23]$/, /^[sd]laqr3$/, /^[sd]laqr2$/
          case name
          when "ldwv"
            anames = ["wv"]
          end
        when /^[sdcz]pbsvx$/, /^[sdcz]gbbrd$/, /^[sdcz]pbequ$/
          case name
          when "ldab"
            anames = ["ab"]
          end
        end
      else
        raise "error #{str} #{name} #{sub_name}"
      end
      ff0 = false
      ff1 = true
      anames.each{|an|
        if (ar = args[an])
          ff0 = true
          if /input/ =~ ar[:intent]
            ff1 = false
            break
          end
        end
      }
      unless ff0
        raise "arg not found [#{anames.join(",")}], #{name}, #{sub_name}"
      end
      if ff1 && str != ""
        if anames.length == 1
          aname = anames[0]
        elsif anames.include?(name.sub(/^ld/,""))
          aname = name.sub(/^ld/,"")
        else
          case sub_name
          when /^[sd]lasda$/
            case name
            when "ldgcol"
              aname = "givcol"
            end
          when /^[sd]lasd6$/
            case name
            when "ldgnum"
              aname = "givnum"
            end
          end
        end
        unless aname
          raise "cannot select anames [#{anames.join(",")}], #{name}, #{sub_name}"
        end
        un = name.upcase
        str.sub!(/\.$/,"")
        str.gsub!(/\.GE\./,"=")
        str.gsub!(/>=/,"=")
        str.gsub!(/=\s*>/,"=")
        str.gsub!(/(It )?must be at least/, "")
        str.gsub!(/#{un}\s*=/," ")
        str.sub!(/(just )?as declared in the (in the )?calling (subroutine|procedure)\./,"")
        str.strip
        if  /^([^;]+); if ([^,]+), ([^;]+); if ([^,]+), (.+)$/ =~ str
          if @@debug
            p 100
          end
          v = get_vname($1)
          cond0 = get_cond($2)
          v0 = get_vname($3)
          cond1 = get_cond($4)
          v1 = get_vname($5)
          str = "#{cond0} ? #{v0} : #{cond1} ? #{v1} : #{v}"
        elsif /^If ([^\,]+)\,\s+([^;]+); if ([^\,]+)\,\s+(.*)$/ =~ str || /^If ([^\,]+)\,\s+([^\.]+)\. If ([^\,]+)\,\s+(.*)$/ =~ str
          if @@debug
            p 110
          end
          v0 = get_vname($1)
          cond0 = get_cond($2)
          v1 = get_vname($3)
          cond1 = get_cond($4)
          str = "#{cond0} ? #{v0} : #{cond1} ? #{v1} : 0"
        elsif /^([^\,]+)\, and if ([^\,]+)\, (?:then )?(.*)$/ =~ str || /^([^;]+); if ([^\,]+)\, (?:then )?(.*)$/ =~ str
          if @@debug
            p 120
          end
          v0 = get_vname($1)
          cond1 = get_cond($2)
          v1 = get_vname($3)
          str = "#{cond1} ? #{v1} : #{v0}"
        elsif /^(.+?) if ([^;]+); (.+?) otherwise$/ =~ str
          if @@debug
            p 130
          end
          v0 = get_vname($1)
          cond0 = get_cond($2)
          v = get_vname($3)
          str = "#{cond0} ? #{v0} : #{v}"
        elsif /^If ([^,]+), then\s+(.+?)\.\s+In any case, (.+)$/ =~ str
          if @@debug
            p 140
          end
          cond0 = get_cond($1)
          v0 = get_vname($2)
          v = get_vname($3)
          str = "#{cond0} ? #{v0} : #{v}"
        elsif  /^([^;]+); and if ([^,]+), (.+)$/ =~ str
          if @@debug
            p 150
          end
          v = get_vname($1)
          cond0 = get_cond($2)
          v0 = get_vname($3)
          str = "#{cond0} ? #{v0} : #{v}"
        elsif /^(.+?) \.LE\. #{un}$/ =~ str
          if @@debug
            p 160
          end
          str = get_vname($1)
        elsif /^If ([^,]+), then\s+(.+)$/ =~ str
          if @@debug
            p 170
          end
          cond0 = get_cond($1)
          v0 = get_vname($2)
          str = "#{cond0} ? #{v0} : 0"
        elsif (/^[sd]bdsdc$/ =~ sub_name && /^ld(u|vt)$/ =~ name)
          if @@debug
            p 180
          end
          str = 'lsame_(&compq,"I") ? MAX(1,n) : 0'
        else
          if @@debug
            p 190
          end
          if /^[sdcz]laqr[23]$/ =~ sub_name && name == "ldwv"
            str = "nw"
          end
          begin
            str = get_vname(str)
          rescue
            raise "error #{str}, #{name}, #{sub_name}"
          end
        end
        if /^[sdcz]larrv$/ =~ sub_name && name == "ldz"
          str = "n"
        end
        subst[name] = str
        sss = (SUBSTS[sub_name] ||= Hash.new)
        ss = (sss[aname] ||= Hash.new)
        ss.update subst
        if (so = sss[:order]) && (soa = so[aname])
          soa += subst.keys
        end
      end
    elsif /^(.*?) array of size (.*)$/ =~ type || /^(.*?) arrays?,?(.*)$/i =~ type
      type = $1.strip
      str = $2.strip
      if /^([A-Z\s]+) work$/ =~ type
        type = $1.strip
      end
      d = DIMS[sub_name]
      dims = d[name] if d

      unless dims
        str.gsub!(/(#{name}) must be at least\s+/,'\\1 = ')
        str.gsub!(/must be at least\s+/,"")
        str.gsub!(/at least\s*/,"")
        str.gsub!(/dimension ([^\(]+) if/, '(\\1) if')
        str.gsub!(/dimension \>= /,"")
        str.gsub!(/the dimension of #{name.upcase}/,"")
        str.gsub!(/.?\s+dimension is/,"")
        str.gsub!(/(of )?dimensions?/,"")
        str.gsub!(/\.GE\.?/,">=")
        str.gsub!(/\.GT\.?/,">")
        str.gsub!(/\.LE\.?/,"<=")
        str.gsub!(/\.LT\.?/,"<")
        str.gsub!(/\.TRUE\.?/,"TRUE")
        str.gsub!(/\.FALSE\.?/,"FALSE")
        str.gsub!(/log_2\s*\(/,"1.0/log(2.0)*log((double)")
        str.gsub!(/INT\s*\(/,"(int)\(")

        str.strip!
        if intent == "input"
          str.sub!(/\s*if .*$/i, "")
        end
        if /[Ii]f/ =~ str || /where/ =~ str || /when/ =~ str
          dims = Array.new
          flag = true
          str.sub!(/;$/,"")
          if /^\([^;]+\); (\(.+?\) if [A-Z\d]+ = .+? or \(.+?\) if [A-Z\d]+ = .+)$/ =~ str
            str = $1
          end
          if @@debug
            p name
            p str
          end
          while (str && str!="")
            if (/^(?:and )?\((.*?)\)\s+if\s+([^\,\;]+)[\,\;]\s+(.*)$/ =~ str) && (dim = $1.strip) && (cond = $2.strip) && (str_tmp = $3.strip) && (/ if .* if / !~ cond) && (/\([^\)]+$/ !~ cond) && (/=/ !~ dim)
              if @@debug
                p 1
              end
              str = str_tmp
            elsif (/^\(([^;]+?)\)\s+if\s+(.+?)\s+(?:or|and)\s+(\(.*\) if.*)$/ =~ str) || (/^\(([^;]+?)\)\s+if\s+(.+?)\s+(\(.*\)\s+if.*)$/ =~ str)
              if @@debug
                p 2
              end
              dim = $1
              cond = $2.strip
              str = $3.strip
            elsif /^If ([^\,]+)\,\s+([^;\.]+)[;\.]?(.*)$/ =~ str
              if @@debug
                p 2.5
              end
              cond = $1.strip
              dim = $2
              str = $3.strip
            elsif /^\((.*?)\)\s+if\s+([^\,]+)\, (\(.*\) otherwise)$/ =~ str
              if @@debug
                p 3
              end
              dim = $1.strip
              cond = $2.strip
              str = $3.strip
            elsif (/^(?:and )?\((.+?)\)\s+if\s+(.*)$/ =~ str) && (dim = $1.strip) && (cond = $2.strip) && (/if/ !~ dim) && (/ if /i !~ cond) && (/^[^\(]+\)/ !~ dim)
              if @@debug
                p 4
              end
              str = nil
            elsif /^\((.*?)\)\s+otherwise$/ =~ str && (dim = $1.strip) && (/ if / !~ dim)
              if @@debug
                p 5
              end
              if cond == "not referenced"
                dims = dims.collect{|d| d + "0"}
              else
                get_dims(dim).each_with_index{|d,i|
                  dims[i] = dims[i] + d
                }
              end
              flag = false
              break
            elsif /^\((.+?)\);?\s+If ([A-Z_\d]+ = [^\,]+)\, (then )?#{name.upcase} is not referenced.?$/ =~ str
              if @@debug
                p 6
              end
              dims = get_dims($1)
              cond = get_cond($2)
              dims = dims.collect{|dim|
                "#{cond} ? 0 : #{dim}"
              }
              flag = false
              break
            elsif /^\((.+?)\); Not referenced if (.*)$/ =~ str
              if @@debug
                p 6.5
              end
              v0 = get_vname($1)
              cond0 = get_cond($2)
              dims = ["#{cond0} ? 0 : #{v0}"]
              flag = false
              break
            elsif /\((.+?)\)\, where ([A-Z\d]+) >= ([A-Z\d]+) when ([^;]+); otherwise, #{name.upcase} is not referenced\.?$/ =~ str
              if @@debug
                p 7
              end
              dims = get_dims($1)
              c0 = $2
              c1 = $3
              cond = get_cond($4)
              subst[get_vname(c0)] = "#{cond} ? #{get_vname(c1)} : 0"
              flag = false
              break
            elsif /\((.+)\)[\.\,] where ([A-Z\d]+) = (.+)\.?$/ =~ str
              if @@debug
                p 8
              end
              dims = get_dims($1)
              c0 = $2
              c1 = $3
              subst[get_vname(c0)] = get_vname(c1)
              flag = false
              break
            elsif /^(?:and )?(not referenced) if (.+).?$/ =~ str
              if @@debug
                p 9
              end
              dim = $1
              cond = $2
              str = nil
            elsif (/^\((.+?)\) where ([A-Z\d_]+) = (.*)$/ =~ str) && (dim = $1) && (c0 = $2) && (c1 = $3) && /=/ !~ c1
              dims = get_dims(dim)
              subst[get_vname(c0)] = get_vname(c1)
              flag = false
              break
            else
              if /^\((.+?)\) (.+?) = (.+?) (?:when|if) ([^,]+), and (.+?) when (.+)$/ =~ str
                if @@debug
                  p 10
                end
                dim = $1
                dim1 = $2
                unless dim == dim1
                  raise "error #{name} #{sub_name}"
                end
                c0 = $3
                cond0 = get_cond($4)
                c1 = $5
                cond1 = get_cond($6)
                dims = get_dims($1)
                subst[get_vname(dim)] = "#{cond0} ? #{get_vname(c0)} : #{cond1} ? #{get_vname(c1)} : 0"
                flag = false
                break
              elsif /^\((.+?)\) (.+?) = (.+?) (?:when|if) ([^,]+), and (.+?) otherwise.?$/ =~ str
                if @@debug
                  p 11
                end
                dim = $1
                dim1 = $2
                unless dim == dim1
                  raise "error #{name} #{sub_name}"
                end
                c0 = $3
                cond0 = get_cond($4)
                c1 = $5
                dims = get_dims(dim)
                subst[get_vname(dim)] = "#{cond0} ? #{get_vname(c0)} : #{get_vname(c1)}"
                flag = false
                break
              elsif /^\((.+?)\);\s+[Ii]f ([A-Z_\d]+ = [^\,]+)\, ([A-Z_\d]+) \>?= ([^\.]+)\.\s+Otherwise\, ([A-Z_\d]+) \>?= (.*)$/ =~ str
                if @@debug
                  p 14
                end
                dim = $1
                cond0 = $2
                c00 = $3.downcase
                c01 = $4.downcase
                c10 = $5.downcase
                c11 = $6.downcase
                unless c00==c10
                  raise "error #{name} #{sub_name}"
                end
                dims = get_dims(dim)
                cond0 = get_cond(cond0)
                subst[get_vname(c00)] = "#{cond0} ? #{get_vname(c01)} : #{get_vname(c11)}"
                flag = false
                break
              end
              if /^\((.+?)\)\.\s+[Ii]f ([A-Z_\d]+ = [^\,]+)\, ([A-Z_\d]+) = ([A-Z_\d]+); if ([A-Z_\d]+ = [^\,]+), ([A-Z_\d]+) = ([A-Z_\d]+)$/ =~ str
                if @@debug
                  p 13
                end
                dim = $1
                cond0 = $2
                c00 = $3.downcase
                c01 = $4.downcase
                cond1 = $5
                c10 = $6.downcase
                c11 = $7.downcase
              elsif /\((.+?)\);? ([A-Z_\d]+) = (.+?) if ([^;]+); ([A-Z_\d]+) = (.+?) if (.+)$/ =~ str
                if @@debug
                  p 13
                end
                dim = $1
                c00 = $2.downcase
                c01 = $3.downcase
                cond0 = $4
                c10 = $5.downcase
                c11 = $6.downcase
                cond1 = $7
              else
                raise "error '#{str}' in #{name} (#{fname})"
              end
              dims = get_dims(dim)
              cond0 = get_cond(cond0)
              cond1 = get_cond(cond1)
              if c00 == c10
                subst[get_vname(c00)] = "#{cond0} ? #{get_vname(c01)} : #{cond1} ? #{get_vname(c11)} : 0"
              else
                raise "error #{name} #{sub_name}"
              end
              flag = false
              break
            end
            if flag
              cond = get_cond(cond)
              if dim == "not referenced"
                dims = dims.collect{|d| d + cond + " ? 0 : "}
              else
                ds = get_dims(dim)
                ds.each_with_index{|d,i|
                  dims[i] ||= ""
                  dims[i] << cond + " ? " + d + " : "
                }
              end
            end
          end
          dims = dims.collect{|d| d + "0"} if flag
        else
          dims = get_dims(str)
        end

        dims.each_with_index{|dim,i|
          if /\?/ =~ dim
            str = dim
            aa = Array.new
            while( /^[^\?]+ \? ([^:]+) : (.+)$/ =~ str )
              aa.push $1
              str = $2
            end
            aa.push str
            if aa.length > 2
              aa = aa.uniq
              if aa.length == 1 || (aa.length == 2 && aa[-1] == "0")
                dims[i] = aa[0]
              end
            elsif aa.length == 2
              if aa.uniq.length == 1
                dims[i] = aa[0]
              end
            else
              raise "error [#{aa.join(",")}] #{dim} #{name} #{sub_name}"
            end
          end
        }
      end
      sss = (SUBSTS[sub_name] ||= Hash.new)
      ss = (sss[name] ||= Hash.new)
      ss.update subst
      if (so = sss[:order]) && (soa = so[name])
        soa += subst.keys
      end
      if /^[sd]lasda$/ =~ sub_name  && name == "difl"
        ss["nlvl"].sub!(/\)$/,"")
      end
      hash[:dims] = dims
    elsif (/^(.+) work array$/ =~ type) || (/^(.+) array$/ =~ type)
      type = $1.strip
      dims = DIMS[sub_name] && DIMS[sub_name][name]
      unless dims
        raise "dimension is not defined (#{name} in #{fname})"
      end
    elsif /^CHARACTER\*(\d+)$/ =~ type || /^CHARACTER\*\((\d+)\)$/ =~ type || /^character string/ =~ type
      type = "CHARACTER"
      if $1 && $1 != "1"
        hash[:dims] = [$1]
      end
    elsif /^CHARACTER\*\(\*\)$/ =~ type
      type = "CHARACTER"
      hash[:dims] = ["*"]
    elsif /^(.*)\,\s*#{name}\s*=\s*>/i =~ type
      type = $1.strip
    elsif /^(.+?)\s+scalar$/ =~ type
      type = $1
    elsif /^(.+?)\s+with value/ =~ type
      type = $1
    end
    if (t = TYPES[sub_name]) && (t = t[name])
      hash[:type] = t
    else
      type.sub!(/ scalar/,"")
      if /^(.+?) FUNCTION of ([^\s]+) (.+?) arguments?$/ =~ type
        hash[:block_type] = CTYPES[$1] || raise("error block type is invalid #{$1} #{name} #{sub_name}")
        hash[:block_arg_num] = {"one"=>1,"two"=>2,"three"=>3}[$2] || raise("error #{$2}")
        hash[:block_arg_type] = CTYPES[$3] || raise("error block arg type is invalid #{$3} #{name} #{sub_name}")
      else
        hash[:type] = CTYPES[type.upcase] || raise("type (#{type}) is not defined in #{name} (#{fname})")
      end
    end
    args[name] = hash
  }
  case sub_name
  when /^[cz]laqr[04]$/
    args["iloz"] = {:type => "integer", :intent => "input"}
    args["ihiz"] = {:type => "integer", :intent => "input"}
  end
  if @@debug
    pp args
  end

  return  {:sub_name => sub_name, :sub_type => sub_type, :func_type => func_type, :arg_names => arg_names, :args => args, :help => help}
end


def create_code(fname)
=begin
  if IGNORE.include?( File.basename(fname).sub(/\.\w+$/,"") )
    warn "skip #{fname}"
    return nil
  end
=end
  hash = parse_file(fname)
  if hash.nil?
    warn "skip #{fname}"
    return nil
  end
  sub_name = hash[:sub_name]
  sub_type = hash[:sub_type]
  func_type = hash[:func_type]
  arg_names = hash[:arg_names]
  args = hash[:args]
  help = hash[:help]

  unless sub_name
    warn "this has no subroutine (#{fname})"
    return nil
  end
  unless arg_names
    raise "no arg_names (#{fname})"
  end
  unless args && !args.empty?
    raise "no args (#{fname})"
  end

  if ss = SUBSTS[sub_name]
    so = ss[:order]
    if so
      so.each{|k,v|
        args[k][:subst_order] = v
      }
    end
  end
  args.each{|name,arg|
    arg[:subst] = (ss && ss[name]) || Hash.new
  }



  inputs = Array.new
  outputs = Array.new
  inouts = Array.new
  workspaces = Array.new
  block = nil
  arg_names.each{|name|
    arg = args[name]
    unless arg
      if ARGS[sub_name] && arg = ARGS[sub_name][name]
        args[name] = arg
      else
        raise "arg #{name} is not defined (#{fname})"
      end
    end
    case arg[:intent]
    when "input"
      inputs.push name
    when "output"
      outputs.push name
    when "input/output", "input or output", "input or input/output"
      inputs.push name
      inouts.push name
    when "workspace"
      workspaces.push name
    when "workspace/output"
      outputs.push name
    when "external procedure"
      if block
        raise "only one block is supported"
      end
      block = name
    else
      raise "intent (#{arg[:intent]}) is invalid (#{name}) #{sub_name}"
    end
  }
  args.each{|key,arg|
    dims = arg[:dims]
    next unless dims
    dims.each{|dim|
      inputs.delete(dim)
    }
  }

  if @@debug
    p "inputs"
    p inputs
    p "outputs"
    p outputs
    p "inouts"
    p inouts
    p "workspaces"
    p workspaces
    p "block"
    p block
  end

  code = ""

  if sub_type == :func
    outputs.push "__out__"
    args["__out__"] = {:type => func_type}
    unless sub_name == "lsame" || sub_name == "ilaenv" || sub_name == "zladiv"
      code += "extern VOID #{sub_name}_(#{func_type} *__out__, #{arg_names.collect{|an|t = args[an][:type];t+' *'+an}.join(', ')});\n"
    end
  end

  code +=<<"EOF"
static VALUE
#{RBPREFIX}#{sub_name}(int argc, VALUE *argv, VALUE self){
EOF

  dimdefs = Array.new

  (inputs+outputs).each{|name|
    arg = args[name]
    type = arg[:type]
    dims = arg[:dims]
    code +=<<"EOF"
  VALUE #{RBPREFIX}#{name};
  #{type} #{dims ? "*" : ""}#{name}; 
EOF
    dimdefs.push name
  }
  inouts.each{|name|
    arg = args[name]
    dims = arg[:dims]
    if dims
      type = arg[:type]
      code +=<<"EOF"
  VALUE #{RBPREFIX}#{name}_out__;
  #{type} *#{name}_out__;
EOF
    end
  }

  workspaces.each{|name|
    arg = args[name]
    type = arg[:type]
    dims = arg[:dims]
    code << "  #{type} #{dims ? "*" : ""}#{name};\n"
  }
  debug_dims = Array.new
  code << "\n"
  (inputs+outputs+workspaces).each{|name|
    arg = args[name]
    if dims = arg[:dims]
      dims.each{|dim|
        if /^[a-z][a-z_\d]*$/ !~ dim
          next
        end
        unless dimdefs.include?(dim)
          code << "  integer #{dim};\n"
          dimdefs.push dim
        end
      }
    end
    if ss = arg[:subst]
      ss.each{|k,v|
        unless dimdefs.include?(k)
          code << "  integer #{k};\n"
          dimdefs.push k
        end
      }
    end
  }


  code << "\n"

  if block
    block_help = "{|" + ['a','b','c'][0...args[block][:block_arg_num]].join(",") + "| ... }"
  else
    block_help = ""
  end

  help =<<"EOF"
USAGE:
  #{(outputs+inouts).join(", ")} = NumRu::Lapack.#{sub_name}( #{inputs.join(", ")})#{block_help}
    or
  NumRu::Lapack.#{sub_name}  # print help


FORTRAN MANUAL
#{help}
EOF
  ilen = inputs.length
  code +=<<"EOF"
  if (argc == 0) {
    printf("%s\\n", "#{help.gsub(/\\/,'\\\\\\').gsub(/\n/,'\n').gsub(/"/,'\"')}");
    return Qnil;
  }
  if (argc != #{ilen})
    rb_raise(rb_eArgError,"wrong number of arguments (%d for #{ilen})", argc);
EOF
  inputs.each_with_index{|arg,i|
    code << "  #{RBPREFIX}#{arg} = argv[#{i}];\n"
  }
  code << "\n"

  case sub_name
  when /^[scdz]latzm$/
    inputs.sort!{|a,b| a=="c1" ? 1 : -1}
  when /^[sdcz]spsv$/, /^[cz]tprfs$/, /^[cz]hpsv$/
    inputs.sort!{|a,b| (a=="ap" && b=="b") ? 1 : -1}
  when /^[cz]upgtr$/, /^[sd]opgtr$/
    inputs.sort!{|a,b| (a=="ap" && b=="tau") ? 1 : -1}
  when /^[sdcz]latps$/
    inputs.sort!{|a,b| (a=="ap" && b=="x") ? 1 : -1}
  when /^[sdcz]ppsvx$/, /^[sdcz]laqsp$/, /^[cz]laqhp$/
    inputs.sort!{|a,b| ((a=="ap"||a=="afp") && b=="s") ? 1 : -1}
  when /^[sdcz]sprfs$/, /^[sdcz]sptrs$/, /^[cz]hpcon$/, /^[cz]hptri$/, /^[sdcz]spcon$/, /^[cz]hpsvx$/, /^[sdcz]sptri$/, /^[sdcz]spsvx$/, /^[cz]hptrs$/, /^[cz]hprfs$/
    inputs.sort!{|a,b| ((a=="ap"||a=="afp") && b=="ipiv") ? 1 : -1}
  when /^[sdcz]pprfs$/
    inputs.sort!{|a,b| ((a=="ap"||a=="afp") && b=="b") ? 1 : -1}
  when /^[sdcz]gtcon$/, /^[sdcz]gttrf$/, /^[sdcz]gtts2$/, /^[sdcz]gttrs$/, /^[sdcz]gtsvx$/, /^[sdcz]gtsv$/, /^[sdcz]gtrfs$/, /^[sdcz]lagtm$/, /^[sdcz]langt$/
    inputs.sort!{|a,b| (a=="dl" && b=="d") ? 1 : -1}
  when /^[sdcz]larzt$/, /^[sdcz]larft$/
    inputs.sort!{|a,b| (a=="v" && (b=="tau"||b=="n")) ? 1 : -1}
  when /^[sd]laeda$/
    inputs.sort!{|a,b| a=="qptr" ? -1 : 1}
  when /^[sd]larrd$/
    inputs.sort!{|a,b| a=="d" ? -1 : 1}
  when /^[cz]syr$/
    inputs.sort!{|a,b| (a=="x" && b=="a") ? 1 : -1}
  end

  dimdefs = Hash.new
  substs = Array.new
  inputs.each_with_index{|name,i|
    arg = args[name]
    type = arg[:type]
    dims = arg[:dims]
    subst = arg[:subst]
    next if dims
    code << get_cobj(name, type);
  }
  inputs.each_with_index{|name,i|
    arg = args[name]
    type = arg[:type]
    dims = arg[:dims]
    subst = arg[:subst]
    next unless dims
    if type == "char"
      code << "  #{name} = StringValueCStr(#{RBPREFIX}#{name});\n"
      next
    end
    code +=<<"EOF"
  if (!NA_IsNArray(#{RBPREFIX}#{name}))
    rb_raise(rb_eArgError, "#{name} (#{i+1}th argument) must be NArray");
  if (NA_RANK(#{RBPREFIX}#{name}) != #{dims.length})
    rb_raise(rb_eArgError, "rank of #{name} (#{i+1}th argument) must be %d", #{dims.length});
EOF
    if sa = arg[:subst_order]
      sa.each{|k|
        v = subst[k]
        code << "  #{k} = #{v};\n"
        substs.push k
      }
    else
      subst.each{|k,v|
        code << "  #{k} = #{v};\n"
        substs.push k
      }
    end

    dims.each_with_index{|dim,j|
      raise "bug: NA_SHAPE? cannot use {#{dim} in #{name}: #{sub_name}" if j>2
      if substs.include?(dim) || dimdefs[dim] == true
        code += <<"EOF"
  if (NA_SHAPE#{j}(#{RBPREFIX}#{name}) != #{dim})
    rb_raise(rb_eRuntimeError, "shape #{j} of #{name} must be #{dim}");
EOF
      elsif (dimdef = dimdefs[dim])
        code += <<"EOF"
  if (NA_SHAPE#{j}(#{RBPREFIX}#{name}) != #{dim})
    rb_raise(rb_eRuntimeError, "shape #{j} of #{name} must be the same as shape #{dimdef[:index]} of #{dimdef[:name]}");
EOF
      elsif /^[a-z][a-z_\d]*$/ !~ dim
        get_vars(dim).each{|d|
          unless dimdefs[d] || substs.include?(d) || ((inputs.include?(d)) && (ar = args[d]) && (ar[:type]=="integer"||ar[:type]=="logical"))
            raise "undefined #{d}  #{name} #{sub_name}"
          end
        }
        code += <<"EOF"
  if (NA_SHAPE#{j}(#{RBPREFIX}#{name}) != (#{dim}))
    rb_raise(rb_eRuntimeError, "shape #{j} of #{name} must be %d", #{dim});
EOF
      else
        code << "  #{dim} = NA_SHAPE#{j}(#{RBPREFIX}#{name});\n"
        dimdefs[dim] = {:name => name, :index => j}
      end
    }
    natype =  NATYPES[type] || raise("na type is not deifned (#{type})")
    code +=<<"EOF"
  if (NA_TYPE(#{RBPREFIX}#{name}) != #{natype})
    #{RBPREFIX}#{name} = na_change_type(#{RBPREFIX}#{name}, #{natype});
  #{name} = NA_PTR_TYPE(#{RBPREFIX}#{name}, #{type}*);
EOF
  }

  outputs.each{|name|
    arg = args[name]
    type = arg[:type]
    if dims = arg[:dims]
      if ss = arg[:subst]
        if sa = arg[:subst_order]
          sa.each{|k|
            v = ss[k]
            code << "  #{k} = #{v};\n"
            substs.push k
          }
        else
          ss.each{|k,v|
            code << "  #{k} = #{v};\n"
            substs.push k
          }
        end
      end
      code +=<<"EOF"
  {
    int shape[#{dims.length}];
EOF
      dims.each_with_index{|dim,k|
        get_vars(dim).each{|d|
          unless dimdefs[d] || substs.include?(d) || ((inputs.include?(d)) && (ar = args[d]) && (ar[:type]=="integer"||ar[:type]=="logical"))
            raise "undefined #{d}  #{name} #{sub_name}"
          end
        }
        code << "    shape[#{k}] = #{dim};\n"
      }
      code +=<<"EOF"
    #{RBPREFIX}#{name} = na_make_object(#{NATYPES[type]}, #{dims.length}, shape, cNArray);
  }
  #{name} = NA_PTR_TYPE(#{RBPREFIX}#{name}, #{type}*);
EOF
    end
  }

  inouts.each{|name|
    arg = args[name]
    type = arg[:type]
    if dims = arg[:dims]
      code +=<<"EOF"
  {
    int shape[#{dims.length}];
EOF
      dims.each_with_index{|dim,k|
        get_vars(dim).each{|d|
          unless dimdefs[d] || substs.include?(d) || ((inputs.include?(d)) && (ar = args[d]) && (ar[:type]=="integer"||ar[:type]=="logical"))
            raise "undefined #{d}  #{name} #{sub_name}"
          end
        }
        code << "    shape[#{k}] = #{dim};\n"
      }
      code +=<<"EOF"
    #{RBPREFIX}#{name}_out__ = na_make_object(#{NATYPES[type]}, #{dims.length}, shape, cNArray);
  }
  #{name}_out__ = NA_PTR_TYPE(#{RBPREFIX}#{name}_out__, #{type}*);
  MEMCPY(#{name}_out__, #{name}, #{type}, NA_TOTAL(#{RBPREFIX}#{name}));
  #{RBPREFIX}#{name} = #{RBPREFIX}#{name}_out__;
  #{name} = #{name}_out__;
EOF
    end
  }

  workspaces.each{|name|
    arg = args[name]
    if dims = arg[:dims]
      type = arg[:type]
      if ss = arg[:subst]
        if sa = arg[:subst_order]
          sa.each{|k|
            v = ss[k]
            code << "  #{k} = #{v};\n"
            substs.push k
          }
        else
          ss.each{|k,v|
            code << "  #{k} = #{v};\n"
            substs.push k
          }
        end
      end
      dims.each{|dim|
        get_vars(dim).each{|d|
          unless dimdefs[d] || substs.include?(d) || ((inputs.include?(d)) && (ar = args[d]) && (ar[:type]=="integer"))
            raise "undefined #{d}  #{name} #{sub_name}"
          end
        }
      }
      len = dims.collect{|dim| "(#{dim})"}.join("*")
      code << "  #{name} = ALLOC_N(#{type}, #{len});\n"
    end
  }
  code << "\n"


  cargs = arg_names.collect{|name|
    block== name ? RBPREFIX+name : args[name][:dims] ? name : "&"+name
  }
  if sub_name == "ilaenv"
    code << "  __out__ = #{sub_name}_(#{cargs.join(", ")}, (ftnlen) strlen(name), (ftnlen)1);\n\n"
  else
    code << "  #{sub_name}_(#{sub_type==:func ? "&__out__, " : ""}#{cargs.join(", ")});\n\n"
  end

  workspaces.each{|name|
    arg = args[name]
    if dims = arg[:dims]
      type = arg[:type]
      len = dims.collect{|dim| "(#{dim})"}.join("*")
      code << "  free(#{name});\n"
    end
  }

  out = outputs + inouts

  out.each{|name|
    arg = args[name]
    if arg[:dims]
      if arg[:type] == "char"
        code << "  #{RBPREFIX}#{name} = rb_str_new2(&#{name});\n"
      end
    else
      code << get_robj(name, arg[:type])
    end
  }

  case out.length
  when 0
    result = "Qnil"
  when 1
    result = RBPREFIX+out[0];
  else
    result = "rb_ary_new3(#{out.length}, #{out.collect{|op| RBPREFIX+op}.join(", ")})"
  end

  code +=<<"EOF"
  return #{result};
}

EOF

  code +=<<"EOF"
void
init_lapack_#{sub_name}(VALUE mLapack){
  rb_define_module_function(mLapack, \"#{sub_name}\", #{RBPREFIX}#{sub_name}, -1);
}
EOF

  code_all = "#include \"#{RBPREFIX}lapack.h\"\n\n"

  if block
    arg = args[block]
    type = arg[:block_type]
    atype = arg[:block_arg_type]
    anum = arg[:block_arg_num]
    cas = Array.new
    ras = Array.new
    anum.times{|n|
      cas.push "#{atype} *arg#{n}"
      ras.push "#{RBPREFIX}arg#{n}"
    }
    code_all += <<EOF
static #{type}
#{RBPREFIX}#{block}(#{cas.join(", ")}){
  VALUE #{ras.join(", ")};

  VALUE #{RBPREFIX}ret;
  #{type} ret;

EOF
    anum.times{|n|
      code_all << get_robj("arg#{n}", atype, true)
    }
    code_all += <<EOF

  rb_ret = rb_yield_values(#{anum}, #{ras.join(", ")});

EOF
    code_all << get_cobj("ret", type)
    code_all += <<EOF
  return ret;
}

EOF
  end

  code_all += code

  return [code_all, sub_name]
end

def generate_code(fnames)
  nfnames = fnames.length
  sub_names = Array.new
  fnames.each_with_index{|fname,i|
    print "#{i+1}/#{nfnames}\n" if (i+1)%100==0
    code, sub_name = create_code(fname)
    if code
      sub_names.push sub_name
      File.open(sub_name+".c","w"){|file|
        file.print code
      }
    end
  }

  File.open("#{RBPREFIX}lapack.h","w"){|file|
    file.print <<"EOF"
#include <string.h>
#include <math.h>
#include "ruby.h"
#include "narray.h"
#include "f2c.h"
#include "clapack.h"

#define MAX(a,b) a > b ? a : b
#define MIN(a,b) a < b ? a : b
#define LG(n) (int)ceil(log((double)n)/log(2.0))

extern logical lsame_(char *ca, char *cb);
extern logical ilaenv_(integer *ispec, char *name__, char *opts, integer *n1, integer *n2, integer *n3, integer *n4, ftnlen name_len, ftnlen opts_len);
extern int cunmtr_(char *side, char *uplo, char *trans, integer *m, integer *n, complex *a, integer *lda, complex *tau, complex *c__, integer *ldc, complex *work, integer *lwork, integer *info);
extern int cunmrz_(char *side, char *trans, integer *m, integer *n, integer *k, integer *l, complex *a, integer *lda, complex *tau, complex *c__, integer *ldc, complex *work, integer *lwork, integer *info);
EOF
  }


  File.open("rb_lapack.c","w"){|file|
    file.print <<"EOF"
#include "ruby.h"

EOF

    sub_names.each{|sname|
      file.print "extern void init_lapack_#{sname}(VALUE mLapack);\n"
    }

    file.print <<"EOF"

void Init_lapack(){
  VALUE mNumRu;
  VALUE mLapack;

  rb_require("narray");

  mNumRu = rb_define_module("NumRu");
  mLapack = rb_define_module_under(mNumRu, "Lapack");

EOF
    sub_names.each{|sname|
      file.print "  init_lapack_#{sname}(mLapack);\n"
    }
    file.print "}\n"
  }
end





@@debug = ARGV.delete("--debug")

dname = ARGV[0]
if File.directory?(dname)
#  fnames = %w(sgees dgees cgees zgees sgeesx dgeesx cgeesx zgeesx sgges dgges cgges zgges sggesx dggesx cggesx zggesx).collect{|n| File.join("../lapack-3.1.1/SRC",n)+".f"}
#  fnames = Dir[ File.join(dname,"*.f") ][0..10]
  fnames = Dir[ File.join(dname,"*.f") ]
elsif File.file?(dname)
  fnames = [dname]
  @@debug = true
end

generate_code(fnames)


