DataTypes = [
             ["S", "REAL"],
             ["D", "DOUBLE PRECISION"],
             ["C", "COMPLEX"],
             ["Z", "COMPLEX*16 or DOUBLE COMPLEX"],
             ["DS", "Data type in double but solving problem using single precision"],
             ["ZC", "Data type in complex*16 but solving problem using complex precision"]
]

MatrixTypes = [
              ["BD", "bidiagonal"],
              ["DI", "diagonal"],
              ["GB", "general band"],
              ["GE", "general (i.e., unsymmetric, in some cases rectangular)"],
              ["GG", "general matrices, generalized problem (i.e., a pair of general matrices)"],
              ["GT", "general tridiagonal"],
              ["HB", "(complex) Hermitian band"],
              ["HE", "(complex) Hermitian"],
              ["HG", "upper Hessenberg matrix, generalized problem (i.e a Hessenberg and a triangular matrix)"],
              ["HP", "(complex) Hermitian, packed storage"],
              ["HS", "upper Hessenberg"],
              ["OP", "(real) orthogonal, packed storage"],
              ["OR", "(real) orthogonal"],
              ["PB", "symmetric or Hermitian positive definite band"],
              ["PO", "symmetric or Hermitian positive definite"],
              ["PP", "symmetric or Hermitian positive definite, packed storage"],
              ["PT", "symmetric or Hermitian positive definite tridiagonal"],
              ["SB", "(real) symmetric band"],
              ["SP", "symmetric, packed storage"],
              ["ST", "(real) symmetric tridiagonal"],
              ["SY", "symmetric"],
              ["TB", "triangular band"],
              ["TG", "triangular matrices, generalized problem (i.e., a pair of triangular matrices)"],
              ["TP", "triangular, packed storage"],
              ["TR", "triangular (or in some cases quasi-triangular)"],
              ["TZ", "trapezoidal"],
              ["UN", "(complex) unitary"],
              ["UP", "(complex) unitary, packed storageBDbidiagonal"]
]


require "numru/lapack"
include NumRu

prefix = File.dirname(__FILE__)+"/../doc"

desc = Hash.new

methods = Lapack.singleton_methods
dts = Hash.new
DataTypes.each{|cdt, dt|
  cdt = cdt.downcase
  dmethods = Array.new
  methods.each{|m|
    dmethods.push m if /^#{cdt}/ =~ m
  }
  dmethods.each do |m|
    methods.delete m
  end
  mts = Array.new
  MatrixTypes.each{|cmt, mt|
    cmt = cmt.downcase
    reg = /^#{cdt}#{cmt}/
    ms = Array.new
    dmethods.each{|m|
      next unless reg =~ m
      ms.push m
    }
    ms.sort!
    unless ms.empty?
      mts.push [cmt,mt]
      dts[cmt] ||= Array.new
      dts[cmt].push [cdt, dt]
      File.open(File.join(prefix,"#{cdt}#{cmt}.html"),"w"){|file|
        file.print <<"EOF"
<HTML>
  <HEAD>
    <TITLE>#{dt} routines for #{mt} matrix</TITLE>
  </HEAD>
  <BODY>
    <A NAME="top"></A>
    <H1>#{dt} routines for #{mt} matrix</H1>
    <UL>
EOF
        ms.each{|m|
          file.print <<"EOF"
      <LI><A HREF=\"##{m}\">#{m}</A></LI>
EOF
        }
        file.print <<"EOF"
    </UL>

EOF
        ms.each{|m|
          file.print <<"EOF"
    <A NAME="#{m}"></A>
    <H2>#{m}</H2>
    <PRE>
EOF
          IO.popen("-") do |io|
            if io # parent
              file.print io.read
            else # child
              Lapack.send(m, :help => true)
            end
          end
          file.print <<"EOF"
    </PRE>
    <A HREF="#top">go to the page top</A>

EOF
        }
        file.print <<"EOF"
    <HR />
    <A HREF="#{cdt}.html">back to matrix types</A><BR>
    <A HREF="#{cdt}.html">back to data types</A>
  </BODY>
</HTML>
EOF
      }
    end
  }

  unless mts.empty?
    File.open(File.join(prefix,"#{cdt}.html"),"w"){|file|
      file.print <<"EOF"
<HTML>
  <HEAD>
    <TITLE>#{dt} routines</TITLE>
  </HEAD>
  <BODY>
    <H1>#{dt} routines</H1>
    <UL>
EOF
      mts.each{|cmt,mt|
        file.print "      <LI><A HREF=\"#{cdt}#{cmt}.html\">#{cmt.upcase}: #{mt}</A></LI>\n"
      }
      file.print <<"EOF"
    </UL>
    <HR />
    <A HREF="index.html">back to index.html</A>
  </BODY>
</HTML>
EOF
    }
  end
}

MatrixTypes.each do |cmt,mt|
  cmt = cmt.downcase
  if dts[cmt]
    File.open(File.join(prefix,"#{cmt}.html"),"w") do |file|
      file.print <<"EOF"
<HTML>
  <HEAD>
    <TITLE>#{mt} routines</TITLE>
  </HEAD>
  <BODY>
    <H1>#{mt} routines</H1>
    <UL>
EOF
      dts[cmt].each{|cdt,dt|
        file.print "      <LI><A HREF=\"#{cdt}#{cmt}.html\">#{cdt.upcase}: #{dt}</A></LI>\n"
      }
      file.print <<"EOF"
    </UL>
    <HR />
    <A HREF="index.html">back to index.html</A>
  </BODY>
</HTML>
EOF
    end
  end
end

if methods.any?
  File.open(File.join(prefix,"others.html"),"w") do |file|
    file.print <<EOF
<HTML>
  <HEAD>
    <TITLE>other routines</TITLE>
  </HEAD>
  <BODY>
    <A NAME="top"></A>
    <H1>other routines</H1>
    <UL>
EOF
    methods.each do |m|
      file.print <<EOF
      <LI><A HREF=\"##{m}\">#{m}</A></LI>
EOF
    end
    file.print <<EOF
    </UL>

EOF
    methods.each do |m|
      file.print <<EOF
    <A NAME="#{m}"></A>
    <H2>#{m}</H2>
    <PRE>
EOF
      IO.popen("-") do |io|
        if io # parent
          file.print io.read
        else # child
          Lapack.send(m, :help => true)
        end
      end
      file.print <<EOF
    </PRE>
    <A HREF="#top">go to the page top</A>

EOF
    end
    file.print <<"EOF"
    <HR />
    <A HREF="index.html">back to index</A>
  </BODY>
</HTML>
EOF
  end
end

File.open(File.join(prefix,"index.html"),"w"){|file|
  file.print <<"EOF"
<HTML>
  <HEAD>
    <TITLE>LAPACK routines</TITLE>
  </HEAD>
  <BODY>
    <H1>Data types</H1>
    <UL>
EOF
  DataTypes.each{|cdt,dt|
    file.print "      <LI><A HREF=\"#{cdt.downcase}.html\">#{cdt}: #{dt}</A></LI>\n"
  }
  file.print <<"EOF"
    </UL>

    <H1>Matrix types</H1>
    <UL>
EOF
  MatrixTypes.each do |cmt,mt|
    if dts[cmt.downcase]
      file.print "     <LI><A HREF=\"#{cmt.downcase}.html\">#{cmt}: #{mt}</A></LI>\n"
    end
  end
  file.print <<"EOF"
    </UL>
EOF
  if methods.any?
    file.print <<EOF

    <H1>others</H1>
    <UL>
      <LI><A HREF=\"others.html\">others</A></LI>
    </UL>
EOF
  end
  file.print <<"EOF"
  </BODY>
</HTML>
EOF
}
