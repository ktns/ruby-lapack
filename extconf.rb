require "mkmf"


def header_not_found(name)
  warn <<EOF
 #{name}.h was not found.
 If you have #{name}.h, try the following:
   % ruby extconf.rb --with-#{name}-include=path
EOF
  exit 1
end

def library_not_found(lname, fname=nil)
  if fname
    warn <<EOF
  #{fname} was not found.
  If you have #{lname} library, try the following:
    % ruby extconf.rb --with-#{lname}-lib=path --with-#{lname}-name=name
  e.g.
    If you have /usr/local/#{lname}/#{fname},
     % ruby extconf.rb --with-#{lname}-lib=/usr/local/#{lname} --with-#{lname}-name=#{fname}
EOF
    exit 1
  else
    warn <<EOF
  lib#{lname}.{a|so} was not found.
  If you have lib#{lname}.{a|so}, try the following:
    % ruby extconf.rb --with-#{lname}-lib=path
EOF
    exit 1
  end
end

def find_library(lib, func, name=nil)
  func = "main" if !func or func.empty?
  ldir = with_config(lib+'-lib')
  ldirs = ldir ? Array === ldir ? ldir : ldir.split(File::PATH_SEPARATOR) : []
  $LIBPATH = ldirs | $LIBPATH
  if /\.(a|so)$/ =~ name
    libs = $libs
    $LIBPATH.each{|path|
      f = File.join(path,name)
      if File.exist?(f)
        libs = f + " " + $libs
        break
      end
    }
  else
    name = LIBARG%lib
    libs = append_library($libs, lib)
  end
  paths = {}
  checking_for "#{func}() in #{name}" do
    libpath = $LIBPATH
    begin
      until r = try_func(func, libs) or paths.empty?
        $LIBPATH = libpath | [paths.shift]
      end
      if r
        $libs = libs
        libpath = nil
      end
    ensure
      $LIBPATH = libpath if libpath
    end
    r
  end
end




dir_config("f2c", "/usr/local")
unless find_header("f2c.h")
  header_not_found("f2c")
end
unless find_library("f2c", "s_copy")
  library_not_found("f2c")
end

dir_config("clapack", "/usr/local")
unless find_header("clapack.h")
  header_not_found("clapack")
end
#unless find_library("lapack","dsyevr_")
#  library_not_found("lapack",nil)

  warn "CLAPACK will be tried to find"
  name = with_config("cblas-name","blas_LINUX.a")
  unless find_library("cblas", "f2c_dcopy", name)
    library_not_found("cblas",name)
  end
  name = with_config("clapack-name","lapack_LINUX.a")
  unless find_library("clapack", "dsyevr_", name)
    library_not_found("lapack",name)
  end
#end

sitearchdir = Config::CONFIG["sitearchdir"]
dir_config("narray", sitearchdir, sitearchdir)
unless find_header("narray.h") && have_header("narray_config.h")
  header_not_found("narray")
end
unless find_library("narray", "na_make_object", "narray.so")
  library_not_found("narray","narray.so")
end

create_makefile("numru/lapack")
