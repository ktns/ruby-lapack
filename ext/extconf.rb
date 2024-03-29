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



dir_config("lapack")
unless find_library("lapack", nil)
  library_not_found("lapack",nil)

  warn "LAPACK will be tried to find"

  name = with_config("blas-name","blas_LINUX.a")
  unless have_library(name)
    lib_path = with_config("blas-lib","/usr/local/lib")
    _libarg = LIBARG
    LIBARG.replace "#{lib_path}/%s"
    unless have_library(name)
      library_not_found("blas",name)
    end
    LIBARG.replace _libarg
  end
  name = with_config("lapack-name","lapack_LINUX.a")
  unless have_library(name)
    lib_path = with_config("lapack-lib","/usr/local/lib")
    _libarg = LIBARG
    LIBARG.replace "#{lib_path}/%s"
    unless have_library(name)
      library_not_found("lapack",name)
    end
    LIBARG.replace _libarg
  end
end

sitearchdir = Config::CONFIG["sitearchdir"]
dir_config("narray", sitearchdir, sitearchdir)
gem_path = nil
begin
  require "rubygems"
  if (spec = Gem.source_index.find_name("narray")).any?
    gem_path = spec.last.full_gem_path
  end
rescue LoadError
end
unless find_header("narray.h",gem_path) && have_header("narray_config.h")
  header_not_found("narray")
end

create_makefile("numru/lapack")
