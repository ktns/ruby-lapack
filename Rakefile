require "rubygems"
require "rake/clean"
require "rake/gempackagetask"
require "rake/testtask"

version = "1.5"
target_prefix = "numru"

# get options
destdir = ENV["DESTDIR"] || ""
libdir = ENV["SITELIBDIR"] || Config::CONFIG["sitelibdir"]
archdir = ENV["SITEARCHDIR"] || Config::CONFIG["sitearchdir"]
config_opts = ENV["CONFIG_OPTIONS"]

NAME = "lapack"
LIBS = FileList["lib/#{target_prefix}/*rb"]
DLLIB = "ext/#{NAME}.so"
so_file = File.join("lib", target_prefix, "#{NAME}.so")


task :default => so_file

desc "Building extensions"
file so_file => DLLIB do
  mkdir_p File.dirname(so_file)
  rm_f so_file
  cp DLLIB, so_file
end
file DLLIB => "ext/Makefile" do
  system("cd ext; make")
end
file "ext/Makefile" => "ext/rb_lapack.h" do
  unless system("cd ext; ruby extconf.rb #{config_opts}")
    warn <<-EOL

To give options to extconf.rb, set the options to CONFIG_OPTIONS
e.g.
% rake CONFIG_OPTIONS="--with-lapack=/opt/lapack"
    EOL
  end
end
file "ext/rb_lapack.h" => "dev/make_csrc.rb" do
  system("ruby dev/make_csrc.rb")
end

desc "Install files to system"
task :install => [:install_so, :install_rb]

task :install_so => DLLIB do
  dst = File.join(destdir, archdir, target_prefix)
  mkdir_p dst
  install DLLIB, dst, :mode => 0755
end

task :install_rb => LIBS do
  dst = File.join(destdir, libdir, target_prefix)
  mkdir_p dst
  LIBS.each do |lib|
    install lib, dst, :mode => 644
  end
end

CLEAN.include("ext/*.o")
CLOBBER.include(DLLIB, so_file)
CLOBBER.include("ext/Makefile")


PKG_FILES = FileList["lib/#{target_prefix}/*rb"]
PKG_FILES.include("ext/rb_lapack.h")
PKG_FILES.include("ext/f2c_minimal.h")
PKG_FILES.include("ext/*.c")
PKG_FILES.include("Rakefile")
PKG_FILES.include("COPYING", "GPL", "README.rdoc")
PKG_FILES.include("doc/*.html", "samples/**/*rb")
PKG_FILES.include("dev/*.rb", "dev/defs/*")
TEST_FILES = FileList["tests/**/*.rb"]

Rake::TestTask.new do |t|
  t.libs << "lib"
  t.libs << "tests"
  t.test_files = TEST_FILES
end

spec = Gem::Specification.new do |s|
  s.name = "ruby-lapack"
  s.version = version
  s.summary = "A Ruby wrapper of Lapack"
  s.description = <<EOL
Ruby-LAPACK is a Ruby wrapper of Lapack, which is a linear algebra package (http://www.netlib.org/lapack/).
EOL
  s.author = "Seiya Nishizawa"
  s.email = "seiya@gfd-dennou.org"
  s.homepage = "http://ruby.gfd-dennou.org/products/ruby-lapack/"
  s.has_rdoc = false
  s.files = PKG_FILES
  s.test_files = TEST_FILES
  s.add_dependency('narray')
  s.extensions = %w(ext/extconf.rb)
end


Rake::GemPackageTask.new(spec) do |pkg|
  pkg.need_tar_gz = true
  pkg.need_tar_bz2 = true
end



binary_pkg = "pkg/#{spec.name}-#{spec.version}-#{Config::CONFIG["arch"]}.gem"
gem_pkg = "pkg/#{spec.name}-#{spec.version}.gem"
desc "Build binary package"
task :binary_package => binary_pkg

file binary_pkg => gem_pkg do
  system "gem compile --fat 1.8:ruby1.8,1.9:ruby1.9 #{gem_pkg}"
end
