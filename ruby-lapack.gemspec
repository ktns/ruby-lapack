Gem::Specification.new do |s|
  s.name = "ruby-lapack"
  s.version = "1.3"
  s.summary = "A Ruby wrapper of Lapack"
  s.description = <<EOL
Ruby-LAPACK is a Ruby wrapper of Lapack, which is a linear algebra package (http://www.netlib.org/lapack/).
EOL
  s.author = "Seiya Nishizawa"
  s.email = "seiya@gfd-dennou.org"
  s.homepage = "http://ruby.gfd-dennou.org/products/ruby-lapack/"
  s.has_rdoc = false
  s.files = %w(lib/lapack.rb COPYING GPL README.rdoc) + Dir.glob("ext/*.c") + Dir.glob("ext/*.h") + Dir.glob("doc/*.html")
  s.test_files = Dir.glob("tests/*.rb") + Dir.glob("tests/*/*/test_*.rb")
  s.add_dependency('narray')
  s.extensions = %w(ext/extconf.rb)
end
