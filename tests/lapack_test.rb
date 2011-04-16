$:.unshift(File.dirname(__FILE__), "..", "lib")
require "test/unit"
require "numru/lapack"

module LapackTest

  I = Complex::I

  def assert_narray(expected, actual, delta=nil, message="")
    unless delta
      case actual.typecode
      when NArray::SFLOAT, NArray::SCOMPLEX
        delta = 1.0e-5
      when NArray::DFLOAT, NArray::DCOMPLEX
        delta = 1.0e-13
      when NArray::INT, NArray::LINT
        delta = 0
      else
        raise "typecode is invalid"
      end
    end
    if message.empty?
      message = <<EOF
<#{expected.inspect}>
and
<#{actual.inspect}>
expected to have maximan differnce <#{(expected-actual).abs.max}> within
<#{delta}>.
EOF
    end
    assert (expected - actual).abs.max <= delta, message
  end

  def get_int(x)
    x = x.real if x.respond_to?(:real)
    x.to_i
  end

  def comp_sign(a, b)
    a = a.real if a.respond_to?(:real)
    b = b.real if b.respond_to?(:real)
    a*b < 0
  end

  def get_rc(x)
    /\A[sd]/ =~ x ? :r : :c
  end
  module_function :get_rc

end
