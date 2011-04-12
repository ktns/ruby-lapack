$:.push File.dirname(__FILE__) + "/../.."
require "lapack_test"

class GtsvTest < Test::Unit::TestCase
  include LapackTest

  def setup
    @du = Hash.new
    @d  = Hash.new
    @dl = Hash.new
    @b = Hash.new
    @b_exp = Hash.new

    @du[:r] = NArray[2.1, -1.0,  1.9,  8.0]
    @d[:r]  = NArray[3.0,  2.3, -5.0, -0.9, 7.1]
    @dl[:r] = NArray[3.4,  3.6,  7.0, -6.0]
    @b[:r] = NArray[[2.7, -0.5,  2.6,  0.6, 2.7]]
    @b_exp[:r] = NArray[[-4.0, 7.0, 3.0, -4.0, -3.0]]

    @du[:c] = NArray[ 2.0-1.0*I,  2.0+1.0*I, -1.0+1.0*I,  1.0-1.0*I]
    @d[:c]  = NArray[-1.3+1.3*I, -1.3+1.3*I, -1.3+3.3*I, -0.3+4.3*I, -3.3+1.3*I]
    @dl[:c] = NArray[ 1.0-2.0*I,  1.0+1.0*I,  2.0-3.0*I,  1.0+1.0*I]
    @b[:c] = NArray[[2.4-5.0*I, 3.4+18.2*I, -14.7+9.7*I, 31.9-7.7*I, -1.0+1.6*I]]
    @b_exp[:c] = NArray[[1.0+1.0*I, 3.0-1.0*I, 4.0+5.0*I, -1.0-2.0*I, 1.0-1.0*I]]
  end

  %w(s d c z).each do |x|
    method = "#{x}gtsv"
    rc = LapackTest.get_rc(x)

    define_method("test_#{method}") do 
      info, dl, d, du, b = NumRu::Lapack.send(method, @dl[rc], @d[rc], @du[rc], @b[rc])
      assert_equal 0, info
      assert_narray @b_exp[rc], b
    end

  end

end
