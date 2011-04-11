$:.push File.dirname(__FILE__) + "/../.."
require "lapack_test"

class GelsTest < Test::Unit::TestCase
  include LapackTest

  def setup
    @a = Hash.new
    @b = Hash.new
    @b_exp = Hash.new

    @a[:r] = NMatrix[[-0.57, -1.28, -0.39,  0.25],
                     [-1.93,  1.08, -0.31, -2.14],
                     [ 2.30,  0.24,  0.40, -0.35],
                     [-1.93,  0.64, -0.66,  0.08],
                     [ 0.15,  0.30,  0.15, -2.13],
                     [-0.02,  1.03, -1.43,  0.50]].to_lm
    @b[:r] = NVector[[-2.67, -0.55, 3.34, -0.77, 0.48, 4.10]]
    @b_exp[:r] = NArray[[1.5339, 1.8707, -1.5241, 0.0392]]

    i = Complex::I
    @a[:c] = NMatrix[[ 0.96-0.81*I, -0.03+0.96*I, -0.91+2.06*I, -0.05+0.41*I],
                     [-0.98+1.98*I, -1.20+0.19*I, -0.66+0.42*I, -0.81+0.56*I],
                     [ 0.62-0.46*I,  1.01+0.02*I,  0.63-0.17*I, -1.11+0.60*I],
                     [-0.37+0.38*I,  0.19-0.54*I, -0.98-0.36*I,  0.22-0.20*I],
                     [ 0.83+0.51*I,  0.20+0.01*I, -0.17-0.46*I,  1.47+1.59*I],
                     [ 1.08-0.28*I,  0.20-0.12*I, -0.07+1.23*I,  0.26+0.26*I]].to_lm
    @b[:c] = NVector[[-2.09+1.93*I, 3.34-3.53*I, -4.94-2.04*I, 0.17+4.23*I, -5.19+3.63*I, 0.98+2.53*I]]
    @b_exp[:c] = NArray[[-0.5044-1.2179*I, -2.4281+2.8574*I, 1.4872-2.1955*I, 0.4537+2.6904*I]]
  end

  %w(s d c z).each do |x|
    method = "#{x}gels"
    rc = LapackTest.get_rc(x)

    define_method("test_#{method}") do
      work, info, a, b = NumRu::Lapack.send(method, "N", @a[rc], @b[rc])
      assert_equal 0, info
      assert_narray @b_exp[rc], b, 1.0e-4
    end

    define_method("test_#{method}_inquiring_lwork") do
      work, info, a, b = NumRu::Lapack.send(method, "N", @a[rc], @b[rc], :lwork => -1)
      assert_equal 0, info
      lwork = get_int(work[0])
      work, info, a, b = NumRu::Lapack.send(method, "N", @a[rc], @b[rc], :lwork => lwork)
      assert_equal 0, info
      assert_equal lwork, get_int(work[0])
      assert_narray @b_exp[rc], b, 1.0e-4
    end

  end

end
