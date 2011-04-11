$:.push File.dirname(__FILE__) + "/../.."
require "lapack_test"

class GelsdTest < Test::Unit::TestCase
  include LapackTest

  def setup
    @a = Hash.new
    @b = Hash.new
    @b_exp = Hash.new
    @s_exp = Hash.new
    @rank_exp = Hash.new

    @a[:r] = NMatrix[[-0.09, -1.56, -1.48, -1.09,  0.08, -1.59],
                     [ 0.14,  0.20, -0.43,  0.84,  0.55, -0.72],
                     [-0.46,  0.29,  0.89,  0.77, -1.13,  1.06],
                     [ 0.68,  1.09, -0.71,  2.11,  0.14,  1.24],
                     [ 1.29,  0.51, -0.96, -1.27,  1.74,  0.34]].to_lm
    @b[:r] = NVector[[7.4, 4.3, -8.1, 1.8, 8.7]]
    @b_exp[:r] = NArray[[1.5938, -0.1180, -3.1501, 0.1554, 2.5529, -1.6730]]
    @s_exp[:r] = NArray[3.9997, 2.9962, 2.0001, 0.9988, 0.0025]
    @rank_exp[:r] = 4

    @a[:c] = NMatrix[[ 0.47-0.34*I, -0.32-0.23*I,  0.35-0.60*I,  0.89+0.71*I, -0.19+0.06*I],
                     [-0.40+0.54*I, -0.05+0.20*I, -0.52-0.34*I, -0.45-0.45*I,  0.11-0.85*I],
                     [ 0.60+0.01*I, -0.26-0.44*I,  0.87-0.11*I, -0.02-0.57*I,  1.44+0.80*I],
                     [ 0.80-1.02*I, -0.43+0.17*I, -0.34-0.09*I,  1.14-0.78*I,  0.07+1.14*I]].to_lm
    @b[:c] = NVector[[2.15-0.20*I, -2.24+1.82*I, 4.45-4.28*I, 5.70-6.25*I]]

    @b_exp[:c] = NArray[[3.9747-1.8377*I, -0.9186+0.8253*I, -0.3105+0.1477*I, 1.0050+0.8626*I, -0.2256-1.9425*I]]
    @s_exp[:c] = NArray[2.9979, 1.9983, 1.0044, 0.0064]
    @rank_exp[:c] = 3

    @rcond = 0.01
  end

  %w(s d c z).each do |x|
    method = "#{x}gelsd"
    rc = LapackTest.get_rc(x)

    define_method("test_#{method}") do
      s, rank, work, info, b = NumRu::Lapack.send(method, @a[rc], @b[rc], @rcond, -1)
      assert_equal 0, info
      lwork = get_int(work[0])
      s, rank, work, info, b = NumRu::Lapack.send(method, @a[rc], @b[rc], @rcond, lwork)
      assert_equal 0, info
      assert_equal lwork, get_int(work[0])
      assert_equal @rank_exp[rc], rank
      assert_narray @b_exp[rc], b, 10e-4
      assert_narray @s_exp[rc], s, 10e-4
    end
  end

end
