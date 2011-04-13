$:.push File.dirname(__FILE__) + "/../.."
require "lapack_test"

class GelssTest < Test::Unit::TestCase
  include LapackTest

  def setup
    @a = Hash.new
    @b = Hash.new
    @b_exp = Hash.new
    @s_exp = Hash.new
    @rank_exp = Hash.new

    @a[:r] = NMatrix[[-0.09,  0.14, -0.46,  0.68,  1.29],
                     [-1.56,  0.20,  0.29,  1.09,  0.51],
                     [-1.48, -0.43,  0.89, -0.71, -0.96],
                     [-1.09,  0.84,  0.77,  2.11, -1.27],
                     [ 0.08,  0.55, -1.13,  0.14,  1.74],
                     [-1.59, -0.72,  1.06,  1.24,  0.34]].to_lm
    @b[:r] = NVector[[7.4, 4.2, -8.3, 1.8, 8.6, 2.1]]
    @rank_exp[:r] = 4
    @b_exp[:r] = NArray[[0.6344, 0.9699, -1.4403, 3.3678, 3.3992]]
    @s_exp[:r] = NArray[3.9997, 2.9962, 2.0001, 0.9988, 0.0025]

    @a[:c] = NMatrix[[ 0.47-0.34*I, -0.40+0.54*I,  0.60+0.01*I,  0.80-1.02*I],
                     [-0.32-0.23*I, -0.05+0.20*I, -0.26-0.44*I, -0.43+0.17*I],
                     [ 0.35-0.60*I, -0.52-0.34*I,  0.87-0.11*I, -0.34-0.09*I],
                     [ 0.89+0.71*I, -0.45-0.45*I, -0.02-0.57*I,  1.14-0.78*I],
                     [-0.19+0.06*I,  0.11-0.85*I,  1.44+0.80*I,  0.07+1.14*I]].to_lm
    @b[:c] = NVector[[-1.08-2.59*I, -2.61-1.49*I, 3.13-3.61*I, 7.33-8.01*I, 9.12+7.63*I]]
    @b_exp[:c] = NArray[[1.1673-3.3222*I, 1.3480+5.5028*I, 4.1762+2.3434*I, 0.6465+0.0105*I]]
    @s_exp[:c] = NArray[2.9979, 1.9983, 1.0044, 0.0064]
    @rank_exp[:c] = 3

    @rcond = 0.01
  end

  %w(s d c z).each do |x|
    method = "#{x}gelss"
    rc = LapackTest.get_rc(x)

    define_method("test_#{method}") do
      s, rank, work, info, a, b = NumRu::Lapack.send(method, @a[rc], @b[rc], @rcond)
      assert_equal 0, info
      assert_narray @b_exp[rc], b, 1e-4
      assert_narray @s_exp[rc], s, 1e-4
      assert @rank_exp[rc], rank
    end

    define_method("test_#{method}_inquiring_lwork") do
      s, rank, work, info, a, b = NumRu::Lapack.send(method, @a[rc], @b[rc], @rcond, :lwork => -1)
      assert_equal 0, info
      lwork = get_int(work[0])
      s, rank, work, info, a, b = NumRu::Lapack.send(method, @a[rc], @b[rc], @rcond, :lwork => lwork)
      assert_equal 0, info
      assert_equal lwork, get_int(work[0])
      assert_narray @b_exp[rc], b, 1e-4
      assert_narray @s_exp[rc], s, 1e-4
      assert @rank_exp[rc], rank
    end

    define_method("test_#{method}_inquiring_lwork_oldargstyle") do
      s, rank, work, info, a, b = NumRu::Lapack.send(method, @a[rc], @b[rc], @rcond, :lwork => -1)
      assert_equal 0, info
      lwork = get_int(work[0])
      s, rank, work, info, a, b = NumRu::Lapack.send(method, @a[rc], @b[rc], @rcond, -1)
      assert_equal 0, info
      assert_equal lwork, get_int(work[0])
    end

  end

end
