$:.push File.dirname(__FILE__) + "/../.."
require "lapack_test"
require "numru/lapack"

class GelssTest < Test::Unit::TestCase
  include LapackTest

  def setup
    @a = NMatrix[[-0.09,  0.14, -0.46,  0.68,  1.29],
                 [-1.56,  0.20,  0.29,  1.09,  0.51],
                 [-1.48, -0.43,  0.89, -0.71, -0.96],
                 [-1.09,  0.84,  0.77,  2.11, -1.27],
                 [ 0.08,  0.55, -1.13,  0.14,  1.74],
                 [-1.59, -0.72,  1.06,  1.24,  0.34]].to_lm
    @b = NVector[[7.4, 4.2, -8.3, 1.8, 8.6, 2.1]]
    @rcond = 0.01

    @bout = NArray[[0.6344, 0.9699, -1.4403, 3.3678, 3.3992]]
    @s = NArray[3.9997, 2.9962, 2.0001, 0.9988, 0.0025]


    i = Complex::I
    @ac = NMatrix[[ 0.47-0.34*i, -0.40+0.54*i,  0.60+0.01*i,  0.80-1.02*i],
                  [-0.32-0.23*i, -0.05+0.20*i, -0.26-0.44*i, -0.43+0.17*i],
                  [ 0.35-0.60*i, -0.52-0.34*i,  0.87-0.11*i, -0.34-0.09*i],
                  [ 0.89+0.71*i, -0.45-0.45*i, -0.02-0.57*i,  1.14-0.78*i],
                  [-0.19+0.06*i,  0.11-0.85*i,  1.44+0.80*i,  0.07+1.14*i]].to_lm
    @bc = NVector[[-1.08-2.59*i, -2.61-1.49*i, 3.13-3.61*i, 7.33-8.01*i, 9.12+7.63*i]]

    @bcout = NArray[[1.1673-3.3222*i, 1.3480+5.5028*i, 4.1762+2.3434*i, 0.6465+0.0105*i]]
    @sc = NArray[2.9979, 1.9983, 1.0044, 0.0064]
  end

  %w(s d).each do |sd|
    method = "#{sd}gelss"
    define_method("test_#{method}") do
      s, rank, work, info, a, b = NumRu::Lapack.send(method, @a, @b, @rcond, :lwork => -1)
      assert_equal 0, info
      lwork = work[0].to_i
      [lwork, nil].each do |lw|
        s, rank, work, info, a, b = NumRu::Lapack.send(method, @a, @b, @rcond, :lwork => lw)
        assert_equal 0, info
        assert_equal lwork, work[0].to_i
        assert_narray @bout, b, 1e-4
        assert_narray @s, s, 1e-4
        assert 4, rank
      end
    end
  end

  %w(c z).each do |cz|
    method = "#{cz}gelss"
    define_method("test_#{method}") do
      s, rank, work, info, a, b = NumRu::Lapack.send(method, @ac, @bc, @rcond, :lwork => -1)
      assert_equal 0, info
      lwork = work[0].real.to_i
      [lwork, nil].each do |lw|
        s, rank, work, info, a, b = NumRu::Lapack.send(method, @ac, @bc, @rcond, :lwork => lw)
        assert_equal 0, info
        assert_equal lwork, work[0].real.to_i
        assert_narray @bcout, b, 1e-4
        assert_narray @sc, s, 1e-4
        assert 4, rank
      end
    end
  end

end
