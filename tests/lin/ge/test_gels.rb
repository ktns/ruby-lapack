require "test/unit"
require "numru/lapack"

class GelsTest < Test::Unit::TestCase

  def setup
    @a = NMatrix[[-0.57, -1.28, -0.39,  0.25],
                 [-1.93,  1.08, -0.31, -2.14],
                 [ 2.30,  0.24,  0.40, -0.35],
                 [-1.93,  0.64, -0.66,  0.08],
                 [ 0.15,  0.30,  0.15, -2.13],
                 [-0.02,  1.03, -1.43,  0.50]].to_lm
    @b = NVector[[-2.67, -0.55, 3.34, -0.77, 0.48, 4.10]]

    @bout = NArray[[1.5339, 1.8707, -1.5241, 0.0392]]


    i = Complex::I
    @ac = NMatrix[[ 0.96-0.81*i, -0.03+0.96*i, -0.91+2.06*i, -0.05+0.41*i],
                  [-0.98+1.98*i, -1.20+0.19*i, -0.66+0.42*i, -0.81+0.56*i],
                  [ 0.62-0.46*i,  1.01+0.02*i,  0.63-0.17*i, -1.11+0.60*i],
                  [-0.37+0.38*i,  0.19-0.54*i, -0.98-0.36*i,  0.22-0.20*i],
                  [ 0.83+0.51*i,  0.20+0.01*i, -0.17-0.46*i,  1.47+1.59*i],
                  [ 1.08-0.28*i,  0.20-0.12*i, -0.07+1.23*i,  0.26+0.26*i]].to_lm
    @bc = NVector[[-2.09+1.93*i, 3.34-3.53*i, -4.94-2.04*i, 0.17+4.23*i, -5.19+3.63*i, 0.98+2.53*i]]

    @bcout = NArray[[-0.5044-1.2179*i, -2.4281+2.8574*i, 1.4872-2.1955*i, 0.4537+2.6904*i]]
  end

  %w(s d).each do |sd|
    method = "#{sd}gels"
    define_method("test_#{method}") do
      work, info, a, b = NumRu::Lapack.send(method, "N", @a, @b)
      assert_equal 0, info
      lwork = work[0].to_i
      work, info, a, b = NumRu::Lapack.send(method, "N", @a, @b, :lwork => -1)
      assert_equal 0, info
      assert_equal lwork, work[0].to_i
      [lwork, nil].each do |lw|
        work, info, a, b = NumRu::Lapack.send(method, "N", @a, @b, :lwork => lw)
        assert_equal 0, info
        assert_equal 12, work[0]
        assert_in_delta 0.0, (@bout - b).abs.max, 1.0e-4
      end
    end
  end

  %w(c z).each do |cz|
    method = "#{cz}gels"
    define_method("test_#{method}") do
      work, info, a, b = NumRu::Lapack.send(method, "N", @ac, @bc)
      assert_equal 0, info
      lwork = work[0].real.to_i
      work, info, a, b = NumRu::Lapack.send(method, "N", @ac, @bc, :lwork => -1)
      assert_equal 0, info
      assert_equal lwork, work[0].real.to_i
      work, info, a, b = NumRu::Lapack.send(method, "N", @ac, @bc, :lwork => lwork)
      assert_equal 0, info
      assert_equal 12, work[0]
      assert_in_delta 0.0, (@bcout - b).abs.max, 1.0e-4
    end
  end

end
