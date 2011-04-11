require "test/unit"
require "numru/lapack"

class GelsdTest < Test::Unit::TestCase

  def setup
    @a = NMatrix[[-0.09, -1.56, -1.48, -1.09,  0.08, -1.59],
                 [ 0.14,  0.20, -0.43,  0.84,  0.55, -0.72],
                 [-0.46,  0.29,  0.89,  0.77, -1.13,  1.06],
                 [ 0.68,  1.09, -0.71,  2.11,  0.14,  1.24],
                 [ 1.29,  0.51, -0.96, -1.27,  1.74,  0.34]].to_lm
    @b = NVector[[7.4, 4.3, -8.1, 1.8, 8.7]]
    @rcond = 0.01

    @bout = NArray[[1.5938, -0.1180, -3.1501, 0.1554, 2.5529, -1.6730]]
    @s = NArray[3.9997, 2.9962, 2.0001, 0.9988, 0.0025]


    i = Complex::I
    @ac = NMatrix[[ 0.47-0.34*i, -0.32-0.23*i,  0.35-0.60*i,  0.89+0.71*i, -0.19+0.06*i],
                  [-0.40+0.54*i, -0.05+0.20*i, -0.52-0.34*i, -0.45-0.45*i,  0.11-0.85*i],
                  [ 0.60+0.01*i, -0.26-0.44*i,  0.87-0.11*i, -0.02-0.57*i,  1.44+0.80*i],
                  [ 0.80-1.02*i, -0.43+0.17*i, -0.34-0.09*i,  1.14-0.78*i,  0.07+1.14*i]].to_lm
    @bc = NVector[[2.15-0.20*i, -2.24+1.82*i, 4.45-4.28*i, 5.70-6.25*i]]

    @bcout = NArray[[3.9747-1.8377*i, -0.9186+0.8253*i, -0.3105+0.1477*i, 1.0050+0.8626*i, -0.2256-1.9425*i]]
    @sc = NArray[2.9979, 1.9983, 1.0044, 0.0064]
  end

  %w(s d).each do |sd|
    method = "#{sd}gelsd"
    define_method("test_#{method}") do
      s, rank, work, info, b = NumRu::Lapack.send(method, @a, @b, @rcond, -1)
      assert_equal 0, info
      lwork = work[0].to_i
      s, rank, work, info, b = NumRu::Lapack.send(method, @a, @b, @rcond, lwork)
      assert_equal 0, info
      assert_equal lwork, work[0].to_i
      assert_equal 4, rank
      assert_in_delta 0.0, (@bout-b).abs.max, 10e-4
      assert_in_delta 0.0, (@s-s).abs.max, 10e-4
    end
  end

  %w(c z).each do |cz|
    method = "#{cz}gelsd"
    define_method("test_#{method}") do
      s, rank, work, info, b = NumRu::Lapack.send(method, @ac, @bc, @rcond, -1)
      assert_equal 0, info
      lwork = work[0].real.to_i
      s, rank, work, info, b = NumRu::Lapack.send(method, @ac, @bc, @rcond, lwork)
      assert_equal 0, info
      assert_equal lwork, work[0].real.to_i
      assert_equal 3, rank
      assert_in_delta 0.0, (@bcout-b).abs.max, 10e-4
      assert_in_delta 0.0, (@sc-s).abs.max, 10e-4
    end
  end

end
