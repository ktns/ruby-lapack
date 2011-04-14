$:.push File.dirname(__FILE__) + "/../.."
require "lapack_test"

class GesvxTest < Test::Unit::TestCase
  include LapackTest

  def setup
    @a = Hash.new
    @b = Hash.new
    @x_exp = Hash.new
    @ipiv_exp = Hash.new
    @rcond_exp = Hash.new
    @rpgf_exp = Hash.new

    @a[:r] = NMatrix[[   1.80,    2.88,  2.05,   -0.89],
                     [ 525.00, -295.00,-95.00, -380.00],
                     [   1.58,   -2.69, -2.90,   -1.04],
                     [  -1.11,   -0.66, -0.59,    0.80]].to_lm
    @b[:r] = NVector[[ 9.52, 2435.00,   0.77, -6.22],
                     [18.47,  225.00, -13.28, -6.21]]
    @x_exp[:r] = NArray[[1.0, -1.0, 3.0, -5.0],
                        [3.0,  2.0, 4.0,  1.0]]
    @rcond_exp[:r] = 1.8e-2
    @rpgf_exp[:r] = 7.4e-1

    @a[:c] = NMatrix[[-1.34 +2.55*I,  0.28+3.17*I, -6.39 -2.20*I,  0.72 -0.92*I],
                     [-1.70-14.10*I, 33.10-1.50*I, -1.50+13.40*I, 12.90+13.80*I],
                     [-3.29 -2.39*I, -1.91+4.42*I, -0.14 -1.35*I,  1.72 +1.35*I],
                     [ 2.41 +0.39*I, -0.56+1.47*I, -0.83 -0.69*I, -1.96 +0.67*I]].to_lm
    @b[:c] = NVector[[26.26+51.78*I,  64.30-86.80*I, -5.75+25.31*I,  1.16+2.57*I],
                     [31.32 -6.70*I, 158.60-14.20*I, -2.15+30.19*I, -2.56+7.55*I]]
    @x_exp[:c] = NArray[[ 1.0+1.0*I, 2.0-3.0*I, -4.0-5.0*I, 0.0+6.0*I],
                        [-1.0-2.0*I, 5.0+1.0*I, -3.0+4.0*I, 2.0-3.0*I]]
    @rcond_exp[:c] = 1.0e-2
    @rpgf_exp[:c] = 8.3e-1
  end

  %w(s d c z).each do |x|
    method = "#{x}gesvx"
    rc = LapackTest.get_rc(x)

    define_method("test_#{method}") do 
      x, rcond, ferr, berr, work, info, a, af, ipiv, equed, r, c, b = NumRu::Lapack.send(method, "E", "N", @a[rc], @b[rc])
      assert_equal 0, info
      assert_narray @x_exp[rc], x, 1.0e-4
      assert_in_delta @rcond_exp[rc], rcond, 1e-3
      assert_in_delta @rpgf_exp[rc], work[0], 1e-2
    end

  end

end
