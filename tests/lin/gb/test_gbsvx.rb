$:.push File.dirname(__FILE__) + "/../.."
require "lapack_test"

class GbsvxTest < Test::Unit::TestCase
  include LapackTest

  def setup
    @ab = Hash.new
    @b = Hash.new
    @x_exp = Hash.new
    @ipiv_exp = Hash.new
    @rcond_exp = Hash.new

    @kl = 1
    @ku = 2
    @rpgf_exp = 1.0

    @ab[:r] = NMatrix[[-0.23,  2.54, -3.66,  0.00],
                      [-6.98,  2.46, -2.73, -2.13],
                      [ 0.00,  2.56,  2.46,  4.07],
                      [ 0.00,  0.00, -4.78, -3.82]].to_lb(@kl, @ku)
    @b[:r] = NVector[[  4.42,  27.13, -6.14,  10.50],
                     [-36.01, -31.67, -1.16, -25.82]]
    @x_exp[:r] = NArray[[-2.0,  3.0, 1.0, -4.0],
                        [ 1.0, -4.0, 7.0, -2.0]]
    @ipiv_exp[:r] = NArray[2, 3, 3, 4]
    @rcond_exp[:r] = 1.8e-2

    @ab[:c] = NMatrix[[-1.65+2.26*I, -2.05-0.85*I,  0.97-2.84*I, 0.00+0.00*I],
                      [ 0.00+6.30*I, -1.48-1.75*I, -3.99+4.01*I, 0.59-0.48*I],
                      [ 0.00+0.00*I, -0.77+2.83*I, -1.06+1.94*I, 3.33-1.04*I],
                      [ 0.00+0.00*I,  0.00+0.00*I,  4.48-1.09*I, -0.46-1.72*I]].to_lb(@kl, @ku)
    @b[:c] = NVector[[-1.06+21.50*I, -22.72-53.90*I, 28.24-38.60*I, -34.56+16.73*I],
                     [12.85+ 2.84*I, -70.22+21.57*I, -20.73- 1.23*I,  26.01+31.97*I]]
    @x_exp[:c] = NArray[[-3.0+2.0*I,  1.0-7.0*I, -5.0+4.0*I,  6.0-8.0*I],
                        [ 1.0+6.0*I, -7.0-4.0*I,  3.0+5.0*I, -8.0+2.0*I]]
    @ipiv_exp[:c] = NArray[2, 3, 3, 4]
    @rcond_exp[:c] = 9.6e-3
  end

  %w(s d c z).each do |x|
    method = "#{x}gbsvx"
    rc = LapackTest.get_rc(x)

    define_method("test_#{method}") do
      x, rcond, ferr, berr, work, info, ab, afb, ipiv, equed, r, c, b = NumRu::Lapack.send(method, "E", "N", @kl, @ku, @ab[rc], @b[rc])
      assert_equal(0, info)
      assert_narray @x_exp[rc], x
      assert_equal @ipiv_exp[rc], ipiv
      assert_in_delta @rcond_exp[rc], rcond, 1.0e-3
      assert_equal @rpgf_exp, work[0]
    end

  end

end
