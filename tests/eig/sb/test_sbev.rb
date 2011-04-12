$:.push File.dirname(__FILE__) + "/../.."
require "lapack_test"

class SbevTest < Test::Unit::TestCase
  include LapackTest

  def setup
    @kd = 2
    @ab =  NMatrix[[1.0, 2.0, 3.0, 0.0, 0.0],
                      [2.0, 2.0, 3.0, 4.0, 0.0],
                      [3.0, 3.0, 3.0, 4.0, 5.0],
                      [0.0, 4.0, 4.0, 4.0, 5.0],
                      [0.0, 0.0, 5.0, 5.0, 5.0]]
    @w_exp = NArray[-3.2474, -2.6633, 1.7511, 4.1599, 14.9997]
    @z_exp = NArray[[-0.0394, -0.5721,  0.4372,  0.4424, -0.5332],
                        [ 0.6238, -0.2575, -0.5900,  0.4308,  0.1039],
                        [ 0.5635, -0.3896,  0.4008, -0.5581,  0.2421],
                        [-0.5165, -0.5955, -0.1470,  0.0470,  0.5956],
                        [-0.1582, -0.3161, -0.5277, -0.5523, -0.5400] ]
  end

  %w(s d).each do |x|
    method = "#{x}sbev"

    %w(U L).each do |uplo|
      define_method("test_#{method}_uplo_#{uplo}") do
        w, z, info, ab = NumRu::Lapack.send(method, "V", uplo, @kd, @ab.to_lsb(uplo, @kd))
        assert_equal 0, info
        assert_narray @w_exp, w, 1.0e-4
        for n in 0...z.shape[1]
          z[true,n] *= -1 if @z_exp[0,n]*z[0,n] < 0
        end
        assert_narray @z_exp, z, 1.0e-4
      end
    end

  end

end
