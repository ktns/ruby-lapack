$:.push File.dirname(__FILE__) + "/../.."
require "lapack_test"

class GgevTest < Test::Unit::TestCase
  include LapackTest

  def setup
    @a = Hash.new
    @b = Hash.new
    @vr_exp = Hash.new

    @a[:r] = NMatrix[[3.9,  12.5, -34.5,  -0.5],
                     [4.3,  21.5, -47.5,   7.5],
                     [4.3,  21.5, -43.5,   3.5],
                     [4.4,  26.0, -46.0,   6.0]].to_lm
    @b[:r] = NMatrix[[1.0, 2.0, -3.0, 1.0],
                     [1.0, 3.0, -5.0, 4.0],
                     [1.0, 3.0, -4.0, 3.0],
                     [1.0, 3.0, -4.0, 4.0]].to_lm

    @evr_exp = NArray[2.0, 3.0,  3.0, 4.0]
    @evi_exp = NArray[0.0, 4.0, -4.0, 0.0]
    @vr_exp[:r] = NArray[[ 1.0000e-0,  5.7143e-3,  6.2857e-2,  6.2857e-2],
                         [-4.3979e-1, -8.7958e-2, -1.4241e-1, -1.4241e-1],
                         [-5.6021e-1, -1.1204e-1,  3.1418e-3,  3.1418e-3],
                         [-1.0000e+0, -1.1111e-2,  3.3333e-2, -1.5556e-1]]

    @a[:c] = NMatrix[[-21.10-22.50*I, 53.50-50.50*I, -34.50+127.50*I,   7.50 +0.50*I],
                     [ -0.46 -7.78*I, -3.50-37.50*I, -15.50 +58.50*I, -10.50 -1.50*I],
                     [  4.30 -5.50*I, 39.70-17.10*I, -68.50 +12.50*I,  -7.50 -3.50*I],
                     [  5.50 +4.40*I, 14.40+43.30*I, -32.50 -46.00*I, -19.00-32.50*I]].to_lm
    @b[:c] = NMatrix[[1.00-5.00*I,  1.60+1.20*I, -3.00+0.00*I,  0.00-1.00*I],
                     [0.80-0.60*I,  3.00-5.00*I, -4.00+3.00*I, -2.40-3.20*I],
                     [1.00+0.00*I,  2.40+1.80*I, -4.00-5.00*I,  0.00-3.00*I],
                     [0.00+1.00*I, -1.80+2.40*I,  0.00-4.00*I,  4.00-5.00*I]].to_lm

    @ev_exp = NArray[3.0-9.0*I, 2.0-5.0*I, 3.0-1.0*I, 4.0-5.0*I]
    @vr_exp[:c] = NArray[[-8.2377e-1-1.7623e-1*I, -1.5295e-1+7.0655e-2*I, -7.0655e-2-1.5295e-1*I,  1.5295e-1-7.0655e-2*I],
                         [ 6.3974e-1+3.6026e-1*I,  4.1597e-3-5.4650e-4*I,  4.0212e-2+2.2645e-2*I, -2.2645e-2+4.0212e-2*I],
                         [ 9.7754e-1+2.2465e-2*I,  1.5910e-1-1.1371e-1*I,  1.2090e-1-1.5371e-1*I,  1.5371e-1+1.2090e-1*I],
                         [-9.0623e-1+9.3766e-2*I, -7.4303e-3+6.8750e-3*I,  3.0208e-2-3.1255e-3*I, -1.4586e-2-1.4097e-1*I]]


  end

  %w(s d).each do |x|
    method = "#{x}ggev"
    rc = :r

    define_method("test_#{method}") do
      alphar, alphai, beta, vl, vr, work, info, a, b = NumRu::Lapack.send(method, "N", "V", @a[rc], @b[rc])
      assert_equal 0, info
      assert_narray @evr_exp, alphar/beta, 1.0e-4
      assert_narray @evi_exp, alphai/beta, 1.0e-4
      vr.shape[1].times do |i|
        vr[true,i] *= -1 if comp_sign(@vr_exp[rc][0,i], vr[0,i])
      end
      assert_narray @vr_exp[rc], vr, 1.0e-4
    end

    define_method("test_#{method}_inquiring_lwork") do
      alphar, alphai, beta, vl, vr, work, info, a, b = NumRu::Lapack.send(method, "N", "V", @a[rc], @b[rc], :lwork => -1)
      assert_equal 0, info
      lwork = get_int(work[0])
      alphar, alphai, beta, vl, vr, work, info, a, b = NumRu::Lapack.send(method, "N", "V", @a[rc], @b[rc], :lwork => lwork)
      assert_equal 0, info
      assert_narray @evr_exp, alphar/beta, 1.0e-4
      assert_narray @evi_exp, alphai/beta, 1.0e-4
      vr.shape[1].times do |i|
        vr[true,i] *= -1 if comp_sign(@vr_exp[rc][0,i], vr[0,i])
      end
      assert_narray @vr_exp[rc], vr, 1.0e-4
    end

    define_method("test_#{method}_inquiring_lwork_oldargstyle") do
      alphar, alphai, beta, vl, vr, work, info, a, b = NumRu::Lapack.send(method, "N", "V", @a[rc], @b[rc], :lwork => -1)
      assert_equal 0, info
      lwork = get_int(work[0])
      alphar, alphai, beta, vl, vr, work, info, a, b = NumRu::Lapack.send(method, "N", "V", @a[rc], @b[rc], -1)
      assert_equal 0, info
      assert_equal lwork, get_int(work[0])
    end

  end

  %w(c z).each do |x|
    method = "#{x}ggev"
    rc = :c

    define_method("test_#{method}") do
      alpha, beta, vl, vr, work, rwork, info, a, b = NumRu::Lapack.send(method, "N", "V", @a[rc], @b[rc])
      assert_equal 0, info
      assert_narray @ev_exp, alpha/beta, 1.0e-4
      vr.shape[1].times do |i|
        vr[true,i] *= -1 if comp_sign(@vr_exp[rc][0,i], vr[0,i])
      end
      assert_narray @vr_exp[rc], vr, 2.0e-2
    end

    define_method("test_#{method}_inquiring_lwork") do
      alpha, beta, vl, vr, work, rwork, info, a, b = NumRu::Lapack.send(method, "N", "V", @a[rc], @b[rc], :lwork => -1)
      assert_equal 0, info
      lwork = get_int(work[0])
      alpha, beta, vl, vr, work, rwork, info, a, b = NumRu::Lapack.send(method, "N", "V", @a[rc], @b[rc], :lwork => lwork)
      assert_equal 0, info
      assert_narray @ev_exp, alpha/beta, 1.0e-4
      vr.shape[1].times do |i|
        vr[true,i] *= -1 if comp_sign(@vr_exp[rc][0,i], vr[0,i])
      end
      assert_narray @vr_exp[rc], vr, 2.0e-2
    end

    define_method("test_#{method}_inquiring_lwork_oldargstyle") do
      alpha, beta, vl, vr, work, rwork, info, a, b = NumRu::Lapack.send(method, "N", "V", @a[rc], @b[rc], :lwork => -1)
      assert_equal 0, info
      lwork = get_int(work[0])
      alpha, beta, vl, vr, work, rwork, info, a, b = NumRu::Lapack.send(method, "N", "V", @a[rc], @b[rc], -1)
      assert_equal 0, info
      assert_equal lwork, get_int(work[0])
    end

  end

end
