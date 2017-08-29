require "./spec_helper"

describe LA::Matrix do
  it "basic flags test" do
    a = LA::Mat.rand(10, 10)
    a.flags.should eq LA::MatrixFlags::None
    a += a.transpose
    a.detect(LA::MatrixFlags::Symmetric).should be_true
    a.transpose.flags.symmetric?.should be_true
    b = a + 5*LA::Mat.diag(*a.size, 1.5) / 2
    b.flags.should eq LA::MatrixFlags::Symmetric
    a.inv.flags.should eq LA::MatrixFlags::Symmetric
  end

  it "correct flags for diag, zeros, ones" do
    a = LA::Mat.diag(10, 4, 1.5)
    a.flags.lower_triangular?.should be_true
    a.flags.lower_triangular?.should be_true
    a.flags.symmetric?.should be_false
    a = LA::Mat.diag(4, 4, 1.5)
    a.flags.symmetric?.should be_true
    a = LA::Mat.ones(4, 10)
    a.flags.should eq LA::MatrixFlags::None
    a = LA::Mat.zeros(10, 4)
    a.flags.symmetric?.should be_false
  end

  it "flags for complex to double conversion" do
    a = LA::GMatComplex[{0, 6.i, 0, 0}, {0, 0, 6, 0}, {0, 0, 0, 6}, {0, 0, 0, 0}]
    a.detect
    a.flags.should eq LA::MatrixFlags::UpperTriangular
    a.to_real.flags.should eq LA::MatrixFlags::UpperTriangular
    a.to_imag.flags.should eq LA::MatrixFlags::UpperTriangular
    b = (1.i*a).expm
    b.to_real.flags.should eq LA::MatrixFlags::UpperTriangular
    b.to_imag.flags.should eq LA::MatrixFlags::UpperTriangular

    a = a + a.conjt
    a.flags.should eq LA::MatrixFlags::None
    a.detect
    a.flags.should eq LA::MatrixFlags::Hermitian
    (1.i*a).flags.should eq LA::MatrixFlags::None
    a.to_real.flags.should eq LA::MatrixFlags::Hermitian | LA::MatrixFlags::Symmetric
    a.to_imag.flags.should eq LA::MatrixFlags::None
    b = (1.i*a).expm
    b.flags.should eq LA::MatrixFlags::None
    b.detect
    b.flags.should eq LA::MatrixFlags::Orthogonal
    b.to_real.flags.should eq LA::MatrixFlags::None
    b.to_imag.flags.should eq LA::MatrixFlags::None
  end

  it "clean unused elements when inverting diagonal matrices" do
    u = LA::GMat[
      [1e+16, -1, 0],
      [0.0, 1e+16, 2],
      [0.0, 0.0, 1e+17],
    ]
    u.assume! LA::MatrixFlags::UpperTriangular | LA::MatrixFlags::LowerTriangular
    u.inv[0, 1].should be_close 0.0, 1e-16
  end
end
