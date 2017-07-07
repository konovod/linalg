require "./spec_helper"

describe Linalg::Matrix do
  it "basic flags test" do
    a = Linalg::Mat.rand(10, 10)
    a.flags.should eq Linalg::MatrixFlags::None
    a += a.transpose
    a.detect(Linalg::MatrixFlags::Symmetric).should be_true
    a.transpose.flags.symmetric?.should be_true
    b = a + 5*Linalg::Mat.diag(*a.size, 1.5) / 2
    b.flags.should eq Linalg::MatrixFlags::Symmetric
    a.inv.flags.should eq Linalg::MatrixFlags::Symmetric
  end

  it "correct flags for diag, zeros, ones" do
    a = Linalg::Mat.diag(10, 4, 1.5)
    a.flags.lower_triangular?.should be_true
    a.flags.lower_triangular?.should be_true
    a.flags.symmetric?.should be_false
    a = Linalg::Mat.diag(4, 4, 1.5)
    a.flags.symmetric?.should be_true
    a = Linalg::Mat.ones(4, 10)
    a.flags.should eq Linalg::MatrixFlags::None
    a = Linalg::Mat.zeros(10, 4)
    a.flags.symmetric?.should be_false
  end

  it "flags for complex to double conversion" do
    a = Linalg::GMatComplex[{0, 6.i, 0, 0}, {0, 0, 6, 0}, {0, 0, 0, 6}, {0, 0, 0, 0}]
    a.detect
    a.flags.should eq Linalg::MatrixFlags::UpperTriangular
    a.to_real.flags.should eq Linalg::MatrixFlags::UpperTriangular
    a.to_imag.flags.should eq Linalg::MatrixFlags::UpperTriangular
    b = (1.i*a).expm
    b.to_real.flags.should eq Linalg::MatrixFlags::UpperTriangular
    b.to_imag.flags.should eq Linalg::MatrixFlags::UpperTriangular

    a = a + a.conjt
    a.flags.should eq Linalg::MatrixFlags::None
    a.detect
    a.flags.should eq Linalg::MatrixFlags::Hermitian
    (1.i*a).flags.should eq Linalg::MatrixFlags::None
    a.to_real.flags.should eq Linalg::MatrixFlags::Hermitian | Linalg::MatrixFlags::Symmetric
    a.to_imag.flags.should eq Linalg::MatrixFlags::None
    b = (1.i*a).expm
    b.flags.should eq Linalg::MatrixFlags::None
    b.detect
    b.flags.should eq Linalg::MatrixFlags::Orthogonal
    b.to_real.flags.should eq Linalg::MatrixFlags::None
    b.to_imag.flags.should eq Linalg::MatrixFlags::None
  end
end
