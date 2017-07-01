require "./spec_helper"

j = Complex.new(0, 1)
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
end
