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
end
