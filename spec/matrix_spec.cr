require "./spec_helper"

include LAPACK
describe LAPACK::Matrix do
  it "can be created with given size" do
    m = Matrix(Float64).new(10, 15)
    m.raw.size.should eq 10*15
    m.raw[0].should eq 0
  end

  it "can be created from array with given dimension" do
    m = Matrix(Float64).new(4, 3, [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12])
    m.n.should eq 4
    m[1, 0].should eq 5
  end

  it "can be created from array of arrays" do
    m = Matrix(Float64).new({ {1, 2, 3}, {4, 5, 6}, {7, 8, 9}, {10, 11, 12} })
    m.n.should eq 3
    m.m.should eq 4
    m[1, 0].should eq 4
  end

  it "can be created from a dimensions and a block" do
    m = Matrix(Float32).new(3, 3) { |i, j| i*3 + j }
    m[1, 1].should eq 4
  end
end
