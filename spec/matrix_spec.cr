require "./spec_helper"

include LAPACK
describe LAPACK::Matrix do
  it "can be created with given size" do
    m = Matrix(Float64).new(10, 15)
    m.raw.size.should eq 10*15
    m.raw[0].should eq 0
  end

  it "can be created from array with given dimension" do
    m = Matrix(Float64).new(3, 4, [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12])
    m.rows.should eq 3
    m[1, 0].should eq 5
  end

  it "can be created from array of arrays" do
    m = Matrix(Float64).new({ {1, 2, 3}, {4, 5, 6}, {7, 8, 9}, {10, 11, 12} })
    m.rows.should eq 4
    m.columns.should eq 3
    m[1, 0].should eq 4
    expect_raises(IndexError) { Matrix(Float64).new({ {1, 2, 3}, {4, 5} }) }
  end

  it "can be created from a dimensions and a block" do
    m = Matrix(Float32).new(3, 3) { |i, j| i*3 + j }
    m[1, 1].should eq 4
  end

  it "can't create matrix of unsupported type" do
    # commented as it causes a compile error, as it should
    # m = Matrix(Int32).new(3, 3)
  end

  it "access to members works and is zero-based" do
    m = Matrix(Float64).new(5, 4)
    m[4, 3] = 1.0
    m[4, 3].should eq 1.0
    expect_raises(IndexError) { m[5, 4] = 1.0 }
    expect_raises(IndexError) { m[0, 6] = 1.0 }
  end

  it "can multiply matrices" do
    m1 = Matrix(Float64).new(3, 4) { |i, j| i + j }
    m2 = Matrix(Float64).new(4, 2) { |i, j| i - j }
    m = m1*m2
    m.should eq Matrix(Float64).new(3, 2, [14.0, 8.0, 20.0, 10.0, 26.0, 12.0])
    expect_raises(ArgumentError) { Matrix(Float64).new(3, 4) * Matrix(Float64).new(3, 4) }
  end

  it "can do sum and scalar multiply" do
    m1 = Matrix(Float64).new(3, 4) { |i, j| i }
    m2 = Matrix(Float64).new(3, 4) { |i, j| j }
    m = m1 + m2*2
    m.should eq Matrix(Float64).new(3, 4) { |i, j| i + 2*j }
    expect_raises(ArgumentError) { Matrix(Float64).new(3, 4) + Matrix(Float64).new(4, 4) }
  end
end
