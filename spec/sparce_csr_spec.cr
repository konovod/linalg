require "./spec_helper"

include LA
include LA::Sparse
describe CSRMatrix do
  it "can be created empty" do
    ms = CSRMatrix(Float64).new(3, 4, 5)
    ms.should eq GMat[
      [0, 0, 0, 0],
      [0, 0, 0, 0],
      [0, 0, 0, 0],
    ]
  end

  it "can be constructed from arrays" do
    ms = CSRMatrix(Float64).new(4, 4, raw_rows: [0, 1, 2, 3, 4], raw_columns: [0, 1, 2, 1], raw_values: [5.0, 8.0, 3.0, 6.0])
    ms.should eq GMat[
      [5, 0, 0, 0],
      [0, 8, 0, 0],
      [0, 0, 3, 0],
      [0, 6, 0, 0],
    ]
    expect_raises(ArgumentError) { CSRMatrix(Float64).new(3, 4, raw_rows: [0, 1, 2, 3, 4], raw_columns: [0, 1, 2, 1], raw_values: [5.0, 8.0, 3.0, 6.0]) }
    expect_raises(ArgumentError) { CSRMatrix(Float64).new(4, 4, raw_rows: [0, 1, 2, 3, 4], raw_columns: [0, 1, 2, 1], raw_values: [5.0, 8.0, 3.0]) }
  end

  it "can be created from dense matrix" do
    m = GMat.eye(5)
    m[3, 1] = 2.0
    ms = CSRMatrix(Float64).new(m)
    ms.should eq m
    ms.nonzeros.should eq 6
  end

  it "can be created from COO matrix" do
    m = GMat.eye(5)
    m[3, 1] = 2.0
    m = COOMatrix(Float64).new(m)
    ms = CSRMatrix(Float64).new(m)
    ms.should eq m
    ms.nonzeros.should eq 6
  end

  it "can be cleared" do
    ms = CSRMatrix(Float64).new(4, 4, raw_rows: [0, 1, 2, 3, 4], raw_columns: [0, 1, 2, 1], raw_values: [5.0, 8.0, 3.0, 6.0])
    ms.clear
    ms.should eq GMat[
      [0, 0, 0, 0],
      [0, 0, 0, 0],
      [0, 0, 0, 0],
      [0, 0, 0, 0],
    ]
  end

  it "can be compared with matrix when empty" do
    m = GMat.eye(5)
    ms = CSRMatrix(Float64).new(*m.shape)
    ms.should_not eq m
    m.should_not eq ms
    ms.should eq GMat.new(5, 5)
    ms.clear
    ms.should_not eq m
    m.should_not eq ms
    ms.should eq GMat.new(5, 5)
  end

  it "elements can be accessed and changed" do
    m = GMat.eye(5)
    m[3, 1] = 2.0
    ms = CSRMatrix(Float64).new(m)
    ms[3, 1].should eq 2
    ms[3, 2].should eq 0
    ms[3, 1] = 0.0
    ms[3, 1].should eq 0
    expect_raises(Exception) { ms[4, 2] = 0 }
  end

  it "can be created using `diag`" do
    ms = CSRMatrix(Float64).eye(3)
    ms.should eq GMat.new [
      [1, 0, 0],
      [0, 1, 0],
      [0, 0, 1],
    ]
    ms = CSRMatrix(Float64).diag(4, 5, [1, 2, 3])
    ms.should eq GMat.new [
      [1, 0, 0, 0, 0],
      [0, 2, 0, 0, 0],
      [0, 0, 3, 0, 0],
      [0, 0, 0, 0, 0],
    ]
    ms.nonzeros.should eq 3
    expect_raises(Exception) { CSRMatrix(Float64).diag(4, 5, [1, 2, 3, 4, 5]) }
  end
end
