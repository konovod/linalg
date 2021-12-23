require "./spec_helper"

include LA
include LA::Sparse
describe COOMatrix do
  it "can be created from dense matrix" do
    m = GMat.eye(5)
    m[3, 1] = 2.0
    ms = COOMatrix(Float64).new(m)
    ms.should eq m
  end

  it "can be constructed from arrays" do
    ms = COOMatrix(Float64).new(3, 4, [0, 1, 2], [1, 1, 2], [5.0, 10.0, 15.0])
    ms.should eq GMat[
      [0, 5, 0, 0],
      [0, 10, 0, 0],
      [0, 0, 15, 0],
    ]
    expect_raises(ArgumentError) { COOMatrix(Float64).new(3, 3, [0, 1, 3], [1, 1, 2], [5.0, 10.0, 15.0]) }
    expect_raises(ArgumentError) { COOMatrix(Float64).new(3, 3, [0, 1, 2], [1, 1], [5.0, 10.0, 15.0]) }
  end

  it "can be created from another COO matrix" do
    m = GMat.eye(5)
    ms1 = COOMatrix(Float32).new(m)
    ms1[3, 1] = 1.0
    ms2 = COOMatrix(Float32).new(ms1)
    ms2.should eq ms1
    ms1[3, 1] = 2.0
    ms2.should_not eq ms1

    ms3 = COOMatrix(Float64).new(ms1)
    ms3.should eq ms1
  end

  it "can be created empty" do
    m = COOMatrix(Float64).new(2, 5)
    m[0, 0] = 1.0
    m[1, 2] = 5.0
    m.should eq GMat[
      [1, 0, 0, 0, 0],
      [0, 0, 5, 0, 0],
    ]
  end
end
