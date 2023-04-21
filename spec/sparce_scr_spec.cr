require "./spec_helper"

include LA
include LA::Sparse
describe COOMatrix do
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
end
