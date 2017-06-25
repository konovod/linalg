require "./spec_helper"

include Linalg
j = Complex.new(0, 1)
describe "constructing of special matrices" do
  it "block_diag matrix" do
    a = GMat.new [[1, 0],
                  [0, 1]]
    b = GMat.new [[3, 4, 5],
                  [6, 7, 8]]
    c = GMat.new [[7]]
    Mat.block_diag(a, b, c).should eq GMat.new [
      [1, 0, 0, 0, 0, 0],
      [0, 1, 0, 0, 0, 0],
      [0, 0, 3, 4, 5, 0],
      [0, 0, 6, 7, 8, 0],
      [0, 0, 0, 0, 0, 7],
    ]
  end

  it "toeplitz matrix" do
    Mat.toeplitz([1, 2, 3], [1, 4, 5, 6]).should eq GMat.new [
      [1, 4, 5, 6],
      [2, 1, 4, 5],
      [3, 2, 1, 4],
    ]
    MatComplex.toeplitz(GMatComplex.new [[1.0, 2 + 3*j, 4 - j]]).should eq GMatComplex.new [
      [1, 2 - 3*j, 4 + j],
      [2 + 3*j, 1, 2 - 3*j],
      [4 - j, 2 + 3*j, 1],
    ]
  end

  it "circulant matrix" do
    Mat.circulant([1, 2, 3]).should eq GMat.new [
      [1, 3, 2],
      [2, 1, 3],
      [3, 2, 1],
    ]
  end
end
