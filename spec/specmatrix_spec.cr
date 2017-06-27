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

  it "leslie matrix" do
    Mat.leslie([0.1, 2.0, 1.0, 0.1], [0.2, 0.8, 0.7]).should eq GMat.new [
      [0.1, 2, 1, 0.1],
      [0.2, 0, 0, 0],
      [0, 0.8, 0, 0],
      [0, 0, 0.7, 0],
    ]
  end

  it "companion matrix" do
    Mat.companion([1, -10, 31, -30]).should eq GMat.new [
      [10, -31, 30],
      [1, 0, 0],
      [0, 1, 0],
    ]
  end
  it "hadamard matrix" do
    expect_raises(ArgumentError) do
      Mat.hadamard(1000)
    end
    Mat.hadamard(4).should eq GMat.new [
      [1, 1, 1, 1],
      [1, -1, 1, -1],
      [1, 1, -1, -1],
      [1, -1, -1, 1],
    ]
  end

  it "hankel matrix" do
    Mat.hankel([1, 2, 3], [1, 4, 5, 6]).should eq GMat.new [
      [1, 2, 3, 4],
      [2, 3, 4, 5],
      [3, 4, 5, 6],
    ]
    MatComplex.hankel(GMatComplex.new [[1.0, 2 + 3*j, 4 - j]]).should eq GMatComplex.new [
      [1, 2 + 3*j, 4 - j],
      [2 + 3*j, 4 - j, 0],
      [4 - j, 0, 0],
    ]
  end

  it "helmert matrix" do
    h1 = GMat.new [
      [0.4472136, 0.4472136, 0.4472136, 0.4472136, 0.4472136],
      [0.70710678, -0.70710678, 0, 0, 0],
      [0.40824829, 0.40824829, -0.81649658, 0, 0],
      [0.28867513, 0.28867513, 0.28867513, -0.8660254, 0],
      [0.2236068, 0.2236068, 0.2236068, 0.2236068, -0.89442719],
    ]
    (Mat.helmert(5, true) - h1).abs.should be_close(0, 1e-6)
    (Mat.helmert(5, false) - h1[1..4, 0..4]).abs.should be_close(0, 1e-6)
  end

  it "hilbert matrix" do
    Mat.hilbert(3).should be_close(GMat.new([
      [1, 0.5, 0.33333333],
      [0.5, 0.33333333, 0.25],
      [0.33333333, 0.25, 0.2],
    ]), 1e-6)
  end

  it "dft matrix" do
    m = MatComplex.dft(8)
    x = GMatComplex.new [[1, 2, 3, 0, 3, 2, 1, 0]]
    (x*m).should almost_eq GMatComplex.new([
      [12, -2 - 2*j, -4*j, -2 + 2*j, 4, -2 - 2*j, 4*j, -2 + 2*j],
    ])
  end
end
