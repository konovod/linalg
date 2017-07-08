require "./spec_helper"

include Linalg
describe "constructing of special matrices" do
  it "block_diag matrix" do
    a = Mat.eye(2)
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
    MatComplex.toeplitz(GMatComplex.new [[1.0, 2 + 3.i, 4 - 1.i]]).should eq GMatComplex.new [
      [1, 2 - 3.i, 4 + 1.i],
      [2 + 3.i, 1, 2 - 3.i],
      [4 - 1.i, 2 + 3.i, 1],
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
    MatComplex.hankel(GMatComplex.new [[1.0, 2 + 3.i, 4 - 1.i]]).should eq GMatComplex.new [
      [1, 2 + 3.i, 4 - 1.i],
      [2 + 3.i, 4 - 1.i, 0],
      [4 - 1.i, 0, 0],
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
      [12, -2 - 2.i, -4.i, -2 + 2.i, 4, -2 - 2.i, 4.i, -2 + 2.i],
    ])
  end

  it "pascal matrix" do
    Mat.pascal(4).should eq GMat[
      [1, 1, 1, 1],
      [1, 2, 3, 4],
      [1, 3, 6, 10],
      [1, 4, 10, 20]]
    MatComplex.pascal(4, PascalKind::Lower).should eq GMatComplex[
      [1, 0, 0, 0],
      [1, 1, 0, 0],
      [1, 2, 1, 0],
      [1, 3, 3, 1]]
    Mat32.pascal(4, PascalKind::Upper).should eq Mat32.pascal(4, PascalKind::Lower).transpose
    Mat.pascal(50)[-1, -1].should be_close(2.547761225898085e28, 1e29)
  end

  it "invpascal matrix" do
    Mat.invpascal(10, PascalKind::Upper).should almost_eq Mat.pascal(10, PascalKind::Upper).inv
    Mat32.invpascal(10, PascalKind::Lower).should almost_eq Mat32.pascal(10, PascalKind::Lower).inv
    MatComplex.invpascal(10).should almost_eq MatComplex.pascal(10).inv
  end
end
