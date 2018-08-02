require "./spec_helper"
include LA

describe LA do
  it "high-level: calculate nonsymmetric eigenvalues" do
    a = GMat32[[-2, 4, 1], [2, -4, 1], [1, 1, 1]]
    vals = a.eigvals
    vals[0].should be_close -6, 1e-3
    vals[1].should be_close -1, 1e-3
    vals[2].should be_close 2, 1e-3
  end

  it "high-level: calculate nonsymmetric eigenvalues (complex result)" do
    a = GMat32[[3, -2], [4, -1]]
    vals = a.eigvals
    vals[0].should be_close 1 + 2.i, 1e-3
    vals[1].should be_close 1 - 2.i, 1e-3
  end
  it "high-level: calculate nonsymmetric eigenvalues (complex argument)" do
    a = GMatComplex[[3, -2], [4, -1]]
    vals = a.eigvals
    vals[0].should be_close 1 + 2.i, 1e-3
    vals[1].should be_close 1 - 2.i, 1e-3
  end

  it "high-level: calculate nonsymmetric eigenvectors" do
    a = GMat32[[-2, 4, 1], [2, -4, 1], [1, 1, 1]]
    vals, vectors = a.eigs
    raise "" if vals.is_a? Array(Complex)
    vals.each { |e| (a*vectors - vectors*e).det.should be_close 0, 1e-4 }

    a = GMatComplex[[-2, 4, 1], [2, -4, 1], [1, 1, 1]]
    vals, vectors = a.eigs(left: true)
    vals.each { |e| (vectors*a - vectors*e).det.should be_close 0, 1e-6 }
  end

  it "high-level: calculate symmetric eigenvalues (real argument)" do
    a = GMat32[
      [1, 2, 3],
      [2, 3, 4],
      [3, 4, 5],
    ]
    a.assume! MatrixFlags::Symmetric
    vals = a.eigvals
    vals[0].should be_close -0.623475, 1e-4
    vals[1].should be_close 0, 1e-4
    vals[2].should be_close 9.62348, 1e-4
  end
  it "high-level: calculate hermitian eigenvalues (complex argument)" do
    a = GMatComplex[
      [1, 2 - 1.i, 3 + 2.i],
      [2 + 1.i, 3, 4],
      [3 - 2.i, 4, 5],
    ]
    a.assume! MatrixFlags::Hermitian
    vals = a.eigvals
    vals[0].should be_close -2.26772, 1e-5
    vals[1].should be_close 1.48798, 1e-5
    vals[2].should be_close 9.77974, 1e-5
  end

  it "high-level: generalized nonsymmetric eigenvalues (real argument)" do
    a = GMat[
      [-2, 4, 1],
      [2, -4, 1],
      [1, 1, 1],
    ]
    b = GMat[
      [-2, 4, 0],
      [2, -3, 1],
      [2, 1, 1],
    ]
    alpha, beta, left, right = LA.eigs(a, b, need_left: true, need_right: true)
    # pp alpha, beta, left, right
    raise "" if alpha.is_a? Array(Complex)
    3.times do |i|
      (a * right.not_nil! - alpha[i] / beta[i] * b * right.not_nil!).det.should be_close(0, 1e-6)
      (left.not_nil!.transpose * a - alpha[i] / beta[i] * left.not_nil!.transpose * b).det.should be_close(0, 1e-6)
    end
  end

  it "high-level: generalized nonsymmetric eigenvalues (complex argument)" do
    a = GMatComplex[
      [-2, 4.i, 1],
      [3, -4.i, 1],
      [1, 0, 1 - 1.i],
    ]
    b = GMatComplex[
      [-2, 4, 0],
      [2, -3, -2.i],
      [1, 1 - 1.i, 1 + 1.i],
    ]
    alpha, beta, left, right = LA.eigs(a, b, need_left: true, need_right: true)
    3.times do |i|
      (a * right.not_nil! - alpha[i] / beta[i] * b * right.not_nil!).det.should be_close(0, 1e-6)
      (left.not_nil!.conjtranspose * a - alpha[i] / beta[i] * left.not_nil!.conjtranspose * b).det.should be_close(0, 1e-6)
    end
  end

  it "high-level: generalized symmetric eigenvalues (real argument)" do
    a = GMat[
      [-2, 4, 1],
      [4, -4, 3],
      [1, 3, 1],
    ]
    b = GMat[
      [2, -1, 0],
      [-1, 2, -1],
      [0, -1, 2],
    ]
    a.assume! MatrixFlags::Symmetric
    b.assume! MatrixFlags::PositiveDefinite
    alpha, beta, left, right = LA.eigs(a, b, need_left: true, need_right: true)
    raise "complex eigenvalue" if alpha.is_a? Array(Complex)
    3.times do |i|
      (a * right.not_nil! - alpha[i] / beta[i] * b * right.not_nil!).det.should be_close(0, 1e-6)
      (left.not_nil!.transpose * a - alpha[i] / beta[i] * left.not_nil!.transpose * b).det.should be_close(0, 1e-6)
    end
  end

  pending "high-level: generalized symmetric eigenvalues (complex argument)" do
    a = GMatComplex[
      [-2, 4.i, 1],
      [-4.i, -4, 0],
      [1, 0, 1],
    ]
    b = GMatComplex[
      [2, -1, 0],
      [-1, 2, -1],
      [0, -1, 2],
    ]
    a.assume! MatrixFlags::Hermitian
    b.assume! MatrixFlags::PositiveDefinite
    alpha, beta, left, right = LA.eigs(a, b, need_left: true, need_right: true)
    3.times do |i|
      (a * right.not_nil! - alpha[i] / beta[i] * b * right.not_nil!).det.should be_close(0, 1e-6)
      (left.not_nil!.conjtranspose * a - alpha[i] / beta[i] * left.not_nil!.conjtranspose * b).det.should be_close(0, 1e-6)
    end
  end
end
