require "./spec_helper"

describe Linalg do
  it "calls LAPACK functions directly" do
    m = [
      1.0, 0.0, 1.0,
      0.0, 4.0, 0.0,
      0.0, 0.0, 1.0,
    ]
    ipiv = Slice(Int32).new(3)
    LibLAPACKE.dgetrf(LibLAPACKE::ROW_MAJOR, 3, 3, m, 3, ipiv).should eq 0
    LibLAPACKE.dgetri(LibLAPACKE::ROW_MAJOR, 3, m, 3, ipiv).should eq 0
    m.should eq [1.0, 0.0, -1.0, 0.0, 0.25, 0.0, 0.0, 0.0, 1.0]
  end

  it "calls functions using matrix class" do
    matrix1 = Mat.new([
      [1, 0, 1],
      [0, 4, 0],
      [0, 0, 1],
    ])
    matrix2 = matrix1*1
    ipiv = Slice(Int32).new(3)
    LibLAPACKE.dgetrf(LibLAPACKE::ROW_MAJOR, 3, 3, matrix2, 3, ipiv).should eq 0
    LibLAPACKE.dgetri(LibLAPACKE::ROW_MAJOR, 3, matrix2, 3, ipiv).should eq 0
    (matrix1*matrix2).should eq Mat.identity(3)
  end

  it "calls functions using high level wrapper" do
    matrix1 = Mat.new([
      [1, 0, 1],
      [0, 4, 0],
      [0, 0, 1],
    ])
    (matrix1*matrix1.inv).should eq Mat.identity(3)
  end

  it "support all types" do
    matrix1 = Mat32.new([
      [1, 0, 1],
      [0, 4, 0],
      [0, 0, 1],
    ])
    (matrix1*matrix1.inv).should eq Mat32.identity(3)

    i = Complex.new(0, 1)
    matrix1 = MatComplex.new({
      {1 + 1*i, 0, 1},
      {0, 4, 0},
      {0, 0, 1},
    })
    (matrix1*matrix1.inv).should eq MatComplex.identity(3)
  end

  it "high-level: solve linear equations" do
    a = Mat32.new(
      [[2, 4],
       [2, 8]]
    )
    b = Mat32.new([[2], [4]])
    Linalg.solve(a, b).should eq (a.inv * b)
  end

  it "high-level: calculate determinant" do
    a = Mat.new(
      [[1, 2, 3],
       [4, 5, 7],
       [-1, 1, -1]]
    )
    a.det.should eq 9
  end
  it "high-level: solve linear least square" do
    a = Mat32.new(
      [[1, 2, 0],
       [0, 4, 3]]
    )
    b = Mat32.new([[8], [18]])
    x = Linalg.lstsq(a, b)
    x_octave = Mat32.new(3, 1, [0.918032, 3.54098, 1.27869])
    x.should be_close(x_octave, 1e-3)
  end

  it "high-level: solve linear least square (complex)" do
    a = MatComplex.new(
      [[1, 2, 0],
       [0, 4, 3]]
    )
    b = MatComplex.new([[8], [18]])
    x = Linalg.lstsq(a, b)
    x_octave = MatComplex.new(3, 1, [0.918032, 3.54098, 1.27869])
    x.should be_close(x_octave, 1e-3)
  end

  # sadly, spec is order-depentent
  it "high-level: calculate nonsymmetric eigenvalues" do
    a = Mat32.new([[-2, 4, 1], [2, -4, 1], [1, 1, 1]])
    vals = a.eigvals
    vals[0].should be_close -6, 1e-3
    vals[1].should be_close -1, 1e-3
    vals[2].should be_close 2, 1e-3
  end

  it "high-level: calculate nonsymmetric eigenvalues (complex result)" do
    a = Mat32.new([[3, -2], [4, -1]])
    vals = a.eigvals
    i = Complex.new(0, 1)
    vals[0].should be_close 1 + 2*i, 1e-3
    vals[1].should be_close 1 - 2*i, 1e-3
  end
  it "high-level: calculate nonsymmetric eigenvalues (complex argument)" do
    a = MatComplex.new([[3, -2], [4, -1]])
    vals = a.eigvals
    i = Complex.new(0, 1)
    vals[0].should be_close 1 + 2*i, 1e-3
    vals[1].should be_close 1 - 2*i, 1e-3
  end

  it "high-level: calculate nonsymmetric eigenvectors" do
    a = Mat32.new([[-2, 4, 1], [2, -4, 1], [1, 1, 1]])
    vals, vectors = a.eigs
    raise "" if vals.is_a? Array(Complex)
    vals.each { |e| (a*vectors - vectors*e).det.should be_close 0, 1e-4 }

    a = MatComplex.new([[-2, 4, 1], [2, -4, 1], [1, 1, 1]])
    vals, vectors = a.eigs(left: true)
    vals.each { |e| (vectors*a - vectors*e).det.should be_close 0, 1e-6 }
  end

  it "high-level: calculate singular value decomposition" do
    a = Mat.new([[1, 2, 3], [4, 5, 6]])
    u, s, vt = Linalg.svd(a)
    (u*Mat.diag(a.rows, a.columns, s)*vt).should be_close a, 1e-6

    ac = MatComplex.new([[1, 2, 3], [4, 5, 6]])
    s1 = ac.svdvals(overwrite_a: true)
    s1[0].should be_close s[0], 1e-6
    s1[1].should be_close s[1], 1e-6
  end

  # TODO - proper spec
  it "high-level: balance matrix" do
    a = Mat32.new([[1, 2, 0], [9, 1, 0.01], [1, 2, 10*Math::PI]])
    b, s = a.balance(separate: true)
    s.should eq Mat32.new([[0.5, 1, 1]])
  end

  it "high-level: LU factorization" do
    a = Mat32.new([[1, 2, 3, 4], [5, 6, 7, 8], [1, 2, 10, 0], [50, 6, 7, 8], [1, 21, 1, 0]])
    p, l, u = a.lu
    (p*l*u).should be_close a, 1e-4
  end

  it "high-level: solve using LU" do
    a = Mat32.new(
      [[2, 4],
       [2, 8]]
    )
    lu = a.lu_factor
    b = Mat32.new([[2], [4]])
    lu.solve(b).should eq (a.inv * b)
  end
end
