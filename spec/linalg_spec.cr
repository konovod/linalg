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
    matrix1 = GMat.new([
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
    matrix1 = GMat.new([
      [1, 0, 1],
      [0, 4, 0],
      [0, 0, 1],
    ])
    (matrix1*matrix1.inv).should eq Mat.identity(3)
  end

  it "ultra-highlevel - calls functions on a virtual matrices" do
    m4 = Mat.ones(4, 4)*5 + Mat.rand(4, 4, Random.new(1))
    m3 = m4[0..2, 0..2]
    (m3*m3.inv).should be_close Mat.identity(3), 1e-9
  end

  it "raises LinAlgError on incorrect data" do
    matrix1 = GMat.new([
      [1, 0, 1],
      [0, 0, 0],
      [0, 0, 1],
    ])
    expect_raises(LinAlgError) do
      matrix1.inv!
    end
  end

  it "support all types" do
    matrix1 = GMat32.new([
      [1, 0, 1],
      [0, 4, 0],
      [0, 0, 1],
    ])
    (matrix1*matrix1.inv).should eq Mat32.identity(3)

    i = Complex.new(0, 1)
    matrix1 = GMatComplex.new({
      {1 + 1*i, 0, 1},
      {0, 4, 0},
      {0, 0, 1},
    })
    (matrix1*matrix1.inv).should eq MatComplex.identity(3)
  end

  it "high-level: solve linear equations" do
    a = GMat32.new(
      [[2, 4],
       [2, 8]]
    )
    b = GMat32.new([[2], [4]])
    Linalg.solve(a, b).should eq (a.inv * b)
  end

  it "high-level: calculate determinant" do
    a = GMat.new(
      [[1, 2, 3],
       [4, 5, 7],
       [-1, 1, -1]]
    )
    a.det.should eq 9
  end
  it "high-level: solve linear least square" do
    a = GMat32.new(
      [[1, 2, 0],
       [0, 4, 3]]
    )
    b = GMat32.new([[8], [18]])
    x = Linalg.lstsq(a, b)
    x_octave = GMat32.new(3, 1, [0.918032, 3.54098, 1.27869])
    x.should be_close(x_octave, 1e-3)
  end

  it "high-level: solve linear least square (complex)" do
    a = GMatComplex.new(
      [[1, 2, 0],
       [0, 4, 3]]
    )
    b = GMatComplex.new([[8], [18]])
    x = Linalg.lstsq(a, b)
    x_octave = GMatComplex.new(3, 1, [0.918032, 3.54098, 1.27869])
    x.should be_close(x_octave, 1e-3)
  end

  # sadly, spec is order-depentent
  it "high-level: calculate nonsymmetric eigenvalues" do
    a = GMat32.new([[-2, 4, 1], [2, -4, 1], [1, 1, 1]])
    vals = a.eigvals
    vals[0].should be_close -6, 1e-3
    vals[1].should be_close -1, 1e-3
    vals[2].should be_close 2, 1e-3
  end

  it "high-level: calculate nonsymmetric eigenvalues (complex result)" do
    a = GMat32.new([[3, -2], [4, -1]])
    vals = a.eigvals
    i = Complex.new(0, 1)
    vals[0].should be_close 1 + 2*i, 1e-3
    vals[1].should be_close 1 - 2*i, 1e-3
  end
  it "high-level: calculate nonsymmetric eigenvalues (complex argument)" do
    a = GMatComplex.new([[3, -2], [4, -1]])
    vals = a.eigvals
    i = Complex.new(0, 1)
    vals[0].should be_close 1 + 2*i, 1e-3
    vals[1].should be_close 1 - 2*i, 1e-3
  end

  it "high-level: calculate nonsymmetric eigenvectors" do
    a = GMat32.new([[-2, 4, 1], [2, -4, 1], [1, 1, 1]])
    vals, vectors = a.eigs
    raise "" if vals.is_a? Array(Complex)
    vals.each { |e| (a*vectors - vectors*e).det.should be_close 0, 1e-4 }

    a = GMatComplex.new([[-2, 4, 1], [2, -4, 1], [1, 1, 1]])
    vals, vectors = a.eigs(left: true)
    vals.each { |e| (vectors*a - vectors*e).det.should be_close 0, 1e-6 }
  end

  it "high-level: calculate singular value decomposition" do
    a = GMat.new([[1, 2, 3], [4, 5, 6]])
    u, s, vt = Linalg.svd(a)
    (u*Mat.diag(a.rows, a.columns, s)*vt).should be_close a, 1e-6

    ac = GMatComplex.new([[1, 2, 3], [4, 5, 6]])
    s1 = ac.svdvals(overwrite_a: true)
    s1[0].should be_close s[0], 1e-6
    s1[1].should be_close s[1], 1e-6
  end

  # TODO - proper spec
  it "high-level: balance matrix" do
    a = GMat32.new([[1, 2, 0], [9, 1, 0.01], [1, 2, 10*Math::PI]])
    b, s = a.balance(separate: true)
    s.should eq GMat32.new([[0.5, 1, 1]])
  end

  it "high-level: LU factorization" do
    a = GMat32.new([[1, 2, 3, 4], [5, 6, 7, 8], [1, 2, 10, 0], [50, 6, 7, 8], [1, 21, 1, 0]])
    p, l, u = a.lu
    (p*l*u).should be_close a, 1e-4
  end

  it "high-level: solve using LU" do
    a = GMat32.new(
      [[2, 4],
       [2, 8]]
    )
    lu = a.lu_factor
    b = GMat32.new([[2], [4]])
    lu.solve(b).should eq (a.inv * b)
  end
end
