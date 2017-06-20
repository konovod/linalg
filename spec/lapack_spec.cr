require "./spec_helper"

describe LAPACK do
  # TODO: Write tests

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
    matrix1 = Matrix(Float64).new([
      [1, 0, 1],
      [0, 4, 0],
      [0, 0, 1],
    ])
    matrix2 = matrix1*1
    ipiv = Slice(Int32).new(3)
    LibLAPACKE.dgetrf(LibLAPACKE::ROW_MAJOR, 3, 3, matrix2, 3, ipiv).should eq 0
    LibLAPACKE.dgetri(LibLAPACKE::ROW_MAJOR, 3, matrix2, 3, ipiv).should eq 0
    (matrix1*matrix2).should eq Matrix(Float64).identity(3)
  end

  it "calls functions using high level wrapper" do
    matrix1 = Matrix(Float64).new([
      [1, 0, 1],
      [0, 4, 0],
      [0, 0, 1],
    ])
    (matrix1*matrix1.inv).should eq Matrix(Float64).identity(3)
  end

  it "support all types" do
    matrix1 = Matrix(Float32).new([
      [1, 0, 1],
      [0, 4, 0],
      [0, 0, 1],
    ])
    (matrix1*matrix1.inv).should eq Matrix(Float32).identity(3)

    i = Complex.new(0, 1)
    matrix1 = Matrix(Complex).new({
      {1 + 1*i, 0, 1},
      {0, 4, 0},
      {0, 0, 1},
    })
    (matrix1*matrix1.inv).should eq Matrix(Complex).identity(3)
  end

  # info = LibLAPACKE.sgeqrf(LibLAPACKE::ROW_MAJOR,
  #   3,
  #   3,
  #   matrix,
  #   3,
  #   out tau)
  # pp info, tau, matrix

  it "high-level: solve linear equations" do
    a = Matrix(Float32).new(
      [[2, 4],
       [2, 8]]
    )
    b = Matrix(Float32).new([[2], [4]])
    solve(a, b).should eq (a.inv * b)
  end

  it "high-level: calculate determinant" do
    a = Matrix(Float64).new(
      [[1, 2, 3],
       [4, 5, 7],
       [-1, 1, -1]]
    )
    a.det.should eq 9
  end
end
