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

  # info = LibLAPACKE.sgeqrf(LibLAPACKE::ROW_MAJOR,
  #   3,
  #   3,
  #   matrix,
  #   3,
  #   out tau)
  # pp info, tau, matrix

  # // solve a system of linear equations
  # var a = [
  # 	[2, 4],
  # 	[2, 8]];
  #
  # var b = [[2],
  # 	[4]];
  #
  # var result = lapack.sgesv(a, b);
  # console.log(result.X);
  # console.log(result.P);

  a = [2.0, 4.0, 2.0, 8.0]
  b = [2.0, 4.0]

  piv = [0, 0]
  info = LibLAPACKE.dgesv(
    LibLAPACKE::ROW_MAJOR,
    2,
    1,
    a,
    2,
    piv,
    b,
    2)
  pp info, b
end
