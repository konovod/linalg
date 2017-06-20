require "./spec_helper"

describe LAPACK do
  # TODO: Write tests

  # var result = lapack.sgeqrf([
  #     [1, 2, 3],
  #     [3, 4, 5],
  #     [5, 6, 7]
  # ]);
  matrix = Matrix(Float32).new([
    [1, 2, 3],
    [3, 4, 5],
    [5, 6, 7],
  ])
  tau = 0.0_f32
  info = LibLAPACKE.sgeqrf(LibLAPACKE::ROW_MAJOR,
    3,
    3,
    matrix,
    3,
    pointerof(tau))
  pp info, tau, matrix

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
