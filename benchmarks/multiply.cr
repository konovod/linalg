require "benchmark"
require "../src/matrix/*"
require "../src/linalg/*"

include Linalg

# Conclusion: blas gemm is significiantly faster, symm doesn't bring any profit

module Linalg::Matrix(T)
  def naive_mult(m : Matrix(T))
    if ncolumns != m.nrows
      raise ArgumentError.new("matrix size should match ([#{nrows}x#{ncolumns}] * [#{m.nrows}x#{m.ncolumns}]")
    end
    result = GeneralMatrix(T).new(nrows, m.ncolumns) do |i, j|
      (0...ncolumns).sum { |k| self[i, k]*m[k, j] }
    end
    result.tap { |r| r.flags = self.flags.mult(m.flags) }
  end

  def ge_mult(m : Matrix(T))
    if ncolumns != m.nrows
      raise ArgumentError.new("matrix size should match ([#{nrows}x#{ncolumns}] * [#{m.nrows}x#{m.ncolumns}]")
    end
    result = Matrix(T).zeros(nrows, m.ncolumns)
    result.inc_mult(self, m)
    result.tap { |r| r.flags = self.flags.mult(m.flags) }
  end

  def ge_mult2(m : Matrix(T))
    if ncolumns != m.nrows
      raise ArgumentError.new("matrix size should match ([#{nrows}x#{ncolumns}] * [#{m.nrows}x#{m.ncolumns}]")
    end
    result = Matrix(T).zeros(nrows, m.ncolumns)
    result.inc_mult(self, m, beta: 0.0)
    result.tap { |r| r.flags = self.flags.mult(m.flags) }
  end
end

def test(n)
  a = Mat.rand(n, n)
  b = Mat.rand(n, n)
  asy = a + a.t
  raise "noo" unless asy.detect MatrixFlags::Symmetric
  puts "*********N = #{n}*************"
  Benchmark.ips do |bench|
    bench.report("naive") { a.naive_mult(b) }
    bench.report("ge_mult") { a.ge_mult(b) }
    # bench.report("ge_mult2") { a.ge_mult2(b) }
    bench.report("sy_mult") { asy.ge_mult(b) }
  end
end

# check methods consistency
a = Mat.rand(10, 10)
a1 = a[0..5, 0..5]
asy = a1 + a1.t
raise "* failure" unless (a1*a1.inv).almost_eq Mat.eye(a1.nrows)
raise "naive failure" unless (a1.naive_mult(a1.inv)).almost_eq Mat.eye(a1.nrows)
raise "ge_mult failure" unless (a1.ge_mult(a1.inv)).almost_eq Mat.eye(a1.nrows)
raise "ge_sym failure" unless (asy.ge_mult2(asy.inv)).almost_eq Mat.eye(asy.nrows)

test(10)
test(50)
test(100)
test(200)
test(500)
# test(1000)

# *********N = 10*************
#   naive  62.95k ( 15.88µs) (±13.81%)  6.64× slower
# ge_mult 417.71k (  2.39µs) (±21.83%)       fastest
# sy_mult 404.63k (  2.47µs) (±22.07%)  1.03× slower
# *********N = 50*************
#   naive 660.28  (  1.51ms) (± 8.65%) 28.44× slower
# ge_mult  12.98k ( 77.06µs) (±16.67%)  1.45× slower
# sy_mult  18.78k ( 53.26µs) (± 3.38%)       fastest
# *********N = 100*************
#   naive  85.75  ( 11.66ms) (± 5.61%) 44.96× slower
# ge_mult   3.83k (260.78µs) (± 2.89%)  1.01× slower
# sy_mult   3.86k (259.37µs) (± 2.22%)       fastest
# *********N = 200*************
#   naive  10.52  ( 95.09ms) (± 5.51%) 55.71× slower
# ge_mult 576.84  (  1.73ms) (± 7.82%)  1.02× slower
# sy_mult 585.87  (  1.71ms) (± 2.14%)       fastest
# *********N = 500*************
#   naive   0.59  (   1.7s ) (± 1.59%) 75.27× slower
# ge_mult  44.35  ( 22.55ms) (± 2.77%)       fastest
# sy_mult  43.98  ( 22.74ms) (± 4.97%)  1.01× slower
