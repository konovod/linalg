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
  puts "*********N = #{n}*************"
  Benchmark.ips do |bench|
    bench.report("naive") { a.naive_mult(b) }
    bench.report("ge_mult") { a.ge_mult(b) }
    bench.report("ge_mult2") { a.ge_mult2(b) }
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
#    naive  35.48k ( 28.19µs) (±13.41%)  9.32× slower
#  ge_mult 330.71k (  3.02µs) (±18.24%)       fastest
# ge_mult2 327.29k (  3.06µs) (±22.18%)  1.01× slower
# *********N = 50*************
#    naive  341.5  (  2.93ms) (± 4.31%) 23.16× slower
#  ge_mult   7.91k (126.43µs) (±14.15%)       fastest
# ge_mult2    7.5k (133.31µs) (±14.16%)  1.05× slower
# *********N = 100*************
#    naive  43.29  (  23.1ms) (± 5.40%) 50.21× slower
#  ge_mult   2.17k (460.11µs) (± 2.90%)       fastest
# ge_mult2   2.15k ( 465.6µs) (± 2.62%)  1.01× slower
# *********N = 200*************
#    naive   5.18  (193.06ms) (± 4.35%) 59.46× slower
#  ge_mult 307.52  (  3.25ms) (± 7.27%)  1.00× slower
# ge_mult2 307.96  (  3.25ms) (± 2.17%)       fastest
# *********N = 500*************
#    naive   0.22  (  4.58s ) (± 0.59%) 103.27× slower
#  ge_mult  22.53  ( 44.39ms) (± 3.33%)        fastest
# ge_mult2  22.46  ( 44.52ms) (± 3.58%)   1.00× slower
