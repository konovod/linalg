require "benchmark"
require "../src/matrix/*"
require "../src/linalg/*"

include Linalg

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
end

def test(n)
  a = Mat.rand(n, n*2)
  b = Mat.rand(n*2, n)
  puts "*********N = #{n}*************"
  Benchmark.ips do |bench|
    bench.report("naive") { a.naive_mult(b) }
    bench.report("ge_mult") { a.ge_mult(b) }
  end
end

# check methods consistency
a = Mat.rand(10, 10)
a1 = a[0..5, 0..5]
raise "failure" unless (a1*a1.inv).almost_eq Mat.eye(a1.nrows)

test(10)
test(50)
test(100)
test(200)
test(500)
test(1000)

# *********N = 10*************
#   naive  35.43k ( 28.22µs) (±11.74%)  9.76× slower
# ge_mult 345.96k (  2.89µs) (±22.20%)       fastest
# *********N = 50*************
#   naive  340.6  (  2.94ms) (± 4.72%) 23.33× slower
# ge_mult   7.95k (125.85µs) (±15.33%)       fastest
# *********N = 100*************
#   naive  43.27  ( 23.11ms) (± 5.38%) 50.06× slower
# ge_mult   2.17k (461.63µs) (± 2.64%)       fastest
# *********N = 200*************
#   naive    5.2  (192.24ms) (± 5.82%) 60.85× slower
# ge_mult 316.53  (  3.16ms) (± 2.31%)       fastest
# *********N = 500*************
#   naive   0.22  (  4.54s ) (± 0.60%) 100.15× slower
# ge_mult  22.06  ( 45.34ms) (± 2.66%)        fastest
# *********N = 1000*************
# GC Warning: Repeated allocation of very large block (appr. size 8003584):
# 	May lead to memory leak and poor performance
#   naive   0.03  ( 33.12s ) (± 0.00%) 96.33× slower
# ge_mult   2.91  (343.82ms) (± 1.61%)       fastest
