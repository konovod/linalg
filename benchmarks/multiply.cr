require "benchmark"
require "../src/matrix/*"
require "../src/linalg/*"

include LA

# Conclusions:
# - openblas gemm is significiantly faster then naive
# - symm is faster only at size ~50
# - trmm is faster at size >= 50

module LA::Matrix(T)
  def naive_mult(m : Matrix(T))
    if ncolumns != m.nrows
      raise ArgumentError.new("matrix size should match ([#{nrows}x#{ncolumns}] * [#{m.nrows}x#{m.ncolumns}]")
    end
    result = GeneralMatrix(T).new(nrows, m.ncolumns, self.flags.mult(m.flags)) do |i, j|
      (0...ncolumns).sum { |k| self[i, k]*m[k, j] }
    end
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

  def tr_mult(m : Matrix(T))
    if ncolumns != m.nrows
      raise ArgumentError.new("matrix size should match ([#{nrows}x#{ncolumns}] * [#{m.nrows}x#{m.ncolumns}]")
    end
    if (square? && flags.triangular?) || (m.square? && m.flags.triangular?)
      if m.square? && m.flags.triangular?
        result = self.clone
        result.tr_mult!(m, left: false)
        result
      else
        result = m.clone
        result.tr_mult!(self, left: true)
        result
      end
    else
      raise "tr_mult inapplicable"
    end
  end
end

def test(n)
  a = Mat.rand(n, n)
  b = Mat.rand(n, n)
  asy = a + a.t
  atr = a.triu
  raise "noo" unless asy.detect MatrixFlags::Symmetric
  puts "*********N = #{n}*************"
  Benchmark.ips do |bench|
    bench.report("naive") { a.naive_mult(b) }
    bench.report("ge_mult") { a.ge_mult(b) }
    # bench.report("ge_mult2") { a.ge_mult2(b) }
    bench.report("sy_mult") { asy.ge_mult(b) }
    bench.report("tr_mult") { atr.tr_mult(b) }
  end
end

# check methods consistency
a = Mat.rand(10, 10)
a1 = a[0..5, 0..5]
asy = a1 + a1.t
atr = a1.tril
raise "* failure" unless (a1*a1.inv).almost_eq Mat.eye(a1.nrows)
raise "naive failure" unless (a1.naive_mult(a1.inv)).almost_eq Mat.eye(a1.nrows)
raise "ge_mult failure" unless (a1.ge_mult(a1.inv)).almost_eq Mat.eye(a1.nrows)
raise "sy_mult failure" unless (asy.ge_mult(asy.inv)).almost_eq Mat.eye(asy.nrows)
raise "tr_mult failure" unless (atr.tr_mult(atr.inv)).almost_eq Mat.eye(atr.nrows)

test(10)
test(50)
test(100)
test(200)
test(500)
# test(1000)

# *********N = 10*************
#   naive   60.2k ( 16.61µs) (±14.10%)  7.13× slower
# ge_mult 429.02k (  2.33µs) (±19.31%)       fastest
# sy_mult 417.76k (  2.39µs) (±20.65%)  1.03× slower
# tr_mult 229.91k (  4.35µs) (±13.02%)  1.87× slower
# *********N = 50*************
#   naive 661.76  (  1.51ms) (± 7.91%) 30.34× slower
# ge_mult  12.87k ( 77.71µs) (±14.92%)  1.56× slower
# sy_mult   17.9k ( 55.87µs) (± 5.02%)  1.12× slower
# tr_mult  20.08k ( 49.81µs) (± 6.75%)       fastest
# *********N = 100*************
#   naive  85.82  ( 11.65ms) (± 5.74%) 55.36× slower
# ge_mult   3.71k (269.74µs) (± 3.58%)  1.28× slower
# sy_mult   3.68k (271.39µs) (± 5.32%)  1.29× slower
# tr_mult   4.75k (210.49µs) (± 3.52%)       fastest
# *********N = 200*************
#   naive  10.45  ( 95.65ms) (± 6.11%) 84.78× slower
# ge_mult  593.3  (  1.69ms) (± 2.07%)  1.49× slower
# sy_mult 586.06  (  1.71ms) (± 6.23%)  1.51× slower
# tr_mult 886.38  (  1.13ms) (± 2.12%)       fastest
# *********N = 500*************
#   naive   0.59  (  1.71s ) (± 1.62%) 125.66× slower
# ge_mult  43.74  ( 22.86ms) (± 3.45%)   1.68× slower
# sy_mult  43.76  ( 22.85ms) (± 1.86%)   1.68× slower
# tr_mult   73.6  ( 13.59ms) (± 3.58%)        fastest
