require "benchmark"
require "../src/matrix/*"
require "../src/linalg/*"

include Linalg

# Conclusion: incremental is predictably faster

def pascal_expm(n)
  upper = GMat.new(n, n) { |r, c| r == c - 1 ? r + 1 : 0 }
  upper.assume! MatrixFlags::UpperTriangular
  lower = GMat.new(n, n) { |r, c| r == c + 1 ? r : 0 }
  lower.assume! MatrixFlags::LowerTriangular
  lower.expm * upper.expm
end

def n_choose_k(n, k)
  (1..n).product(1.0) / ((1..k).product(1.0) * (1..n - k).product(1.0))
end

def pascal_nchoosek(n)
  GMat.new(n, n) { |row, col| n_choose_k(row + col, row) }
end

def pascal_incremental(n)
  m = Mat.ones(n, n)
  m.each_index do |i, j|
    next if i == 0 || j == 0
    m.unsafe_set(i, j, m.unsafe_at(i - 1, j) + m.unsafe_at(i, j - 1))
  end
  m
end

pp pascal_expm(5), pascal_nchoosek(5), pascal_incremental(5)
pp pascal_expm(50)[-1, -1], pascal_nchoosek(50)[-1, -1], pascal_incremental(50)[-1, -1]

Benchmark.ips do |bench|
  bench.report("expm") { pascal_expm(50) }
  bench.report("pascal_nchoosek") { pascal_nchoosek(50) }
  bench.report("pascal_incremental") { pascal_incremental(50) }
end

# pascal_expm(5)        # => Linalg::GeneralMatrix(Float64) (5x5, None):
# [1.0, 1.0, 1.0, 1.0, 1.0]
# [1.0, 2.0, 3.0, 4.0, 5.0]
# [1.0, 3.0, 6.0, 10.0, 15.0]
# [1.0, 4.0, 10.0, 20.0, 35.0]
# [1.0, 5.0, 15.0, 35.0, 70.0]
#
#
# pascal_nchoosek(5)    # => Linalg::GeneralMatrix(Float64) (5x5, None):
# [1.0, 1.0, 1.0, 1.0, 1.0]
# [1.0, 2.0, 3.0, 4.0, 5.0]
# [1.0, 3.0, 6.0, 10.0, 15.0]
# [1.0, 4.0, 10.0, 20.0, 35.0]
# [1.0, 5.0, 15.0, 35.0, 70.0]
#
#
# pascal_incremental(5) # => Linalg::GeneralMatrix(Float64) (5x5, None):
# [1.0, 1.0, 1.0, 1.0, 1.0]
# [1.0, 2.0, 3.0, 4.0, 5.0]
# [1.0, 3.0, 6.0, 10.0, 15.0]
# [1.0, 4.0, 10.0, 20.0, 35.0]
# [1.0, 5.0, 15.0, 35.0, 70.0]
#
#
# (pascal_expm(50))[-1, -1]        # => 0.99999999999999645
# (pascal_nchoosek(50))[-1, -1]    # => 2.5477612258980845e+28
# (pascal_incremental(50))[-1, -1] # => 2.5477612258980867e+28
#               expm  20.26  ( 49.36ms) (±12.95%) 1468.52× slower
#    pascal_nchoosek   2.36k (423.77µs) (± 7.05%)   12.61× slower
# pascal_incremental  29.75k ( 33.62µs) (±13.14%)         fastest
