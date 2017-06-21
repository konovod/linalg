require "benchmark"
require "../src/lapack"

METHODS = {
  LAPACK::LSMethod::QR,
  LAPACK::LSMethod::Orthogonal,
  LAPACK::LSMethod::SVD,
}

def test(n)
  basea = LAPACK::Matrix(Float32).rand(n/2, n/2)*10
  until basea.det.abs > 0.01
    basea = LAPACK::Matrix(Float32).rand(n/2, n/2)*10
  end
  a = LAPACK::Matrix(Float32).new(n/2, n) { |row, column| column < n/2 ? basea[row, column] : rand }
  b = LAPACK::Matrix(Float32).rand(n/2, 1)
  puts "*********N = #{n}*************"
  Benchmark.ips do |bench|
    METHODS.each do |method|
      bench.report(method.to_s) { LAPACK.lstsq(a, b, method) }
    end
  end
end

test 10
test 100
test 500
test 1000
test 2000

# *********N = 10*************8
#         QR 179.05k (  5.59µs) (±16.34%)       fastest
# Orthogonal  87.19k ( 11.47µs) (±13.85%)  2.05× slower
#        SVD  52.21k ( 19.15µs) (±11.53%)  3.43× slower
# *********N = 100*************8
#         QR   2.26k (442.79µs) (±11.14%)       fastest
# Orthogonal   1.37k (732.02µs) (± 6.69%)  1.65× slower
#        SVD 696.62  (  1.44ms) (± 6.72%)  3.24× slower
# *********N = 500*************8
#         QR  22.86  ( 43.74ms) (± 5.45%)       fastest
# Orthogonal  14.39  ( 69.51ms) (± 1.58%)  1.59× slower
#        SVD  10.82  ( 92.46ms) (± 2.59%)  2.11× slower
# *********N = 1000*************8
#         QR   3.11  (321.32ms) (± 1.44%)       fastest
# Orthogonal   1.93  (517.45ms) (± 0.70%)  1.61× slower
#        SVD   1.57  (637.89ms) (± 0.20%)  1.99× slower
# *********N = 2000*************8
#         QR   0.36  (  2.78s ) (± 1.75%)       fastest
# Orthogonal   0.22  (  4.51s ) (± 0.23%)  1.62× slower
#        SVD    0.2  (  5.04s ) (± 0.00%)  1.81× slower
