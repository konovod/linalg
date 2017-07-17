require "benchmark"
require "../src/linalg"

METHODS = {
  LA::LSMethod::QR,
  LA::LSMethod::Orthogonal,
  LA::LSMethod::SVD,
}

def test(n)
  basea = LA::Mat.rand(n/2, n/2)*10
  until basea.det.abs > 0.01
    basea = LA::Mat.rand(n/2, n/2)*10
  end
  a = LA::GMat.new(n/2, n) { |row, column| column < n/2 ? basea[row, column] : rand }
  b = LA::Mat.rand(n/2, 1)
  puts "*********N = #{n}*************"
  Benchmark.ips do |bench|
    METHODS.each do |method|
      bench.report(method.to_s) { LA.lstsq(a, b, method) }
    end
  end
end

test 10
test 100
test 500
test 1000
test 2000

# with OpenBLAS
# *********N = 10*************
#         QR  98.49k ( 10.15µs) (±13.47%)       fastest
# Orthogonal  51.43k ( 19.45µs) (± 9.68%)  1.92× slower
#        SVD  22.79k ( 43.87µs) (± 6.46%)  4.32× slower
# *********N = 100*************
#         QR   2.28k (439.27µs) (±12.25%)       fastest
# Orthogonal   1.28k (778.68µs) (±11.98%)  1.77× slower
#        SVD 292.27  (  3.42ms) (±10.73%)  7.79× slower
# *********N = 500*************
#         QR  40.09  ( 24.94ms) (± 3.20%)       fastest
# Orthogonal  26.61  ( 37.58ms) (± 3.59%)  1.51× slower
#        SVD  11.88  ( 84.18ms) (± 5.46%)  3.38× slower
# *********N = 1000*************
#         QR    7.9  (126.63ms) (± 5.62%)       fastest
# Orthogonal   4.13  ( 241.9ms) (± 1.29%)  1.91× slower
#        SVD   2.94  (339.71ms) (± 1.21%)  2.68× slower
# *********N = 2000*************
# GC Warning: Repeated allocation of very large block (appr. size 16003072):
# 	May lead to memory leak and poor performance
#         QR   1.57  ( 638.6ms) (± 0.66%)       fastest
# Orthogonal   0.56  (   1.8s ) (± 1.10%)  2.81× slower
#        SVD   0.53  (  1.89s ) (± 0.29%)  2.96× slower
