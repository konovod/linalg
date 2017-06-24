require "benchmark"
require "../src/linalg"

METHODS = {
  Linalg::LSMethod::QR,
  Linalg::LSMethod::Orthogonal,
  Linalg::LSMethod::SVD,
}

def test(n)
  basea = Linalg::Mat.rand(n/2, n/2)*10
  until basea.det.abs > 0.01
    basea = Linalg::Mat.rand(n/2, n/2)*10
  end
  a = Linalg::GMat.new(n/2, n) { |row, column| column < n/2 ? basea[row, column] : rand }
  b = Linalg::Mat.rand(n/2, 1)
  puts "*********N = #{n}*************"
  Benchmark.ips do |bench|
    METHODS.each do |method|
      bench.report(method.to_s) { Linalg.lstsq(a, b, method) }
    end
  end
end

test 10
test 100
test 500
test 1000
test 2000

# *********N = 10*************
#         QR 141.47k (  7.07µs) (±17.05%)       fastest
# Orthogonal  72.27k ( 13.84µs) (±17.57%)  1.96× slower
#        SVD  38.63k ( 25.88µs) (± 8.40%)  3.66× slower
# *********N = 100*************
#         QR   1.88k (531.03µs) (±13.97%)       fastest
# Orthogonal   1.07k (933.43µs) (± 9.57%)  1.76× slower
#        SVD 482.19  (  2.07ms) (± 6.12%)  3.91× slower
# *********N = 500*************
#         QR  20.99  ( 47.65ms) (± 2.68%)       fastest
# Orthogonal  13.39  ( 74.68ms) (± 1.28%)  1.57× slower
#        SVD   9.63  (103.79ms) (± 2.39%)  2.18× slower
# *********N = 1000*************
#         QR   2.72  (367.56ms) (± 2.49%)       fastest
# Orthogonal   1.68  (595.04ms) (± 0.32%)  1.62× slower
#        SVD   1.39  (719.77ms) (± 0.89%)  1.96× slower
# *********N = 2000*************
#         QR   0.34  (  2.98s ) (± 0.27%)       fastest
# Orthogonal    0.2  (  5.05s ) (± 0.00%)  1.70× slower
#        SVD   0.17  (  5.73s ) (± 0.00%)  1.93× slower
