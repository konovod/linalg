require "benchmark"
require "../src/matrix/*"

# Check if `##triu` method really needed or clone.triu! is enough.
# Conclusion: there is a difference only for high N( = 1000 ), but it still exists, so triu stays
{2, 5, 10, 100, 1000}.each do |n|
  a = Linalg::Mat.rand(n, n)
  puts "*********N = #{n}*************"
  Benchmark.ips do |bench|
    bench.report("clone.tril!") { a.clone.tril! }
    bench.report("tril") { a.tril }
    bench.report("triu") { a.triu }
    bench.report("clone.triu!") { a.clone.triu! }
  end
end
#
# *********N = 2*************
# clone.tril!    3.7M (270.39ns) (±27.64%)  1.14× slower
#        tril   3.86M (259.26ns) (±26.72%)  1.09× slower
#        triu   4.04M ( 247.8ns) (±28.74%)  1.04× slower
# clone.triu!    4.2M (237.84ns) (±29.58%)       fastest
# *********N = 5*************
# clone.tril!   1.06M (946.96ns) (±21.12%)  1.27× slower
#        tril   1.08M (929.73ns) (±22.17%)  1.24× slower
#        triu    1.1M (909.13ns) (±22.18%)  1.22× slower
# clone.triu!   1.34M (747.55ns) (±22.26%)       fastest
# *********N = 10*************
# clone.tril! 351.25k (  2.85µs) (±19.38%)       fastest
#        tril 327.48k (  3.05µs) (±19.20%)  1.07× slower
#        triu 322.97k (   3.1µs) (±19.90%)  1.09× slower
# clone.triu!  325.8k (  3.07µs) (±18.53%)  1.08× slower
# *********N = 100*************
# clone.tril!   5.43k ( 184.2µs) (±18.67%)       fastest
#        tril   5.29k (189.09µs) (±12.30%)  1.03× slower
#        triu   4.95k (202.02µs) (±14.86%)  1.10× slower
# clone.triu!   5.25k (190.41µs) (±14.78%)  1.03× slower
# *********N = 1000*************
# clone.tril!  81.02  ( 12.34ms) (± 3.17%)  1.18× slower
#        tril  95.81  ( 10.44ms) (± 2.30%)       fastest
#        triu  93.74  ( 10.67ms) (± 4.52%)  1.02× slower
# clone.triu!  81.04  ( 12.34ms) (± 2.59%)  1.18× slower
