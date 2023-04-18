require "benchmark"
require "simplepool"
require "../src/matrix/*"
include LA

# Conclusion: using Set only hurts performance

abstract class LA::Sparse::Matrix(T)
  private TEST_POOL = SimplePool(Set(LA::Utils::RowColumn)).new

  def compare_all(other : Sparse::Matrix(T))
    return false unless nrows == other.nrows && ncolumns == other.ncolumns
    a, b = self, other
    a, b = b, a if a.nonzeros < b.nonzeros
    a.each_with_index(all: false) do |v, i, j|
      return false if b.unsafe_fetch(i, j) != v
    end
    b.each_with_index(all: false) do |v, i, j|
      return false if a.unsafe_fetch(i, j) != v
    end
    true
  end

  def compare_with_set(other : Sparse::Matrix(T))
    return false unless nrows == other.nrows && ncolumns == other.ncolumns
    a, b = self, other
    a, b = b, a if a.nonzeros < b.nonzeros
    TEST_POOL.use do |set|
      set.clear
      a.each_with_index(all: false) do |v, i, j|
        return false if b.unsafe_fetch(i, j) != v
        set.add({i, j})
      end
      b.each_index(all: false) do |i, j|
        next if set.includes?({i, j})
        return false if a.unsafe_fetch(i, j) != self.unsafe_fetch(i, j)
      end
    end
    true
  end
end

def test_matrix(size, count, rng)
  m = LA::Sparse::COOMatrix(Float64).new(size, size)
  count.times do
    m[rng.rand(size), rng.rand(size)] = rng.rand
  end
  m
end

def test(size)
  rng = Random::PCG32.new(123)
  m1 = test_matrix(size, size//10, rng)
  m1[size - 1, size - 1] = 2.0
  m2 = m1.clone
  m_almost1 = m1.clone
  m_almost1[size - 1, size - 1] = 3.0
  rng = Random::PCG32.new(123)
  m_similar = test_matrix(size, size//5, rng)
  m_other = test_matrix(size, size//5, rng)

  puts "*********N = #{size}*************"
  puts "Same:"
  Benchmark.ips do |bench|
    bench.report("no set") { raise "" unless m1.compare_all(m2) }
    bench.report("set") { raise "" unless m1.compare_with_set(m2) }
  end
  puts "Almost same:"
  Benchmark.ips do |bench|
    bench.report("no set") { raise "" if m1.compare_all(m_almost1) }
    bench.report("set") { raise "" if m1.compare_with_set(m_almost1) }
  end
  puts "Similar"
  Benchmark.ips do |bench|
    bench.report("no set") { raise "" if m1.compare_all(m_similar) }
    bench.report("set") { raise "" if m1.compare_with_set(m_similar) }
  end
  puts "Different"
  Benchmark.ips do |bench|
    bench.report("no set") { raise "" if m1.compare_all(m_other) }
    bench.report("set") { raise "" if m1.compare_with_set(m_other) }
  end
end

test(500)

# *********N = 50*************
# Same:
# no set  11.87M ( 84.23ns) (± 1.02%)  0.0B/op        fastest
#    set   6.80M (147.03ns) (± 1.11%)  0.0B/op   1.75× slower
# Almost same:
# no set  22.71M ( 44.03ns) (± 0.74%)  0.0B/op        fastest
#    set   9.19M (108.84ns) (± 1.77%)  0.0B/op   2.47× slower
# Similar
# no set  23.01M ( 43.47ns) (± 0.84%)  0.0B/op        fastest
#    set   9.37M (106.74ns) (± 1.04%)  0.0B/op   2.46× slower
# Different
# no set  80.41M ( 12.44ns) (± 2.24%)  0.0B/op        fastest
#    set  23.40M ( 42.73ns) (± 1.05%)  0.0B/op   3.44× slower
