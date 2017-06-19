module LAPACK
  # general matrix, heap-allocated version
  class Matrix(T)
    getter n : Int32
    getter m : Int32
    getter raw : Array(T)

    def initialize(@n, @m)
      @raw = Array(T).new(n*m, T.new(0))
    end

    def initialize(@n, @m, values)
      @raw = Array(T).new(n*m) { |i| T.new(values[i]) }
    end

    def initialize(values)
      @m = values.size
      @n = values[0].size
      @raw = Array(T).new(n*m) do |index|
        i = index / @n
        j = index % @n
        T.new(values[i][j])
      end
    end

    def initialize(@n, @m, &block)
      @raw = Array(T).new(n*m) do |index|
        i = index / @n
        j = index % @n
        T.new(yield(i, j))
      end
    end

    def [](i, j)
      @raw[i*n + j]
    end

    def []=(i, j, value)
      @raw[i*n + j] = value
    end

    def to_unsafe
      @raw.to_unsafe
    end
  end
end
