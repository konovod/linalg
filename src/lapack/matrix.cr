require "complex"

module LAPACK
  SUPPORTED_TYPES = {Float32, Float64, Complex}

  # general matrix, heap-allocated
  class Matrix(T)
    getter rows : Int32
    getter columns : Int32
    getter raw : Array(T)

    def check_type
      {% unless T == Float32 || T == Float64 || T == Complex %}
        {% raise "Wrong matrix members type: #{T}. Types supported by LAPACK are: #{SUPPORTED_TYPES}" %}
      {% end %}
    end

    def initialize(@rows, @columns)
      check_type
      @raw = Array(T).new(rows*columns, T.new(0))
    end

    def initialize(@rows, @columns, values)
      check_type
      @raw = Array(T).new(rows*columns) { |i| T.new(values[i]) }
    end

    def initialize(values)
      check_type
      @rows = values.size
      @columns = values[0].size
      @raw = Array(T).new(rows*columns) do |index|
        i = index / @columns
        j = index % @columns
        raise IndexError.new("All rows should have same size") if j == 0 && values[i].size != @columns
        T.new(values[i][j])
      end
    end

    def initialize(@rows, @columns, &block)
      check_type
      @raw = Array(T).new(@rows*@columns) do |index|
        i = index / @columns
        j = index % @columns
        T.new(yield(i, j))
      end
    end

    def [](i, j)
      @raw[i*columns + j]
    end

    def []=(i, j, value)
      @raw[i*columns + j] = value
    end

    def to_unsafe
      @raw.to_unsafe
    end
  end
end
