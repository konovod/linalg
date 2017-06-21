require "complex"

struct Complex
  def self.new(value)
    case value
    when Complex
      new(value.real, value.imag)
    else
      new(value, 0.0)
    end
  end
end

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
      # i isn't checked as underlying array will check it anyway
      if j >= 0 && j < @columns
        @raw[i*columns + j]
      else
        raise IndexError.new("access to [#{i}, #{j}] in matrix with size #{@rows}x#{@columns}")
      end
    end

    def dup
      Matrix(T).new(@rows, @columns, @raw)
    end

    def clone
      dup
    end

    def []=(i, j, value)
      # i isn't checked as underlying array will check it anyway
      if j >= 0 && j < @columns
        @raw[i*columns + j] = value
      else
        raise IndexError.new("access to [#{i}, #{j}] in matrix with size #{@rows}x#{@columns}")
      end
    end

    def to_unsafe
      {% if T == Complex %}
        @raw.to_unsafe.as(LibLAPACKE::DoubleComplex*)
      {% else %}
        @raw.to_unsafe
      {% end %}
    end

    def *(m : self)
      if @columns != m.rows
        raise ArgumentError.new("matrix size should match ([#{@rows}x#{@columns}] * [#{m.rows}x#{m.columns}]")
      end
      result = Matrix(T).new(@rows, m.columns) do |i, j|
        (0...@columns).sum { |k| self.[i, k]*m[k, j] }
      end
    end

    def *(k : Int | Float)
      result = Matrix(T).new(@rows, @columns) do |i, j|
        self.[i, j]*k
      end
    end

    def ==(other : self)
      @rows == other.rows && @columns == other.columns && @raw == other.raw
    end

    def +(m : self)
      if @columns != m.columns || @rows != m.rows
        raise ArgumentError.new("matrix size should match ([#{@rows}x#{@columns}] + [#{m.rows}x#{m.columns}]")
      end
      result = Matrix(T).new(@rows, m.columns) do |i, j|
        self.[i, j] + m[i, j]
      end
    end

    def -(m : self)
      if @columns != m.columns || @rows != m.rows
        raise ArgumentError.new("matrix size should match ([#{@rows}x#{@columns}] - [#{m.rows}x#{m.columns}]")
      end
      result = Matrix(T).new(@rows, m.columns) do |i, j|
        self.[i, j] - m[i, j]
      end
    end

    def self.identity(n)
      new(n, n) { |i, j| i == j ? 1 : 0 }
    end

    def square?
      rows == columns
    end

    def self.rand(rows, columns, rng = Random::DEFAULT)
      new(rows, columns) { |i, j| rng.rand }
    end

    def self.zeros(rows, columns)
      new(rows, columns)
    end

    def self.ones(rows, columns)
      new(rows, columns) { |i, j| 1 }
    end

    def repmat(arows, acolumns)
      new(rows*arows, columns*acolumns) do |i, j|
        self[i % rows, j % columns]
      end
    end

    def self.repmat(a : self, rows, columns)
      a.repmat(rows, columns)
    end

    # TODO - proper norms
    # TODO - iteration on cols\rows
    def abs
      (0...@rows).map { |r| (0...@columns).sum { |c| self[r, c].abs } }.max
    end
  end
end
