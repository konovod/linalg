require "./matrix"

module Linalg
  # generic matrix, heap-allocated
  # TODO - constructing from matlab-like string "[1,2,3;3,4,6;1,1,3]" (check regexps?)
  class GeneralMatrix(T)
    include Matrix(T)
    getter rows : Int32
    getter columns : Int32
    getter raw : Slice(T)
    property flags : MatrixFlags = MatrixFlags.new(0)

    private def check_type
      {% unless T == Float32 || T == Float64 || T == Complex %}
      {% raise "Wrong matrix members type: #{T}. Types supported by Linalg are: #{SUPPORTED_TYPES}" %}
    {% end %}
    end

    def initialize(@rows, @columns)
      check_type
      @raw = Slice(T).new(rows*columns, T.new(0))
    end

    def initialize(@rows, @columns, values, @flags = MatrixFlags.new(0))
      check_type
      @raw = Slice(T).new(rows*columns) { |i| T.new(values[i]) }
    end

    def initialize(values)
      check_type
      @rows = values.size
      @columns = values[0].size
      @raw = Slice(T).new(rows*columns) do |index|
        i = index / @columns
        j = index % @columns
        raise IndexError.new("All rows should have same size") if j == 0 && values[i].size != @columns
        T.new(values[i][j])
      end
    end

    def self.from(matrix)
      new(matrix.rows, matrix.columns, matrix.raw)
    end

    def initialize(@rows, @columns, &block)
      check_type
      @raw = Slice(T).new(@rows*@columns) do |index|
        i = index / @columns
        j = index % @columns
        T.new(yield(i, j))
      end
    end

    def unsafe_at(i, j)
      @raw.unsafe_at(i*columns + j)
    end

    def unsafe_set(i, j, value)
      @raw[i*columns + j] = T.new(value)
    end

    # TODO - benchmark is it faster?
    def dup
      GeneralMatrix(T).new(@rows, @columns, @raw)
    end

    def clone
      dup
    end

    def to_unsafe
      {% if T == Complex %}
      @raw.to_unsafe.as(LibLAPACKE::DoubleComplex*)
    {% else %}
      @raw.to_unsafe
    {% end %}
    end

    def ==(other : self)
      @rows == other.rows && @columns == other.columns && @raw == other.raw
    end

    # transposes matrix inplace
    def transpose!
      if square?
        (0..@rows - 2).each do |i|
          (i + 1..@columns - 1).each do |j|
            a = self[i, j]
            self[i, j] = self[j, i]
            self[j, i] = a
          end
        end
      else
        # TODO https://en.wikipedia.org/wiki/In-place_matrix_transposition
        raise "not implemented yet"
      end
    end

    def reshape!(arows, acolumns)
      raise ArgumentError.new("number of elements should not change") if arows*acolumns != @raw.size
      @rows = arows
      @columns = acolumns
      self
    end

    def to_a
      @raw.to_a
    end

    def to_aa
      Array(Array(T)).new(@rows) do |i|
        Array(T).build(@columns) do |pointer|
          pointer.to_slice(@columns).copy_from(@raw[i*@columns, @columns])
          @columns
        end
      end
    end
  end

  alias GMat = GeneralMatrix(Float64)
  alias GMat32 = GeneralMatrix(Float32)
  alias GMatComplex = GeneralMatrix(Complex)
end
