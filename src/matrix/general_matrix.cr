require "./matrix"
require "./submatrix"

module Linalg
  # generic matrix, heap-allocated
  # TODO - constructing from matlab-like string "[1,2,3;3,4,6;1,1,3]" (check regexps?)
  class GeneralMatrix(T)
    include Matrix(T)
    getter nrows : Int32
    getter ncolumns : Int32
    getter raw : Slice(T)
    property flags : MatrixFlags = MatrixFlags.new(0)

    def initialize(@nrows, @ncolumns)
      check_type
      @raw = Slice(T).new(nrows*ncolumns, T.new(0))
    end

    def initialize(@nrows, @ncolumns, values, @flags = MatrixFlags.new(0))
      check_type
      @raw = Slice(T).new(nrows*ncolumns) { |i| T.new(values[i]) }
    end

    def initialize(values)
      check_type
      @nrows = values.size
      @ncolumns = values[0].size
      @raw = Slice(T).new(nrows*ncolumns) do |index|
        i = index / @ncolumns
        j = index % @ncolumns
        raise IndexError.new("All rows must have same size") if j == 0 && values[i].size != @ncolumns
        T.new(values[i][j])
      end
    end

    def self.from(matrix)
      new(matrix.nrows, matrix.ncolumns, matrix.raw)
    end

    def initialize(@nrows, @ncolumns, &block)
      check_type
      @raw = Slice(T).new(@nrows*@ncolumns) do |index|
        i = index / @ncolumns
        j = index % @ncolumns
        T.new(yield(i, j))
      end
    end

    def unsafe_at(i, j)
      @raw.unsafe_at(i*ncolumns + j)
    end

    def unsafe_set(i, j, value)
      clear_flags # TODO - not always?
      @raw[i*ncolumns + j] = T.new(value)
    end

    def dup
      GeneralMatrix(T).new(@nrows, @ncolumns, @raw).tap { |it| it.flags = flags }
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
      @nrows == other.nrows && @ncolumns == other.ncolumns && @raw == other.raw
    end

    # transposes matrix inplace
    def transpose!
      return self if flags.symmetric?
      if square?
        (0..@nrows - 2).each do |i|
          (i + 1..@ncolumns - 1).each do |j|
            a = self[i, j]
            self[i, j] = self[j, i]
            self[j, i] = a
          end
        end
        @flags ^= MatrixFlags::Lower if flags.triangular?
        self
      else
        # TODO https://en.wikipedia.org/wiki/In-place_matrix_transposition
        raise "not implemented yet"
      end
    end

    # transposes matrix inplace
    def conjtranspose!
      {% if T != Complex %}
        return transpose!
      {% end %}
      return self if flags.hermitian?
      each_with_index { |v, i, j| unsafe_set(i, j, v.conj) }
      transpose!
    end

    def reshape!(anrows, ancolumns)
      raise ArgumentError.new("number of elements must not change") if anrows*ancolumns != @raw.size
      @nrows = anrows
      @ncolumns = ancolumns
      self
    end

    def to_a
      @raw.to_a
    end

    def to_aa
      Array(Array(T)).new(@nrows) do |i|
        Array(T).build(@ncolumns) do |pointer|
          pointer.to_slice(@ncolumns).copy_from(@raw[i*@ncolumns, @ncolumns])
          @ncolumns
        end
      end
    end
  end

  alias GMat = GeneralMatrix(Float64)
  alias GMat32 = GeneralMatrix(Float32)
  alias GMatComplex = GeneralMatrix(Complex)
end
