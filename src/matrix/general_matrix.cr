require "./matrix"
require "./submatrix"

module LA
  # generic matrix, heap-allocated
  # TODO - constructing from matlab-like string "[1,2,3;3,4,6;1,1,3]" (check regexps?)
  class GeneralMatrix(T)
    include Matrix(T)
    getter nrows : Int32
    getter ncolumns : Int32
    getter raw : Slice(T)
    property flags : MatrixFlags = MatrixFlags.new(0)

    def initialize(@nrows, @ncolumns, *, @flags = MatrixFlags.new(0))
      check_type
      @raw = Slice(T).new(nrows*ncolumns, T.new(0))
    end

    def initialize(@nrows, @ncolumns, values, *, @flags = MatrixFlags.new(0))
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

    def initialize(@nrows, @ncolumns, *, @flags = MatrixFlags.new(0), &block)
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
      GeneralMatrix(T).new(@nrows, @ncolumns, @raw, flags: @flags)
    end

    def clone
      dup
    end

    def to_unsafe
      {% if T == Complex %}
      @raw.to_unsafe.as(LibCBLAS::ComplexDouble*)
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
      elsif nrows == 1 || ncolumns == 1
        @nrows, @ncolumns = @ncolumns, @nrows
      else
        # TODO https://en.wikipedia.org/wiki/In-place_matrix_transposition
        newraw = Slice(T).new(nrows*ncolumns, T.new(0))
        each_with_index { |v, r, c| newraw[c*nrows + r] = v }
        @raw = newraw
        @nrows, @ncolumns = @ncolumns, @nrows
      end
      self.flags = flags.transpose
      self
    end

    # transposes matrix inplace
    def conjtranspose!
      {% if T != Complex %}
        return transpose!
      {% end %}
      return self if flags.hermitian?
      map! &.conj
      transpose!
    end

    def reshape!(anrows, ancolumns)
      raise ArgumentError.new("number of elements must not change") if anrows*ancolumns != @raw.size
      clear_flags unless anrows == nrows && ancolumns == ncolumns
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

    def cat!(other : Matrix(T), dimension)
      # TODO - better implementation?
      anew = self.cat(other, dimension)
      @nrows = anew.nrows
      @ncolumns = anew.ncolumns
      @raw = anew.raw
      clear_flags
      self
    end

    def vcat!(other)
      cat! other, 0
    end

    def hcat!(other)
      cat! other, 1
    end

    def self.rows(*args)
      new(args)
    end

    def self.columns(*args)
      new(args).transpose!
    end

    def self.[](*args)
      new(args)
    end

    def map!(&block)
      each_with_index { |v, i, j| unsafe_set(i, j, yield(v)) }
      self
    end

    def map_with_index!(&block)
      each_with_index { |v, i, j| unsafe_set(i, j, yield(v, i, j)) }
      self
    end

    # like a tril in scipy - remove all elements above k-diagonal
    def tril!(k = 0)
      map_with_index! { |v, i, j| i < j - k ? 0 : v }
      self.flags = MatrixFlags::LowerTriangular if k <= 0
      self
    end

    # like a triu in scipy - remove all elements below k-diagonal
    def triu!(k = 0)
      map_with_index! { |v, i, j| i > j - k ? 0 : v }
      self.flags = MatrixFlags::UpperTriangular if k >= 0
      self
    end
  end

  alias GMat = GeneralMatrix(Float64)
  alias GMat32 = GeneralMatrix(Float32)
  alias GMatComplex = GeneralMatrix(Complex)
end
