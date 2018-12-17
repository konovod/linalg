require "./matrix"
require "./submatrix"

module LA
  # generic matrix, heap-allocated
  class GeneralMatrix(T) < Matrix(T)
    getter raw : Slice(T)
    getter nrows : Int32
    getter ncolumns : Int32
    property flags = MatrixFlags::None

    def initialize(@nrows, @ncolumns, @flags = MatrixFlags::None)
      check_type
      @raw = Slice(T).new(nrows*ncolumns, T.new(0))
    end

    def initialize(@nrows, @ncolumns, @flags = MatrixFlags::None, &block)
      check_type
      @raw = Slice(T).new(@nrows*@ncolumns) do |index|
        i = index % @nrows
        j = index / @nrows
        T.new(yield(i, j))
      end
    end

    def initialize(values : Indexable)
      check_type
      @nrows = values.size
      @ncolumns = values[0].size
      @raw = Slice(T).new(nrows*ncolumns) do |index|
        i = index % @nrows
        j = index / @nrows
        raise IndexError.new("All rows must have same size") if j == 0 && values[i].size != @ncolumns
        T.new(values[i][j])
      end
    end

    def self.new(matrix : GeneralMatrix(T))
      new(matrix.nrows, matrix.ncolumns, matrix.raw, true, matrix.flags)
    end

    def self.new(matrix : Matrix)
      new(matrix.nrows, matrix.ncolumns, matrix.flags) do |i, j|
        matrix.unsafe_fetch(i, j)
      end
    end

    def initialize(@nrows, @ncolumns, values : Indexable, col_major = false, @flags = MatrixFlags::None)
      check_type
      if col_major
        @raw = Slice(T).new(nrows*ncolumns) { |i| T.new(values[i]) }
      else
        @raw = Slice(T).new(nrows*ncolumns) do |i|
          row = i % nrows
          col = i / nrows
          T.new(values[col + row*@ncolumns])
        end
      end
    end

    def unsafe_fetch(i, j)
      @raw.unsafe_fetch(i + j*nrows)
    end

    def unsafe_set(i, j, value)
      clear_flags # TODO - not always?
      @raw[i + j*nrows] = T.new(value)
    end

    def dup
      GeneralMatrix(T).new(@nrows, @ncolumns, @raw, true, @flags)
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
            a = unsafe_fetch(i, j)
            unsafe_set(i, j, unsafe_fetch(j, i))
            unsafe_set(j, i, a)
          end
        end
      elsif nrows == 1 || ncolumns == 1
        @nrows, @ncolumns = @ncolumns, @nrows
      else
        # TODO https://en.wikipedia.org/wiki/In-place_matrix_transposition
        newraw = Slice(T).new(nrows*ncolumns, T.new(0))
        each_with_index { |v, r, c| newraw[c + r*ncolumns] = v }
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

    # changes nrows and ncolumns of matrix (total number of elements must not change)
    def reshape!(anrows, ancolumns, col_major = false)
      return self if anrows == nrows && ancolumns == ncolumns
      raise ArgumentError.new("number of elements must not change") if anrows*ancolumns != @raw.size
      unless col_major
        anew = reshape(anrows, ancolumns, col_major)
        @raw = anew.raw
      end
      @nrows = anrows
      @ncolumns = ancolumns
      self
    end

    # changes nrows and ncolumns of matrix (total number of elements must not change)
    def reshape(anrows, ancolumns, col_major = false)
      if col_major
        clone.reshape!(anrows, ancolumns, col_major)
      else
        raise ArgumentError.new("number of elements must not change") if anrows*ancolumns != @raw.size
        GeneralMatrix(T).new(anrows, ancolumns) do |i, j|
          row_index = i*ancolumns + j
          arow = row_index / @ncolumns
          acol = row_index % @ncolumns
          unsafe_fetch(arow, acol)
        end
      end
    end

    def resize!(anrows, ancolumns)
      # TODO - better implementation?
      return self if anrows == @nrows && ancolumns == @ncolumns
      anew = GeneralMatrix(T).new(anrows, ancolumns) do |i, j|
        if j >= 0 && j < ncolumns && i >= 0 && i < nrows
          unsafe_fetch(i, j)
        else
          T.new(0)
        end
      end
      @nrows = anrows
      @ncolumns = ancolumns
      @raw = anew.raw
      clear_flags
      self
    end

    def to_a(col_major = false)
      if col_major
        @raw.to_a
      else
        Array(T).new(@ncolumns*@nrows) do |i|
          row = i / @ncolumns
          col = i % @ncolumns
          unsafe_fetch(row, col)
        end
      end
    end

    def to_aa(col_major = false)
      if col_major
        Array(Array(T)).new(@ncolumns) do |i|
          Array(T).build(@nrows) do |pointer|
            pointer.to_slice(@nrows).copy_from(@raw[i*@nrows, @nrows])
            @nrows
          end
        end
      else
        Array(Array(T)).new(@nrows) do |i|
          Array(T).new(@ncolumns) do |j|
            unsafe_fetch(i, j)
          end
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

    def self.diag(nrows, ncolumns, values)
      new(nrows, ncolumns, MatrixFlags.for_diag(nrows == ncolumns)) do |i, j|
        i == j ? values[i] : 0
      end
    end

    def self.diag(nrows, ncolumns, &block)
      new(nrows, ncolumns, MatrixFlags.for_diag(nrows == ncolumns)) do |i, j|
        i == j ? yield(i) : 0
      end
    end
  end

  alias GMat = GeneralMatrix(Float64)
  alias GMat32 = GeneralMatrix(Float32)
  alias GMatComplex = GeneralMatrix(Complex)
end
