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
        i = index/@ncolumns
        j = index % @ncolumns
        T.new(yield(i, j))
      end
    end

    def initialize(values : Indexable)
      check_type
      @nrows = values.size
      @ncolumns = values[0].size
      @raw = Slice(T).new(nrows*ncolumns) do |index|
        i = index/@ncolumns
        j = index % @ncolumns
        raise IndexError.new("All rows must have same size") if j == 0 && values[i].size != @ncolumns
        T.new(values[i][j])
      end
    end

    def self.new(matrix : GeneralMatrix(T))
      new(matrix.nrows, matrix.ncolumns, matrix.raw, matrix.flags)
    end

    def self.new(matrix : Matrix)
      new(matrix.nrows, matrix.ncolumns, matrix.flags) do |i, j|
        matrix.unsafe_at(i, j)
      end
    end

    def initialize(@nrows, @ncolumns, values : Indexable, @flags = MatrixFlags::None)
      check_type
      @raw = Slice(T).new(nrows*ncolumns) { |i| T.new(values[i]) }
    end

    def unsafe_at(i, j)
      @raw.unsafe_at(i*ncolumns + j)
    end

    def unsafe_set(i, j, value)
      clear_flags # TODO - not always?
      @raw[i*ncolumns + j] = T.new(value)
    end

    def dup
      GeneralMatrix(T).new(@nrows, @ncolumns, @raw, @flags)
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
            a = unsafe_at(i, j)
            unsafe_set(i, j, unsafe_at(j, i))
            unsafe_set(j, i, a)
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

    # changes nrows and ncolumns of matrix (total number of elements must not change)
    def reshape!(anrows, ancolumns)
      raise ArgumentError.new("number of elements must not change") if anrows*ancolumns != @raw.size
      clear_flags unless anrows == nrows && ancolumns == ncolumns
      @nrows = anrows
      @ncolumns = ancolumns
      self
    end

    # changes nrows and ncolumns of matrix (total number of elements must not change)
    def reshape(anrows, ancolumns)
      clone.reshape!(anrows, ancolumns)
    end

    def resize!(anrows, ancolumns)
      # TODO - better implementation?
      anew = GeneralMatrix(T).new(anrows, ancolumns) do |i, j|
        if j >= 0 && j < ncolumns && i >= 0 && i < nrows
          unsafe_at(i, j)
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
