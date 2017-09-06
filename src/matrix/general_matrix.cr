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
        i = index / @ncolumns # TODO fixing highlight in sublime /
        j = index % @ncolumns
        T.new(yield(i, j))
      end
    end

    def initialize(values : Indexable)
      check_type
      @nrows = values.size
      @ncolumns = values[0].size
      @raw = Slice(T).new(nrows*ncolumns) do |index|
        i = index / @ncolumns # TODO fixing highlight in sublime /
        j = index % @ncolumns
        raise IndexError.new("All rows must have same size") if j == 0 && values[i].size != @ncolumns
        T.new(values[i][j])
      end
    end

    def self.from(matrix)
      new(matrix.nrows, matrix.ncolumns, matrix.raw)
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
      oldflags = flags
      map_with_index! { |v, i, j| i < j - k ? 0 : v }
      self.flags = oldflags.tril(k <= 0, square?)
      self
    end

    # like a triu in scipy - remove all elements below k-diagonal
    def triu!(k = 0)
      oldflags = flags
      map_with_index! { |v, i, j| i > j - k ? 0 : v }
      self.flags = oldflags.triu(k >= 0, square?)
      self
    end
  end

  alias GMat = GeneralMatrix(Float64)
  alias GMat32 = GeneralMatrix(Float32)
  alias GMatComplex = GeneralMatrix(Complex)
end
