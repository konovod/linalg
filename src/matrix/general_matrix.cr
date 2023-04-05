require "./matrix"
require "./submatrix"

module LA
  # generic matrix, heap-allocated
  #
  # Data are stored in column-major as this is a storage used by LAPACK
  #
  # See `SUPPORTED_TYPES` for supported types
  class GeneralMatrix(T) < Matrix(T)
    # Pointer to a raw data
    getter raw : Slice(T)
    # Count of rows in matrix
    getter nrows : Int32
    # Count of columns in matrix
    getter ncolumns : Int32
    # Matrix flags (see `MatrixFlags` for description)
    property flags : MatrixFlags = MatrixFlags::None

    # Creates zero-initialized matrix of given size
    #
    # Example: `LA::GMat.new(4,4)`
    def initialize(@nrows, @ncolumns, @flags = MatrixFlags::None)
      check_type
      @raw = Slice(T).new(nrows*ncolumns, T.new(0))
    end

    # Creates matrix of given size and then call block to initialize each element
    #
    # Example: `LA::GMat.new(4,4){|i,j| i+j }`
    def initialize(@nrows, @ncolumns, @flags = MatrixFlags::None, &block)
      check_type
      @raw = Slice(T).new(@nrows*@ncolumns) do |index|
        i = index % @nrows
        j = index // @nrows
        T.new(yield(i, j))
      end
    end

    # Creates matrix from any `Indexable` of `Indexable`s
    #
    # Example:
    # ```
    # m = GMat.new([
    #   [1, 2, 3],
    #   [4, 5, 6],
    #   [7, 8, 9],
    #   [10, 11, 12],
    # ])
    # ```
    def initialize(values : Indexable)
      check_type
      @nrows = values.size
      @ncolumns = values[0].size
      @raw = Slice(T).new(nrows*ncolumns) do |index|
        i = index % @nrows
        j = index // @nrows
        raise IndexError.new("All rows must have same size") if j == 0 && values[i].size != @ncolumns
        T.new(values[i][j])
      end
    end

    # :nodoc:
    def self.new(matrix : GeneralMatrix(T))
      new(matrix.nrows, matrix.ncolumns, matrix.raw, true, matrix.flags)
    end

    # Creates matrix with same content as another matrix
    def self.new(matrix : Matrix)
      new(matrix.nrows, matrix.ncolumns, matrix.flags) do |i, j|
        matrix.unsafe_fetch(i, j)
      end
    end

    # Creates matrix with given size and populate elements from `values`
    #
    # if `col_major` is true, `values` content is just copied to `#raw`,
    # otherwise a conversion from row-major form is performed
    # Example:
    # ```
    # values = [1, 2, 3, 4]
    # a = GMat.new(2, 2, values)
    # a.to_aa # => [[1,2],[3,4]]
    # b = GMat.new(2, 2, values, col_major: true)
    # b.to_aa # => [[1,3],[2,4]]
    # ```
    def initialize(@nrows, @ncolumns, values : Indexable, col_major = false, @flags = MatrixFlags::None)
      check_type
      if col_major
        @raw = Slice(T).new(nrows*ncolumns) { |i| T.new(values[i]) }
      else
        @raw = Slice(T).new(nrows*ncolumns) do |i|
          row = i % nrows
          col = i // nrows
          T.new(values[col + row*@ncolumns])
        end
      end
    end

    # returns element at row i and column j, without performing any checks
    def unsafe_fetch(i, j)
      @raw.unsafe_fetch(i + j*nrows)
    end

    # sets element at row i and column j to value, without performing any checks
    def unsafe_set(i, j, value)
      clear_flags # TODO - not always?
      @raw[i + j*nrows] = T.new(value)
    end

    # returns copy of matrix
    def dup
      GeneralMatrix(T).new(@nrows, @ncolumns, @raw, true, @flags)
    end

    # :ditto:
    def clone
      dup
    end

    # Returns pointer to raw data, suitable to e.g. use with LAPACK
    def to_unsafe
      {% if T == Complex %}
        @raw.to_unsafe.as(LibCBLAS::ComplexDouble*)
      {% else %}
        @raw.to_unsafe
      {% end %}
    end

    # see `LA::Matrix#==`
    def ==(other : self)
      @nrows == other.nrows && @ncolumns == other.ncolumns && @raw == other.raw
    end

    # transposes matrix inplace
    #
    # Currently, transpose of non-square matrix still allocates temporary buffer
    def transpose!
      return self if flags.symmetric?
      if square?
        each_upper(diagonal: false) do |v, i, j|
          unsafe_set(i, j, unsafe_fetch(j, i))
          unsafe_set(j, i, v)
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

    # Conjurgate transposes matrix inplace
    #
    # Currently, transpose of non-square matrix still allocates temporary buffer
    def conjtranspose!
      {% if T != Complex %}
        return transpose!
      {% end %}
      return self if flags.hermitian?
      map! &.conj
      transpose!
    end

    # Changes nrows and ncolumns of matrix (total number of elements must not change)
    #
    # if `col_major` is true, just nrows and ncolumns are changed, data kept the same
    # Otherwise, elements are reordered to emulate row-major storage
    # Example:
    # ```
    # a = GMat[[1, 2, 3], [4, 5, 6]]
    # a.reshape!(2, 3)
    # a.to_aa # => [[1.0, 2.0], [3.0, 4.0], [5.0, 6.0]]
    # b = GMat[[1, 2, 3], [4, 5, 6]]
    # b.reshape!(2, 3, col_major: true)
    # this is because in-memory order of values is (1, 4, 2, 5, 3, 6)
    # b.to_aa # => [[1.0, 5.0], [4.0, 3.0], [2.0, 6.0]]
    # ```
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

    # Returns a matrix with different nrows and ncolumns but same elements (total number of elements must not change)
    #
    # if `col_major` is true, just nrows and ncolumns are changed, data kept the same
    # Otherwise, elements are reordered to emulate row-major storage
    # Example:
    # ```
    # a = GMat[[1, 2, 3], [4, 5, 6]]
    # a.reshape(2, 3).to_aa # => [[1.0, 2.0], [3.0, 4.0], [5.0, 6.0]]
    # b = GMat[[1, 2, 3], [4, 5, 6]]
    # this is because in-memory order of values is (1, 4, 2, 5, 3, 6)
    # b.reshape(2, 3, col_major: true).to_aa # => [[1.0, 5.0], [4.0, 3.0], [2.0, 6.0]]
    # ```
    def reshape(anrows, ancolumns, col_major = false)
      if col_major
        clone.reshape!(anrows, ancolumns, col_major)
      else
        raise ArgumentError.new("number of elements must not change") if anrows*ancolumns != @raw.size
        GeneralMatrix(T).new(anrows, ancolumns) do |i, j|
          row_index = i*ancolumns + j
          arow = row_index // @ncolumns
          acol = row_index % @ncolumns
          unsafe_fetch(arow, acol)
        end
      end
    end

    # Change number of rows and columns in matrix.
    #
    # if new number is higher zero elements are added,
    # if new number is lower, exceeding elements are lost
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

    # Converts matrix to plain array of elements
    # if `col_major` is true, elements are returned as stored inplace,
    # otherwise row-major storage is emulated
    def to_a(col_major = false)
      if col_major
        @raw.to_a
      else
        Array(T).new(@ncolumns*@nrows) do |i|
          row = i // @ncolumns
          col = i % @ncolumns
          unsafe_fetch(row, col)
        end
      end
    end

    # Converts matrix to array of array of elements
    # if `col_major` is true, elements are returned as stored inplace,
    # otherwise row-major storage is emulated
    # Example:
    # ```
    # a = GMat[[1, 2], [3, 4]]
    # a.to_aa                  # => [[1.0, 2.0],[3.0, 4.0]]
    # a.to_aa(col_major: true) # => [[1.0, 3.0],[2.0, 4.0]]
    # ```
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

    # Concatenate matrix adding another matrix by `dimension` `Axis::Rows` (horizontal) or `Axis::Columns` (vertical)
    def cat!(other : Matrix(T), dimension)
      # TODO - better implementation?
      anew = self.cat(other, dimension)
      @nrows = anew.nrows
      @ncolumns = anew.ncolumns
      @raw = anew.raw
      clear_flags
      self
    end

    # Concatenate matrix adding another matrix vertically (so they form a column)
    def vcat!(other)
      cat! other, 0
    end

    # Concatenate matrix adding another matrix horizontally (so they form a row)
    def hcat!(other)
      cat! other, 1
    end

    # Creates matrix from a number of rows
    #
    # Example:
    # ```
    # a = GMat.rows([1, 2, 3, 4], [5, 6, 7, 8])
    # a.to_aa # => [[1,2,3,4], [5,6,7,8]]
    # ```
    def self.rows(*args)
      new(args)
    end

    # Creates matrix from a number of columns
    #
    # Example:
    # ```
    # a = GMat.columns([1, 2, 3, 4], [5, 6, 7, 8])
    # a.to_aa # => [[1, 5], [2, 6], [3, 7], [4, 8]]
    # ```
    def self.columns(*args)
      new(args).transpose!
    end

    # Alias for `#new`
    def self.[](*args)
      new(args)
    end

    # Returns diagonal matrix of given size with diagonal elements taken from array `values`
    def self.diag(nrows, ncolumns, values)
      new(nrows, ncolumns, MatrixFlags.for_diag(nrows == ncolumns)) do |i, j|
        i == j ? values[i] : 0
      end
    end

    # Returns diagonal matrix of given size with diagonal elements equal to block value
    def self.diag(nrows, ncolumns, &block)
      new(nrows, ncolumns, MatrixFlags.for_diag(nrows == ncolumns)) do |i, j|
        i == j ? yield(i) : 0
      end
    end
  end

  module Aliases
    alias GMat = GeneralMatrix(Float64)
    alias GMat32 = GeneralMatrix(Float32)
    alias GMatComplex = GeneralMatrix(Complex)
  end
end
