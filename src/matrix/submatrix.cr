require "./matrix"

module LA
  # tuple of coordinates - row and column
  alias Utils::RowColumn = {Int32, Int32}

  # it's like Slice, but for matrices.
  # Lightweight class that stores only pointer to basic matrix, size and offset
  # But allow all operations that allow `Matrix(T)`
  class SubMatrix(T) < Matrix(T)
    include DenseMatrix
    # index of a submatrix first row and column
    getter offset : Utils::RowColumn
    # number of rows in submatrix
    getter nrows : Int32
    # number of columns in submatrix
    getter ncolumns : Int32

    # :nodoc:
    def flags : MatrixFlags
      MatrixFlags::None
    end

    # :nodoc:
    protected def flags=(flags : MatrixFlags)
      raise "cannot set flags of submatrix"
    end

    # Creates submatrix from a given `base`, `offset` and `size`
    def initialize(@base : Matrix(T), @offset, size : RowColumn)
      raise IndexError.new("submatrix offset can't be negative") if @offset.any? &.<(0)
      if @offset[0] + size[0] > @base.nrows || @offset[1] + size[1] > @base.ncolumns
        raise IndexError.new("submatrix size exceeds matrix size")
      end
      @nrows = size[0]
      @ncolumns = size[1]
    end

    # :nodoc:
    def unsafe_set(x, y, value)
      @base.unsafe_set(@offset[0] + x, @offset[1] + y, value)
    end

    # :nodoc:
    def unsafe_fetch(x, y)
      @base.unsafe_fetch(@offset[0] + x, @offset[1] + y)
    end

    # Returns submatrix pointing to same area
    def dup
      SubMatrix(T).new(@base, @offset, {@nrows, @ncolumns})
    end

    # Converts a submatrix to a GeneralMatrix(T)
    def clone
      to_general
    end

    # :nodoc:
    def transpose!
      # TODO - maybe it is possible for square?
      raise "impossible for submatrix"
    end

    # :nodoc:
    def conjtranspose!
      # TODO - maybe it is possible for square?
      raise "impossible for submatrix"
    end
  end
end
