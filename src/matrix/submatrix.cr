require "./matrix"

module LA
  alias RowColumn = {Int32, Int32}

  # it's like Slice, but for matrices.
  # made class not struct to avoid compilation issues
  class SubMatrix(T) < Matrix(T)
    getter offset
    getter nrows : Int32
    getter ncolumns : Int32

    def flags : MatrixFlags
      MatrixFlags::None
    end

    protected def flags=(flags : MatrixFlags)
      raise "cannot set flags of submatrix"
    end

    def initialize(@base : Matrix(T), @offset : RowColumn, size : RowColumn)
      raise IndexError.new("submatrix offset can't be negative") if @offset.any? &.<(0)
      if @offset[0] + size[0] > @base.nrows || @offset[1] + size[1] > @base.ncolumns
        raise IndexError.new("submatrix size exceeds matrix size")
      end
      @nrows = size[0]
      @ncolumns = size[1]
    end

    def unsafe_set(x, y, value)
      @base.unsafe_set(@offset[0] + x, @offset[1] + y, value)
    end

    def unsafe_fetch(x, y)
      @base.unsafe_fetch(@offset[0] + x, @offset[1] + y)
    end

    def dup
      SubMatrix(T).new(@base, @offset, {@nrows, @ncolumns})
    end

    def clone
      to_general
    end

    def transpose!
      raise "impossible for submatrix"
    end

    def conjtranspose!
      raise "impossible for submatrix"
    end
  end
end
