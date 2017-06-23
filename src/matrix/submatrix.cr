require "./matrix"

module Linalg
  alias RowColumn = {Int32, Int32}

  # it's like Slice, but for matrices.
  # in future will be improved to provide same interface as matrix
  struct SubMatrix(T)
    include Matrix(T)
    getter offset
    getter size

    def flags
      MatrixFlags::Virtual
    end

    def initialize(@base : Matrix(T), @offset : RowColumn, @size : RowColumn)
      raise IndexError.new("submatrix offset can't be negative") if @offset.any? &.<(0)
      if @offset[0] + @size[0] > @base.rows || @offset[1] + @size[1] > @base.columns
        raise IndexError.new("submatrix size exceeds matrix size")
      end
    end

    def rows
      @size[0]
    end

    def columns
      @size[1]
    end

    def unsafe_set(x, y, value)
      @base.unsafe_set(@offset[0] + x, @offset[1] + y, value)
    end

    def unsafe_at(x, y)
      @base.unsafe_at(@offset[0] + x, @offset[1] + y)
    end
  end
end
