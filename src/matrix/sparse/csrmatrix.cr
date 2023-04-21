require "./sparse_matrix.cr"

# TODO - inline docs

module LA::Sparse
  class CSRMatrix(T) < Matrix(T)
    protected getter raw_columns : Array(Int32)
    protected getter raw_rows : Array(Int32)
    protected getter raw_values : Array(T)

    def initialize(@nrows, @ncolumns, capacity = 0)
      @raw_rows = Array(Int32).new(nrows + 1)
      @raw_columns = Array(Int32).new(capacity)
      @raw_values = Array(T).new(capacity)
      @flags = MatrixFlags.for_diag(@nrows == @ncolumns)
    end

    def initialize(@nrows, @ncolumns, raw_rows : Array(Int32), raw_columns : Array(Int32), raw_values : Array(T), @flags = MatrixFlags::None, *, dont_clone : Bool = false)
      if raw_rows.size != @nrows + 1
        raise ArgumentError.new("Can't construct CSR Matrix from arrays of different size: rows.size(#{raw_rows.size}) != nrows+1 (#{@nrows + 1}")
      end
      if raw_columns.size != raw_values.size
        raise ArgumentError.new("Can't construct CSR Matrix from arrays of different size: columns.size(#{raw_columns.size}) != values.size(#{values.size})")
      end
      if dont_clone
        @raw_rows = rows
        @raw_columns = columns
        @raw_values = values
      else
        @raw_rows = rows.dup
        @raw_columns = columns.dup
        @raw_values = values.dup
      end
    end

    def nonzeros : Int32
      @raw_values.size
    end

    def self.new(matrix : CSRMatrix(T))
      new(matrix.nrows, matrix.ncolumns, matrix.raw_rows, matrix.raw_columns, matrix.raw_values, flags: matrix.flags)
    end

    def self.new(matrix : CSRMatrix)
      new(matrix.nrows, matrix.ncolumns, matrix.raw_rows.dup, matrix.raw_columns.dup, matrix.raw_values.map { |v| T.new(v) }, dont_clone: true, flags: matrix.flags)
    end

    # def self.new(matrix : LA::Matrix) TODO

    private def ij2index(i, j) : Int32?
      row_start = @raw_rows[i]
      row_end = @raw_rows[i + 1]
      index = (row_start...row_end).bsearch { |n| @raw_columns[n] >= j }
      @raw_columns[n] == j ? index : nil
    end

    private def index2ij(index) : {Int32, Int32}?
      # TODO
      return nil
    end

    def unsafe_fetch(i, j)
      @raw_values[ij2index(i, j)]
    end

    # def unsafe_set(i, j, value)
    # def self.diag(nrows, ncolumns, values) TODO
    # def self.diag(nrows, ncolumns, &block) TODO
    def each_index(*, all = false, &block)
      @nrows.times do |i|
        row_start = @raw_rows[i]
        row_end = @raw_rows[i + 1]
        (row_start...row_end).each do |index|
          yield(i, @row_columns[index])
        end
      end
    end

    def each_with_index(*, all = false, &block)
      @nrows.times do |i|
        row_start = @raw_rows[i]
        row_end = @raw_rows[i + 1]
        (row_start...row_end).each do |index|
          yield(@row_values[index], i, @row_columns[index])
        end
      end
    end

    def map_with_index!(&block)
      @nrows.times do |i|
        row_start = @raw_rows[i]
        row_end = @raw_rows[i + 1]
        (row_start...row_end).each do |index|
          @row_values[index] = yield(@row_values[index], i, @row_columns[index])
        end
      end
    end

    def map_with_index(&block)
      values = @raw_values.map_with_index { |v, i| T.new(yield(v, @raw_rows[i], @raw_columns[i])) }
      CSRMatrix(T).new(@nrows, @ncolumns, @raw_rows.dup, @raw_columns.dup, values, dont_clone: true)
    end

    def map_with_index_f64(&block)
      values = @raw_values.map_with_index { |v, i| Float64.new(yield(v, @raw_rows[i], @raw_columns[i])) }
      CSRMatrix(Float64).new(@nrows, @ncolumns, @raw_rows.dup, @raw_columns.dup, values, dont_clone: true)
    end

    def map_with_index_complex(&block)
      values = @raw_values.map_with_index { |v, i| Complex.new(yield(v, @raw_rows[i], @raw_columns[i])) }
      CSRMatrix(Complex).new(@nrows, @ncolumns, @raw_rows.dup, @raw_columns.dup, values, dont_clone: true)
    end

    # def transpose! TODO
    # def transpose TODO
    # def conjtranspose TODO
    # def add!(k : Number | Complex, m : Sparse::Matrix) TODO
    # def clear TODO
    # def triu!(k = 0) TODO
    # def tril!(k = 0) TODO
    # def tril(k = 0) TODO
    # def triu(k = 0) TODO
    # def self.rand(nrows, ncolumns, *, nonzero_elements, rng : Random = Random::DEFAULT)
    # def self.rand(nrows, ncolumns, *, fill_factor, rng : Random = Random::DEFAULT)
    # def select!(& : T -> Bool)
    #   def select_with_index!(& : (T, Int32, Int32) -> Bool)
    #   def select_index!(& : (Int32, Int32) -> Bool)
    # def resize!(anrows, ancolumns)

  end
end
