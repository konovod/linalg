require "./sparse_matrix.cr"

# TODO - inline docs

module LA::Sparse
  class CSRMatrix(T) < Matrix(T)
    protected getter raw_columns : Array(Int32)
    protected getter raw_rows : Array(Int32)
    protected getter raw_values : Array(T)

    def initialize(@nrows, @ncolumns, capacity = 0)
      @raw_rows = Array(Int32).new(nrows + 1, 0)
      @raw_columns = Array(Int32).new(capacity)
      @raw_values = Array(T).new(capacity)
      @flags = MatrixFlags.for_diag(@nrows == @ncolumns)
    end

    def initialize(@nrows, @ncolumns, raw_rows : Array(Int32), raw_columns : Array(Int32), raw_values : Array(T), @flags = MatrixFlags::None, *, dont_clone : Bool = false)
      if raw_rows.size != @nrows + 1
        raise ArgumentError.new("Can't construct CSR Matrix from arrays of different size: rows.size(#{raw_rows.size}) != nrows+1 (#{@nrows + 1}")
      end
      if raw_columns.size != raw_values.size
        raise ArgumentError.new("Can't construct CSR Matrix from arrays of different size: columns.size(#{raw_columns.size}) != values.size(#{raw_values.size})")
      end
      if dont_clone
        @raw_rows = raw_rows
        @raw_columns = raw_columns
        @raw_values = raw_values
      else
        @raw_rows = raw_rows.dup
        @raw_columns = raw_columns.dup
        @raw_values = raw_values.dup
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

    def self.new(matrix : LA::Matrix)
      if matrix.is_a? Sparse::Matrix
        nonzeros = matrix.nonzeros
      else
        nonzeros = matrix.count { |v| !v.zero? }
      end
      result = new(matrix.nrows, matrix.ncolumns, nonzeros)
      index = 0
      result.nrows.times do |row|
        result.ncolumns.times do |column|
          v = matrix.unsafe_fetch(row, column)
          next if v.zero?
          result.raw_values << T.new(v)
          result.raw_columns << column
          index += 1
        end
        result.raw_rows[row + 1] = index
      end
      result
    end

    private def ij2index(i, j) : Int32?
      row_start = @raw_rows[i]
      row_end = @raw_rows[i + 1]
      index = (row_start...row_end).bsearch { |n| @raw_columns[n] >= j }
      return nil unless index
      @raw_columns[index] == j ? index : nil
    end

    def unsafe_fetch(i, j)
      if index = ij2index(i, j)
        @raw_values.unsafe_fetch(index)
      else
        T.new(0.0)
      end
    end

    def unsafe_set(i, j, value)
      if index = ij2index(i, j)
        @raw_values.unsafe_put(index, T.new(value))
      else
        raise "cannot add values to CSR matrix"
      end
    end

    # Returns diagonal matrix of given size with diagonal elements taken from array `values`
    #
    # Raises if `values.size > {nrows, ncolumns}.min`
    def self.diag(nrows, ncolumns, values)
      raise ArgumentError.new("Too much elements (#{values.size}) for diag matrix [#{nrows}x#{ncolumns}]") if values.size > {nrows, ncolumns}.min
      n = values.size
      new(nrows, ncolumns,
        Array(Int32).new(nrows + 1) { |i| i <= n ? i : n },
        Array(Int32).new(n) { |i| i },
        Array(T).new(n) { |i| T.new(values[i]) },
        MatrixFlags.for_diag(nrows == ncolumns),
        dont_clone: true
      )
    end

    # Returns diagonal matrix of given size with diagonal elements equal to block value
    def self.diag(nrows, ncolumns, &block)
      n = {nrows, ncolumns}.min
      new(nrows, ncolumns,
        Array(Int32).new(nrows + 1) { |i| i <= n ? i : n },
        Array(Int32).new(n) { |i| i },
        Array(T).new(n) { |i| T.new(yield(i)) },
        MatrixFlags.for_diag(nrows == ncolumns),
        dont_clone: true
      )
    end

    def each_index(*, all = false, &block)
      if all
        super(all: true) { |i, j| yield(i, j) }
      else
        @nrows.times do |i|
          row_start = @raw_rows[i]
          row_end = @raw_rows[i + 1]
          (row_start...row_end).each do |index|
            yield(i, @raw_columns[index])
          end
        end
      end
    end

    def each_with_index(*, all = false, &block)
      if all
        super(all: true) { |v, i, j| yield(v, i, j) }
      else
        @nrows.times do |i|
          row_start = @raw_rows[i]
          row_end = @raw_rows[i + 1]
          (row_start...row_end).each do |index|
            yield(@raw_values[index], i, @raw_columns[index])
          end
        end
      end
    end

    def map_with_index!(&block)
      @nrows.times do |i|
        row_start = @raw_rows[i]
        row_end = @raw_rows[i + 1]
        (row_start...row_end).each do |index|
          @raw_values[index] = T.new(yield(@raw_values[index], i, @raw_columns[index]))
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

    def clear
      @raw_rows.fill(0)
      @raw_values.clear
      @raw_columns.clear
    end

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
