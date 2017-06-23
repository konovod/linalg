require "complex"

struct Complex
  def self.new(value)
    case value
    when Complex
      new(value.real, value.imag)
    else
      new(value, 0.0)
    end
  end
end

module Linalg
  # TODO - Complex64?
  SUPPORTED_TYPES = {Float32, Float64, Complex}

  @[Flags]
  enum MatrixFlags
    Symmetric
    Hermitian
    PositiveDefinite
    Hessenberg
    Band
    Diagonal
    Bidiagonal
    Tridiagonal
    Triangular
    Orthogonal

    Virtual # compute values instead of storing
  end

  # class that provide all utility matrix functions
  # TODO - iteration on cols\rows (create row\column object to prevent allocations?)
  # TODO - sums on cols\rows, check numpy for more (require previous point?)
  # TODO - saving/loading to files (what formats? csv?)
  module AbstractMatrix(T)
    # used in constructors to limit T at compile-time
    protected def check_type
      {% unless T == Float32 || T == Float64 || T == Complex %}
          {% raise "Wrong matrix members type: #{T}. Types supported by Linalg are: #{SUPPORTED_TYPES}" %}
        {% end %}
    end

    # creates actual matrix with same content. Useful for virtual matrices
    def clone
      Matrix(T).new(rows, columns) do |i, j|
        unsafe_at(i, j)
      end
    end

    # matrix product to given m
    def *(m : self)
      if columns != m.rows
        raise ArgumentError.new("matrix size should match ([#{rows}x#{columns}] * [#{m.rows}x#{m.columns}]")
      end
      result = Matrix(T).new(rows, m.columns) do |i, j|
        (0...columns).sum { |k| self[i, k]*m[k, j] }
      end
    end

    # multiplies at scalar
    def *(k : Int | Float | Complex)
      result = Matrix(T).new(rows, columns) do |i, j|
        self[i, j]*k
      end
    end

    # returns element-wise sum
    def +(m : self)
      if columns != m.columns || rows != m.rows
        raise ArgumentError.new("matrix size should match ([#{rows}x#{columns}] + [#{m.rows}x#{m.columns}]")
      end
      result = Matrix(T).new(rows, columns) do |i, j|
        self[i, j] + m[i, j]
      end
    end

    # returns element-wise subtract
    def -(m : self)
      if columns != m.columns || rows != m.rows
        raise ArgumentError.new("matrix size should match ([#{rows}x#{columns}] - [#{m.rows}x#{m.columns}]")
      end
      result = Matrix(T).new(rows, columns) do |i, j|
        self[i, j] - m[i, j]
      end
    end

    # returns matrix norm
    # TODO - proper norms
    def abs
      (0...rows).map { |r| (0...columns).sum { |c| self[r, c].abs } }.max
    end

    # returns transposed matrix
    def transpose
      Matrix(T).new(columns, rows) do |i, j|
        self[j, i]
      end
    end

    # retunrs
    def kron(b : self)
      Matrix(T).kron(self, b)
    end

    def tril(k = 0)
      Matrix(T).new(rows, columns) do |i, j|
        i >= j - k ? self[i, j] : 0
      end
    end

    def triu(k = 0)
      Matrix(T).new(rows, columns) do |i, j|
        i <= j - k ? self[i, j] : 0
      end
    end

    def to_s(io)
      io << "\n"
      rows.times do |i|
        io << "["
        columns.times do |j|
          io << ", " unless j == 0
          io << self[i, j]
        end
        io << "]\n"
      end
      io << "\n"
    end

    def reshape(arows, acolumns)
      clone.reshape!(arows, acolumns)
    end

    def square?
      rows == columns
    end

    def repmat(arows, acolumns)
      Matrix(T).new(rows*arows, columns*acolumns) do |i, j|
        self[i % rows, j % columns]
      end
    end

    def [](rows : Range(Int32, Int32), columns : Range(Int32, Int32))
      nrows = rows.end + (rows.excludes_end? ? 0 : 1) - rows.begin
      ncols = columns.end + (columns.excludes_end? ? 0 : 1) - columns.begin
      SubMatrix(T).new(self, {rows.begin, columns.begin}, {nrows, ncols})
    end

    def row(i)
      SubMatrix(T).new(self, {i, 0}, {1, columns})
    end

    def column(i)
      SubMatrix(T).new(self, {0, i}, {rows, 1})
    end

    def [](i, j)
      if j >= 0 && j < columns && i >= 0 && i < rows
        unsafe_at(i, j)
      else
        raise IndexError.new("access to [#{i}, #{j}] in matrix with size #{rows}x#{columns}")
      end
    end

    def []=(i, j, value)
      # i isn't checked as underlying array will check it anyway
      if j >= 0 && j < columns && i >= 0 && i < rows
        unsafe_set(i, j, value)
      else
        raise IndexError.new("access to [#{i}, #{j}] in matrix with size #{rows}x#{columns}")
      end
    end
  end

  # generic matrix, heap-allocated
  # TODO - constructing from matlab-like string "[1,2,3;3,4,6;1,1,3]" (check regexps?)
  class Matrix(T)
    include AbstractMatrix(T)
    getter rows : Int32
    getter columns : Int32
    getter raw : Slice(T)
    property flags : MatrixFlags = MatrixFlags.new(0)

    private def check_type
      {% unless T == Float32 || T == Float64 || T == Complex %}
        {% raise "Wrong matrix members type: #{T}. Types supported by Linalg are: #{SUPPORTED_TYPES}" %}
      {% end %}
    end

    def initialize(@rows, @columns)
      check_type
      @raw = Slice(T).new(rows*columns, T.new(0))
    end

    def initialize(@rows, @columns, values, @flags = MatrixFlags.new(0))
      check_type
      @raw = Slice(T).new(rows*columns) { |i| T.new(values[i]) }
    end

    def initialize(values)
      check_type
      @rows = values.size
      @columns = values[0].size
      @raw = Slice(T).new(rows*columns) do |index|
        i = index / @columns
        j = index % @columns
        raise IndexError.new("All rows should have same size") if j == 0 && values[i].size != @columns
        T.new(values[i][j])
      end
    end

    def self.from(matrix)
      new(matrix.rows, matrix.columns, matrix.raw)
    end

    def initialize(@rows, @columns, &block)
      check_type
      @raw = Slice(T).new(@rows*@columns) do |index|
        i = index / @columns
        j = index % @columns
        T.new(yield(i, j))
      end
    end

    def unsafe_at(i, j)
      @raw.unsafe_at(i*columns + j)
    end

    def unsafe_set(i, j, value)
      @raw[i*columns + j] = T.new(value)
    end

    # TODO - benchmark is it faster?
    def dup
      Matrix(T).new(@rows, @columns, @raw)
    end

    def clone
      dup
    end

    def to_unsafe
      {% if T == Complex %}
        @raw.to_unsafe.as(LibLAPACKE::DoubleComplex*)
      {% else %}
        @raw.to_unsafe
      {% end %}
    end

    def ==(other : self)
      @rows == other.rows && @columns == other.columns && @raw == other.raw
    end

    # transposes matrix inplace
    def transpose!
      if square?
        (0..@rows - 2).each do |i|
          (i + 1..@columns - 1).each do |j|
            a = self[i, j]
            self[i, j] = self[j, i]
            self[j, i] = a
          end
        end
      else
        # TODO https://en.wikipedia.org/wiki/In-place_matrix_transposition
        raise "not implemented yet"
      end
    end

    def reshape!(arows, acolumns)
      raise ArgumentError.new("number of elements should not change") if arows*acolumns != @raw.size
      @rows = arows
      @columns = acolumns
      self
    end

    def self.rand(rows, columns, rng = Random::DEFAULT)
      new(rows, columns) { |i, j| rng.rand }
    end

    def self.zeros(rows, columns)
      new(rows, columns)
    end

    def self.ones(rows, columns)
      new(rows, columns) { |i, j| 1 }
    end

    def self.repmat(a : self, rows, columns)
      a.repmat(rows, columns)
    end

    def self.diag(arows, acolumns, value : Int | Float | Complex)
      diag(arows, acolumns) { value }
    end

    def self.diag(arows, acolumns, values)
      new(arows, acolumns) do |i, j|
        i == j ? values[i] : 0
      end
    end

    def self.diag(values)
      diag(values.size, values.size, values)
    end

    def self.diag(arows, acolumns, &block)
      new(arows, acolumns) do |i, j|
        i == j ? yield(i) : 0
      end
    end

    def self.kron(a, b)
      new(a.rows*b.rows, a.columns*b.columns) do |i, j|
        a[i / b.rows, j / b.columns] * b[i % b.rows, j % b.columns]
      end
    end

    def self.tri(rows, columns, k = 0)
      new(rows, columns) do |i, j|
        i >= j - k ? 1 : 0
      end
    end

    def self.identity(n)
      new(n, n) { |i, j| i == j ? 1 : 0 }
    end

    def to_a
      @raw.to_a
    end

    def to_aa
      Array(Array(T)).new(@rows) do |i|
        Array(T).build(@columns) do |pointer|
          pointer.to_slice(@columns).copy_from(@raw[i*@columns, @columns])
          @columns
        end
      end
    end
  end

  alias RowColumn = {Int32, Int32}

  # it's like Slice, but for matrices.
  # in future will be improved to provide same interface as matrix
  struct SubMatrix(T)
    include AbstractMatrix(T)
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

  alias Mat = Matrix(Float64)
  alias Mat32 = Matrix(Float32)
  alias MatComplex = Matrix(Complex)
end
