require "complex"

module Linalg
  # TODO - Complex64?
  SUPPORTED_TYPES = {Float32, Float64, Complex}

  @[Flags]
  enum MatrixFlags
    Symmetric
    Hermitian
    PositiveDefinite
    # Hessenberg
    # Band
    # Diagonal
    # Bidiagonal
    # Tridiagonal
    Triangular
    Orthogonal
    Unitary    = Orthogonal
    Lower
  end

  # class that provide all utility matrix functions
  # TODO - sums on cols\rows, check numpy for more (require previous point?)
  # TODO - saving/loading to files (what formats? csv?)
  # TODO - replace [] to unsafe at most places
  module Matrix(T)
    # used in constructors to limit T at compile-time
    protected def check_type
      {% unless T == Float32 || T == Float64 || T == Complex %}
        {% raise "Wrong matrix members type: #{T}. Types supported by Linalg are: #{SUPPORTED_TYPES}" %}
      {% end %}
    end

    # to_unsafe method raises at runtime and is overriden by matrix that actually have pointer
    def to_unsafe
      raise ArgumentError.new("Virtual matrix can't be passed unsafe!")
    end

    def size
      {rows, columns}
    end

    # creates generic matrix with same content. Useful for virtual matrices
    def clone
      GeneralMatrix(T).new(rows, columns) do |i, j|
        unsafe_at(i, j)
      end.tap { |it| it.flags = flags }
    end

    # matrix product to given m
    def *(m : Matrix(T))
      if columns != m.rows
        raise ArgumentError.new("matrix size should match ([#{rows}x#{columns}] * [#{m.rows}x#{m.columns}]")
      end
      result = GeneralMatrix(T).new(rows, m.columns) do |i, j|
        (0...columns).sum { |k| self[i, k]*m[k, j] }
      end
    end

    # multiplies at scalar
    def *(k : Number | Complex)
      result = GeneralMatrix(T).new(rows, columns) do |i, j|
        self[i, j]*k
      end
    end

    # divides at scalar
    def /(k : Number | Complex)
      result = GeneralMatrix(T).new(rows, columns) do |i, j|
        self[i, j] / k
      end
    end

    # returns element-wise sum
    def +(m : Matrix(T))
      if columns != m.columns || rows != m.rows
        raise ArgumentError.new("matrix size should match ([#{rows}x#{columns}] + [#{m.rows}x#{m.columns}]")
      end
      result = GeneralMatrix(T).new(rows, columns) do |i, j|
        self[i, j] + m[i, j]
      end
    end

    # returns element-wise subtract
    def -(m : Matrix(T))
      if columns != m.columns || rows != m.rows
        raise ArgumentError.new("matrix size should match ([#{rows}x#{columns}] - [#{m.rows}x#{m.columns}]")
      end
      result = GeneralMatrix(T).new(rows, columns) do |i, j|
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
      return clone if flags.symmetric?
      GeneralMatrix(T).new(columns, rows) do |i, j|
        self[j, i]
      end.tap do |m|
        m.flags = flags
        if flags.triangular?
          m.flags ^= MatrixFlags::Lower
        end
      end
    end

    # returns transposed matrix
    def conjtranspose
      {% raise "Matrix must be Complex for conjtranspose" unless T == Complex %}
      return clone if flags.hermitian?
      GeneralMatrix(T).new(columns, rows) do |i, j|
        self[j, i].conj
      end.tap do |m|
        m.flags = flags
        if flags.triangular?
          m.flags ^= MatrixFlags::Lower
        end
      end
    end

    # returns kroneker product with matrix b
    def kron(b : Matrix(T))
      Matrix(T).kron(self, b)
    end

    # same as tril in scipy - returns lower triangular or trapezoidal part of matrix
    def tril(k = 0)
      x = GeneralMatrix(T).new(rows, columns) do |i, j|
        i >= j - k ? self[i, j] : 0
      end
      if k >= 0
        x.assume! MatrixFlags::Triangular | MatrixFlags::Lower
      end
      x
    end

    # same as triu in scipy - returns upper triangular or trapezoidal part of matrix
    def triu(k = 0)
      x = GeneralMatrix(T).new(rows, columns) do |i, j|
        i <= j - k ? self[i, j] : 0
      end
      if k >= 0
        x.assume! MatrixFlags::Triangular
      end
      x
    end

    # like a tril in scipy - remove all elements above k-diagonal
    def tril!(k = 0)
      (rows*columns).times do |index|
        i = index / columns
        j = index % columns
        @raw[index] = T.new(0) if i < j - k
      end
      if k <= 0
        self.assume! MatrixFlags::Triangular | MatrixFlags::Lower
      end
      self
    end

    # like a triu in scipy - remove all elements below k-diagonal
    def triu!(k = 0)
      (rows*columns).times do |index|
        i = index / columns
        j = index % columns
        @raw[index] = T.new(0) if i > j - k
      end
      if k >= 0
        self.assume! MatrixFlags::Triangular
      end
      self
    end

    # converts to string, with linefeeds before and after matrix:
    # [1, 2, 3, .... 10]
    # [11, 12, 13, .... 20]
    # ...
    # [91, 92, 93, .... 100]
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

    def each(&block)
      rows.times do |row|
        columns.times do |column|
          yield self.unsafe_at(row, column)
        end
      end
    end

    def each_index(&block)
      rows.times do |row|
        columns.times do |column|
          yield row, column
        end
      end
    end

    def each_with_index(&block)
      rows.times do |row|
        columns.times do |column|
          yield self.unsafe_at(row, column), row, column
        end
      end
    end

    def ==(other)
      return false unless rows == other.rows && columns == other.columns
      each_with_index do |value, row, column|
        return false if other.unsafe_at(row, column) != value
      end
      true
    end

    # changes number of rows and columns of matrix (total number of elements must not change)
    def reshape(arows, acolumns)
      clone.reshape!(arows, acolumns)
    end

    # returns True if matrix is square and False otherwise
    def square?
      rows == columns
    end

    # return matrix repeated `arows` times by vertical and `acolumns` times by horizontal
    def repmat(arows, acolumns)
      GeneralMatrix(T).new(rows*arows, columns*acolumns) do |i, j|
        self[i % rows, j % columns]
      end
    end

    # return submatrix over given ranges.
    def [](rows : Range(Int32, Int32), columns : Range(Int32, Int32))
      nrows = rows.end + (rows.excludes_end? ? 0 : 1) - rows.begin
      ncols = columns.end + (columns.excludes_end? ? 0 : 1) - columns.begin
      SubMatrix(T).new(self, {rows.begin, columns.begin}, {nrows, ncols})
    end

    def [](row : Int32, columns : Range(Int32, Int32))
      self[row..row, columns]
    end

    def [](rows : Range(Int32, Int32), column : Int32)
      self[rows, column..column]
    end

    def []=(rows : Range(Int32, Int32), columns : Range(Int32, Int32), value)
      if value.is_a? Matrix
        rows.each do |row|
          columns.each do |column|
            self[row, column] = value[row - rows.begin, column - columns.begin]
          end
        end
      else
        rows.each do |row|
          columns.each do |column|
            self[row, column] = value
          end
        end
      end
    end

    def []=(row : Int32, columns : Range(Int32, Int32), value)
      self[row..row, columns] = value
    end

    def []=(rows : Range(Int32, Int32), column : Int32, value)
      self[rows, column..column] = value
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
      if j >= 0 && j < columns && i >= 0 && i < rows
        unsafe_set(i, j, value)
      else
        raise IndexError.new("access to [#{i}, #{j}] in matrix with size #{rows}x#{columns}")
      end
    end

    def assume!(flag : MatrixFlags)
      @flags |= flag
    end

    # TODO - check for all flags
    # def detect(flag : MatrixFlags)
    #   @flags |=
    # end

    def clear_flags
      @flags = MatrixFlags.new(0)
    end

    def self.rand(rows, columns, rng = Random::DEFAULT)
      GeneralMatrix(T).new(rows, columns) { |i, j| rng.rand }
    end

    def self.zeros(rows, columns)
      GeneralMatrix(T).new(rows, columns)
    end

    def self.ones(rows, columns)
      GeneralMatrix(T).new(rows, columns) { |i, j| 1 }
    end

    def self.repmat(a : Matrix(T), rows, columns)
      a.repmat(rows, columns)
    end

    def self.diag(arows, acolumns, value : Number | Complex)
      diag(arows, acolumns) { value }
    end

    def self.diag(arows, acolumns, values)
      GeneralMatrix(T).new(arows, acolumns) do |i, j|
        i == j ? values[i] : 0
      end
    end

    def self.diag(values)
      diag(values.size, values.size, values)
    end

    def self.diag(arows, acolumns, &block)
      GeneralMatrix(T).new(arows, acolumns) do |i, j|
        i == j ? yield(i) : 0
      end
    end

    def self.kron(a, b)
      GeneralMatrix(T).new(a.rows*b.rows, a.columns*b.columns) do |i, j|
        a[i / b.rows, j / b.columns] * b[i % b.rows, j % b.columns]
      end
    end

    def self.tri(rows, columns, k = 0)
      GeneralMatrix(T).new(rows, columns) do |i, j|
        i >= j - k ? 1 : 0
      end
    end

    def self.identity(n)
      GeneralMatrix(T).new(n, n) { |i, j| i == j ? 1 : 0 }
    end

    def self.eye(n)
      self.identity(n)
    end

    def t
      transpose
    end

    def t!
      transpose!
    end

    def conjt
      conjtranspose
    end

    def conjt!
      conjtranspose!
    end
  end

  alias Mat = Matrix(Float64)
  alias Mat32 = Matrix(Float32)
  alias MatComplex = Matrix(Complex)
end
