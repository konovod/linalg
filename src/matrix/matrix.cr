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
    Orthogonal
    UpperTriangular
    LowerTriangular
    Triangular      = UpperTriangular | LowerTriangular

    def triangular?
      self.upper_triangular? || self.lower_triangular?
    end

    def diagonal?
      self.upper_triangular? && self.lower_triangular?
    end

    def sum(f2 : MatrixFlags)
      self & f2 & (MatrixFlags::UpperTriangular |
                   MatrixFlags::LowerTriangular |
                   MatrixFlags::Symmetric |
                   MatrixFlags::Hermitian)
    end

    def mult(f2 : MatrixFlags)
      self & f2 & (MatrixFlags::UpperTriangular |
                   MatrixFlags::LowerTriangular |
                   MatrixFlags::Symmetric)
    end

    def transpose
      result = self
      if triangular?
        wasup, waslo = upper_triangular?, lower_triangular?
        if wasup
          result |= MatrixFlags::LowerTriangular
        else
          result &= ~MatrixFlags::LowerTriangular
        end
        if waslo
          result |= MatrixFlags::UpperTriangular
        else
          result &= ~MatrixFlags::UpperTriangular
        end
      end
      result
    end

    def scale
      self
    end

    def self.for_diag(square : Bool)
      if square
        Symmetric | UpperTriangular | LowerTriangular | Hermitian
      else
        UpperTriangular | LowerTriangular
      end
    end
  end

  # class that provide all utility matrix functions
  # TODO - sums on cols\nrows, check numpy for more (require previous point?)
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
      {nrows, ncolumns}
    end

    # creates generic matrix with same content. Useful for virtual matrices
    def clone
      GeneralMatrix(T).new(nrows, ncolumns) do |i, j|
        unsafe_at(i, j)
      end.tap { |it| it.flags = flags }
    end

    # matrix product to given m
    def *(m : Matrix(T))
      if ncolumns != m.nrows
        raise ArgumentError.new("matrix size should match ([#{nrows}x#{ncolumns}] * [#{m.nrows}x#{m.ncolumns}]")
      end
      result = GeneralMatrix(T).new(nrows, m.ncolumns) do |i, j|
        (0...ncolumns).sum { |k| self[i, k]*m[k, j] }
      end
      result.tap { |r| r.flags = self.flags.mult(m.flags) }
    end

    # multiplies at scalar
    def *(k : Number | Complex)
      result = GeneralMatrix(T).new(nrows, ncolumns) do |i, j|
        self[i, j]*k
      end
      result.tap { |r| r.flags = self.flags.scale }
    end

    def -
      result = self*(-1)
    end

    # divides at scalar
    def /(k : Number | Complex)
      result = GeneralMatrix(T).new(nrows, ncolumns) do |i, j|
        self[i, j] / k
      end
      result.tap { |r| r.flags = self.flags.scale }
    end

    # returns element-wise sum
    def +(m : Matrix(T))
      if ncolumns != m.ncolumns || nrows != m.nrows
        raise ArgumentError.new("matrix size should match ([#{nrows}x#{ncolumns}] + [#{m.nrows}x#{m.ncolumns}]")
      end
      result = GeneralMatrix(T).new(nrows, ncolumns) do |i, j|
        self[i, j] + m[i, j]
      end
      result.tap { |r| r.flags = self.flags.sum(m.flags) }
    end

    # returns element-wise subtract
    def -(m : Matrix(T))
      if ncolumns != m.ncolumns || nrows != m.nrows
        raise ArgumentError.new("matrix size should match ([#{nrows}x#{ncolumns}] - [#{m.nrows}x#{m.ncolumns}]")
      end
      result = GeneralMatrix(T).new(nrows, ncolumns) do |i, j|
        self[i, j] - m[i, j]
      end
      result.tap { |r| r.flags = self.flags.sum(m.flags) }
    end

    # returns transposed matrix
    def transpose
      return clone if flags.symmetric?
      GeneralMatrix(T).new(ncolumns, nrows) do |i, j|
        self[j, i]
      end.tap do |m|
        m.flags = flags.transpose
      end
    end

    # returns transposed matrix
    def conjtranspose
      {% if T != Complex %}
        return transpose
      {% end %}
      return clone if flags.hermitian?
      GeneralMatrix(T).new(ncolumns, nrows) do |i, j|
        self[j, i].conj
      end.tap do |m|
        m.flags = flags.transpose
      end
    end

    # returns kroneker product with matrix b
    def kron(b : Matrix(T))
      Matrix(T).kron(self, b)
    end

    # same as tril in scipy - returns lower triangular or trapezoidal part of matrix
    def tril(k = 0)
      x = GeneralMatrix(T).new(nrows, ncolumns) do |i, j|
        i >= j - k ? self[i, j] : 0
      end
      if k >= 0
        x.assume! MatrixFlags::LowerTriangular
      end
      x
    end

    # same as triu in scipy - returns upper triangular or trapezoidal part of matrix
    def triu(k = 0)
      x = GeneralMatrix(T).new(nrows, ncolumns) do |i, j|
        i <= j - k ? self[i, j] : 0
      end
      if k >= 0
        x.assume! MatrixFlags::UpperTriangular
      end
      x
    end

    # converts to string, with linefeeds before and after matrix:
    # [1, 2, 3, .... 10]
    # [11, 12, 13, .... 20]
    # ...
    # [91, 92, 93, .... 100]
    def to_s(io)
      io << "\n"
      nrows.times do |i|
        io << "["
        ncolumns.times do |j|
          io << ", " unless j == 0
          io << self[i, j]
        end
        io << "]\n"
      end
      io << "\n"
    end

    def inspect(io)
      io << self.class << " (" << nrows << "x" << ncolumns << ", " << flags << "):"
      to_s(io)
    end

    def each(&block)
      nrows.times do |row|
        ncolumns.times do |column|
          yield self.unsafe_at(row, column)
        end
      end
    end

    def each_index(&block)
      nrows.times do |row|
        ncolumns.times do |column|
          yield row, column
        end
      end
    end

    def each_with_index(&block)
      nrows.times do |row|
        ncolumns.times do |column|
          yield self.unsafe_at(row, column), row, column
        end
      end
    end

    def ==(other)
      return false unless nrows == other.nrows && ncolumns == other.ncolumns
      each_with_index do |value, row, column|
        return false if other.unsafe_at(row, column) != value
      end
      true
    end

    # changes nrows and ncolumns of matrix (total number of elements must not change)
    def reshape(anrows, ancolumns)
      clone.reshape!(anrows, ancolumns)
    end

    # returns True if matrix is square and False otherwise
    def square?
      nrows == ncolumns
    end

    # return matrix repeated `arows` times by vertical and `acolumns` times by horizontal
    def repmat(arows, acolumns)
      GeneralMatrix(T).new(nrows*arows, ncolumns*acolumns) do |i, j|
        self[i % nrows, j % ncolumns]
      end
    end

    def [](i, j)
      if j >= 0 && j < ncolumns && i >= 0 && i < nrows
        unsafe_at(i, j)
      else
        raise IndexError.new("access to [#{i}, #{j}] in matrix with size #{nrows}x#{ncolumns}")
      end
    end

    def []=(i, j, value)
      if j >= 0 && j < ncolumns && i >= 0 && i < nrows
        unsafe_set(i, j, value)
      else
        raise IndexError.new("access to [#{i}, #{j}] in matrix with size #{nrows}x#{ncolumns}")
      end
    end

    # return submatrix over given ranges.
    def [](arows : Range(Int32, Int32), acolumns : Range(Int32, Int32))
      start_row = arows.begin + (arows.begin < 0 ? @nrows : 0)
      start_col = acolumns.begin + (acolumns.begin < 0 ? @ncolumns : 0)
      total_rows = arows.end + (arows.excludes_end? ? 0 : 1) - start_row + (arows.end < 0 ? @nrows : 0)
      total_cols = acolumns.end + (acolumns.excludes_end? ? 0 : 1) - start_col + (acolumns.end < 0 ? @ncolumns : 0)
      SubMatrix(T).new(self, {start_row, start_col}, {total_rows, total_cols})
    end

    def [](row : Int32, acolumns : Range(Int32, Int32))
      self[row..row, acolumns]
    end

    def [](arows : Range(Int32, Int32), column : Int32)
      self[arows, column..column]
    end

    def []=(arows : Range(Int32, Int32), acolumns : Range(Int32, Int32), value)
      submatrix = self[arows, acolumns]
      if value.is_a? Matrix
        raise IndexError.new("submatrix size must match assigned value") unless submatrix.size == value.size
        submatrix.each_index { |i, j| submatrix.unsafe_set i, j, value.unsafe_at(i, j) }
      else
        submatrix.each_index { |i, j| submatrix.unsafe_set i, j, value }
      end
    end

    def []=(row : Int32, acolumns : Range(Int32, Int32), value)
      self[row..row, acolumns] = value
    end

    def []=(nrows : Range(Int32, Int32), column : Int32, value)
      self[nrows, column..column] = value
    end

    def assume!(flag : MatrixFlags, value : Bool = true)
      if value
        @flags |= flag
      else
        @flags &= ~flag
      end
    end

    def tolerance
      amax = 0.0
      # TODO - better iteration?
      each do |value|
        vabs = {% if T == Complex %}value.real.abs + value.imag.abs{% else %}value.abs{% end %}
        amax = vabs if amax < vabs
      end
      # TODO - better eps?
      eps = {% if T == Float32 %}2e-7{% else %}2e-16{% end %}
      amax * nrows * ncolumns * 10*eps
    end

    def almost_eq(other : Matrix(T), eps)
      each_with_index do |value, row, column|
        return false if (value - other.unsafe_at(row, column)).abs > eps
      end
      true
    end

    def almost_eq(other : Matrix(T))
      almost_eq other, {self.tolerance, other.tolerance}.min
    end

    # TODO - check for all flags
    private def check_single(flag : MatrixFlags, eps = tolerance)
      case flag
      when .symmetric?
        # TODO - 2 times less work
        return false unless square?
        each_with_index do |value, row, column|
          return false if row < column && (value - unsafe_at(row, column)).abs > eps
        end
        return true
      when .hermitian?
        {% if T == Complex %}
          return false unless square?
          each_with_index do |value, row, column|
            return false if row < column && (value.conj - unsafe_at(row, column)).abs > eps
          end
          return true
        {% else %}
          return check_single(MatrixFlags::Symmetric, eps)
        {% end %}
      when .positive_definite?
        # TODO - cleaner detection?
        return false unless square? && check_single(MatrixFlags::Hermitian, eps)
        begin
          cholesky(dont_clean: true)
          return true
        rescue LinAlgError
          return false
        end
      when .orthogonal?
        return square? && (self*self.conjt).almost_eq Matrix(T).identity(nrows), eps
      when .upper_triangular?
        each_with_index do |value, row, column|
          return false if row > column && (value).abs > eps
        end
        return true
      when .lower_triangular?
        each_with_index do |value, row, column|
          return false if row < column && (value).abs > eps
        end
        return true
      else
        return false
      end
      return false
    end

    private def detect_single(flag : MatrixFlags, eps = tolerance)
      check_single(flag, eps).tap do |ok|
        assume!(flag, ok)
      end
    end

    def detect(aflags : MatrixFlags = MatrixFlags::All, eps = tolerance)
      result = true
      {MatrixFlags::Symmetric,
       MatrixFlags::Hermitian,
       MatrixFlags::PositiveDefinite,
       MatrixFlags::Orthogonal,
       MatrixFlags::LowerTriangular,
       MatrixFlags::UpperTriangular}.each do |f|
        if aflags & f != MatrixFlags::None
          result = false unless detect_single(f, eps)
        end
      end
      if aflags.triangular?
        result = false unless detect_single(MatrixFlags::LowerTriangular, eps) || detect_single(MatrixFlags::UpperTriangular, eps)
      end
      result
    end

    def assume(aflags : MatrixFlags, value : Bool = true)
      if value
        detect(aflags)
      else
        assume! aflags, false
      end
    end

    def clear_flags
      @flags = MatrixFlags::None
    end

    def self.rand(nrows, ncolumns, rng = Random::DEFAULT)
      GeneralMatrix(T).new(nrows, ncolumns) { |i, j| rng.rand }
    end

    def self.zeros(nrows, ncolumns)
      GeneralMatrix(T).new(nrows, ncolumns).tap { |m| m.flags = MatrixFlags.for_diag(m.square?) }
    end

    def self.ones(nrows, ncolumns)
      GeneralMatrix(T).new(nrows, ncolumns) { |i, j| 1 }.tap do |m|
        if m.square?
          m.flags = MatrixFlags::Symmetric | MatrixFlags::Hermitian
        end
      end
    end

    def self.repmat(a : Matrix(T), nrows, ncolumns)
      a.repmat(nrows, ncolumns)
    end

    def self.diag(nrows, ncolumns, value : Number | Complex)
      diag(nrows, ncolumns) { value }
    end

    def self.diag(nrows, ncolumns, values)
      GeneralMatrix(T).new(nrows, ncolumns) do |i, j|
        i == j ? values[i] : 0
      end.tap { |m| m.flags = MatrixFlags.for_diag(m.square?) }
    end

    def self.diag(values)
      diag(values.size, values.size, values)
    end

    def self.column(values)
      GeneralMatrix(T).new(values.size, 1, values)
    end

    def self.row(values)
      GeneralMatrix(T).new(1, values.size, values)
    end

    def self.diag(nrows, ncolumns, &block)
      GeneralMatrix(T).new(nrows, ncolumns) do |i, j|
        i == j ? yield(i) : 0
      end.tap { |m| m.flags = MatrixFlags.for_diag(m.square?) }
    end

    def self.kron(a, b)
      GeneralMatrix(T).new(a.nrows*b.nrows, a.ncolumns*b.ncolumns) do |i, j|
        a[i / b.nrows, j / b.ncolumns] * b[i % b.nrows, j % b.ncolumns]
      end
    end

    def self.identity(n)
      result = GeneralMatrix(T).new(n, n) { |i, j| i == j ? 1 : 0 }
      result.tap { |m| m.flags = MatrixFlags.for_diag(true) | MatrixFlags::PositiveDefinite }
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

    def cat(other : Matrix(T), dimension)
      raise ArgumentError.new("only dimesion = 0 or 1 supported") unless {0, 1}.includes? dimension
      if self.size[1 - dimension] != other.size[1 - dimension]
        raise ArgumentError.new("matrix size along other dimension should match for concatenation")
      end
      if dimension == 0
        GeneralMatrix(T).new(nrows + other.nrows, ncolumns) do |row, column|
          row < nrows ? unsafe_at(row, column) : other.unsafe_at(row - nrows, column)
        end
      else
        GeneralMatrix(T).new(nrows, ncolumns + other.ncolumns) do |row, column|
          column < ncolumns ? unsafe_at(row, column) : other.unsafe_at(row, column - ncolumns)
        end
      end
    end

    def vcat(other)
      cat other, 0
    end

    def hcat(other)
      cat other, 1
    end

    def map(&block)
      GeneralMatrix(T).new(nrows, ncolumns) { |i, j| yield(unsafe_at(i, j)) }
    end

    def map_with_index(&block)
      GeneralMatrix(T).new(nrows, ncolumns) { |i, j| yield(unsafe_at(i, j), i, j) }
    end
  end

  alias Mat = Matrix(Float64)
  alias Mat32 = Matrix(Float32)
  alias MatComplex = Matrix(Complex)
end
