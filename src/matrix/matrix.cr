require "complex"

module LA
  module Utils
    # Types supported by LAPACK
    #
    #  TODO - Complex64?
    SUPPORTED_TYPES = {Float32, Float64, Complex}
  end

  enum Enums::Axis
    Columns
    Rows
  end

  # class that provide all utility matrix functions
  abstract class Matrix(T)
    include Enumerable(T)

    # Returns number of rows in matrix
    abstract def nrows : Int32
    # Returns number of columns in matrix
    abstract def ncolumns : Int32
    # Returns flags of matrix (see `MatrixFlags`)
    abstract def flags : MatrixFlags

    # Returns flags of matrix (see `MatrixFlags`)
    protected abstract def flags=(flags : MatrixFlags)

    # :nodoc:
    def self.zero
      T.new(0)
    end

    # :nodoc:
    def self.multiplicative_identity
      T.new(1.0)
    end

    # Used in constructors to limit T at compile-time
    protected def check_type
      {% unless T == Float32 || T == Float64 || T == Complex %}
        {% raise "Wrong matrix members type: #{T}. Types supported by linalg are: #{SUPPORTED_TYPES}" %}
      {% end %}
    end

    private macro of_real_type(container, size)
      {% if T == Complex %} {{container}}(Float64).new({{size}}, 0.0) {% else %} {{container}}(T).new({{size}}, T.new(0)) {% end %}
    end

    private macro of_real_type(value)
      {% if T == Complex %} Float64.new({{value}}) {% else %} T.new({{value}}) {% end %}
    end

    private macro real_type_const(constant)
      {% if T == Complex %} Float64::{{constant}} {% else %} T::{{constant}} {% end %}
    end

    # Returns pointer to underlying data
    #
    # Storage format depends of matrix type
    # This method raises at runtime if matrix doesn't have raw pointer
    def to_unsafe
      raise ArgumentError.new("#{self.class} can't be passed unsafe!")
    end

    # Returns shape of matrix in a form of tuple {nrows, ncolumns}
    def size : Tuple(Int32, Int32)
      {nrows, ncolumns}
    end

    # :ditto:
    def shape
      size
    end

    # Creates general matrix with same content. Useful for banded\sparse matrices
    def to_general
      GeneralMatrix(T).new(nrows, ncolumns, flags) do |i, j|
        unsafe_fetch(i, j)
      end
    end

    # Converts complex matrix to real part
    def to_real
      {% unless T == Complex %}
        {% raise "Only complex matrices have #to_real" %}
      {% end %}
      new_flags = flags.real
      map_f64(&.real).tap { |r| r.flags = new_flags }
    end

    # Converts complex matrix to imaginary part
    def to_imag
      {% unless T == Complex %}
        {% raise "Only complex matrices have #to_imag" %}
      {% end %}
      new_flags = flags.imag
      map_f64(&.imag).tap { |r| r.flags = new_flags }
    end

    # matrix product to given m
    # def *(m : Matrix(T))
    #   if ncolumns != m.nrows
    #     raise ArgumentError.new("matrix size should match ([#{nrows}x#{ncolumns}] * [#{m.nrows}x#{m.ncolumns}]")
    #   end
    #   result = GeneralMatrix(T).new(nrows, m.ncolumns) do |i, j|
    #     (0...ncolumns).sum { |k| self.unsafe_fetch(i, k)*m.unsafe_fetch(k, j) }
    #   end
    #   result.tap { |r| r.flags = self.flags.mult(m.flags) }
    # end

    # Multiplies matrix to scalar
    def *(k : Number)
      new_flags = self.flags.scale(false)
      map { |v| v*k }.tap { |result| result.flags = new_flags }
    end

    # :ditto:
    def *(k : Complex)
      new_flags = self.flags.scale(k.imag != 0)
      map_complex(&.*(k)).tap { |r| r.flags = new_flags }
    end

    # Multiplies matrix to -1
    def -
      result = self*(-1)
    end

    # Adds scalar value to every element
    def +(k : Number)
      return clone if k.zero?
      map { |v| v + k }.tap { |result| result.flags = MatrixFlags::None }
    end

    # :ditto:
    def +(k : Complex)
      return clone if k.zero?
      map_complex { |v| v + k }.tap { |result| result.flags = MatrixFlags::None }
    end

    # Substracts scalar value from every element
    def -(k : (Number | Complex))
      self + (-k)
    end

    # Divides matrix to scalar
    def /(k : Number | Complex)
      self*(T.new(1.0) / k)
    end

    protected def assert_same_size(m : Matrix)
      if ncolumns != m.ncolumns || nrows != m.nrows
        raise ArgumentError.new("matrix size doesn't match ([#{nrows}x#{ncolumns}] and [#{m.nrows}x#{m.ncolumns}]")
      end
    end

    # Returns element-wise sum
    #
    # This method raises if another matrix doesn't have same size
    def +(m : Matrix(T))
      assert_same_size(m)
      result = GeneralMatrix(T).new(nrows, ncolumns, flags.sum(m.flags)) do |i, j|
        self.unsafe_fetch(i, j) + m.unsafe_fetch(i, j)
      end
    end

    # Returns element-wise substract
    #
    # This method raises if another matrix doesn't have same size
    def -(m : Matrix(T))
      assert_same_size(m)
      result = GeneralMatrix(T).new(nrows, ncolumns, flags.sum(m.flags)) do |i, j|
        self.unsafe_fetch(i, j) - m.unsafe_fetch(i, j)
      end
    end

    # Returns transposed matrix
    def transpose
      return clone if flags.symmetric?
      GeneralMatrix(T).new(ncolumns, nrows, flags.transpose) do |i, j|
        unsafe_fetch(j, i)
      end
    end

    # Returns conjurgate transposed matrix
    #
    # result is same as `#transpose` for real matrices
    def conjtranspose
      {% if T != Complex %}
        return transpose
      {% end %}
      return clone if flags.hermitian?
      GeneralMatrix(T).new(ncolumns, nrows, flags.transpose) do |i, j|
        unsafe_fetch(j, i).conj
      end
    end

    # Returns kroneker product with matrix b
    def kron(b : Matrix(T))
      Matrix(T).kron(self, b)
    end

    # Same as tril in scipy - returns lower triangular or trapezoidal part of matrix
    #
    # Returns a matrix with all elements above k-th diagonal zeroed
    def tril(k = 0)
      x = GeneralMatrix(T).new(nrows, ncolumns, flags.tril(k <= 0, square?)) do |i, j|
        i >= j - k ? unsafe_fetch(i, j) : 0
      end
      x
    end

    # Same as triu in scipy - returns upper triangular or trapezoidal part of matrix
    #
    # Returns a matrix with all elements below k-th diagonal zeroed
    def triu(k = 0)
      x = GeneralMatrix(T).new(nrows, ncolumns, flags.triu(k >= 0, square?)) do |i, j|
        i <= j - k ? unsafe_fetch(i, j) : 0
      end
      x
    end

    # Converts matrix to string, with linefeeds before and after matrix:
    #
    # Output looks like:
    # ```
    # [1, 2, 3, .... 10]
    # [11, 12, 13, .... 20]
    # ...
    # [91, 92, 93, .... 100]
    # ```
    def to_s(io)
      to_custom(io, "\n[", ", ", "]\n[", "]\n\n")
    end

    # Converts matrix to string for inspection
    #
    # Output looks like:
    # ```
    # GeneralMatrix(Float64) (10x10, MatrixFlags::None)
    # [1, 2, 3, .... 10]
    # [11, 12, 13, .... 20]
    # ...
    # [91, 92, 93, .... 100]
    # ```
    def inspect(io)
      io << self.class << " (" << nrows << "x" << ncolumns << ", " << flags << "):"
      to_s(io)
    end

    # Yields every index
    #
    # `all` argument controls whether to yield all or non-empty elements for banded\sparse matrices
    # Example:
    # `m.each_index { |i, j| m[i, j] = -m[i, j] }`
    def each_index(*, all = false, &block)
      nrows.times do |row|
        ncolumns.times do |column|
          yield row, column
        end
      end
    end

    # Yields every element of matrix
    #
    # `all` argument controls whether to yield all or non-empty elements for banded\sparse matrices
    # Example:
    # `m.each { |v| raise "negative element found" if v < 0 }`
    def each(*, all = false, &block)
      each_index(all: all) { |i, j| yield(unsafe_fetch(i, j)) }
    end

    # Yields every element of matrix with corresponding row and column
    #
    # `all` argument controls whether to yield all or non-empty elements for banded\sparse matrices
    # Example:
    # m.each_with_index { |v, i, j| m2[i, j] = v }
    def each_with_index(*, all = false, &block)
      each_index(all: all) { |i, j| yield(unsafe_fetch(i, j), i, j) }
    end

    # For every element of matrix that is above main diagonal it
    # yields a block with value and with corresponding row and column
    #
    # This method is useful for symmetric matrices and similar cases
    # `include_diagonal` argument controls whether to include elements on main diagonal
    # `all` argument controls whether to yield all or non-empty elements for banded\sparse matrices
    def each_upper(*, diagonal = true, all = false, &block)
      nrows.times do |row|
        range = diagonal ? (row...ncolumns) : (row + 1...ncolumns)
        range.each do |column|
          yield unsafe_fetch(row, column), row, column
        end
      end
    end

    # For every element of matrix that is below main diagonal it
    # yields a block with value and with corresponding row and column
    #
    # This method is useful for symmetric matrices and similar cases
    # `include_diagonal` argument controls whether to include elements on main diagonal
    # `all` argument controls whether to yield all or non-empty elements for banded\sparse matrices
    def each_lower(*, diagonal = true, all = false, &block)
      nrows.times do |row|
        range = diagonal ? (0..row) : (0...row)
        range.each do |column|
          yield unsafe_fetch(row, column), row, column
        end
      end
    end

    # Compare with another matrix
    #
    # Returns True only if all element are exactly equal.
    # Use `#almost_eq` If you need approximate equality
    def ==(other)
      return false unless nrows == other.nrows && ncolumns == other.ncolumns
      each_with_index(all: true) do |value, row, column|
        return false if other.unsafe_fetch(row, column) != value
      end
      true
    end

    # Returns True if matrix is square and False otherwise
    def square?
      nrows == ncolumns
    end

    # Return matrix repeated `arows` times by vertical and `acolumns` times by horizontal
    def repmat(arows, acolumns)
      GeneralMatrix(T).new(nrows*arows, ncolumns*acolumns) do |i, j|
        unsafe_fetch(i % nrows, j % ncolumns)
      end
    end

    # Return element of matrix on row i and column j
    #
    # If i or j negative they are counted from end of matrix
    # If i>=nrows or j>=ncolumns exception is raised
    # Use `#unsafe_fetch` if you need to skip these checks
    def [](i : Int32, j : Int32)
      i += nrows if i < 0
      j += ncolumns if j < 0
      if j >= 0 && j < ncolumns && i >= 0 && i < nrows
        unsafe_fetch(i, j)
      else
        raise IndexError.new("access to [#{i}, #{j}] in matrix with size #{nrows}x#{ncolumns}")
      end
    end

    # Assign element of matrix on row i and column j
    #
    # If i or j negative they are counted from end of matrix
    # If i>=nrows or j>=ncolumns exception is raised
    # Use `#unsafe_set` if you need to skip these checks
    # Note this method reset all matrix flags
    def []=(i : Int32, j : Int32, value)
      i += nrows if i < 0
      j += ncolumns if j < 0
      if j >= 0 && j < ncolumns && i >= 0 && i < nrows
        unsafe_set(i, j, value)
      else
        raise IndexError.new("access to [#{i}, #{j}] in matrix with size #{nrows}x#{ncolumns}")
      end
    end

    # Return submatrix over given ranges.
    #
    # See `SubMatrix(T)` for details on submatrices
    # Example:
    # ```
    # m = GMat32.new([[1, 2, 3, 4], [5, 6, 7, 8], [9, 10, 11, 12], [13, 14, 15, 16]])
    # m1 = m[1..3, 2..3] # => [[7.0, 8.0], [11.0, 12.0], [15.0, 16.0]]
    # m2 = m[..3, 2]     # => [[3.0], [7.0], [11.0], [15.0]]
    # ```
    def [](arows : Range(Int32 | Nil, Int32 | Nil), acolumns : Range(Int32 | Nil, Int32 | Nil))
      arows_begin = arows.begin || 0
      acolumns_begin = acolumns.begin || 0

      arows_end = arows.end || nrows
      if arows.end # compensate if user wanted the endpoint only if not nil
        arows_end += arows.excludes_end? ? 0 : 1
      end

      acolumns_end = acolumns.end || ncolumns
      if acolumns.end
        acolumns_end += acolumns.excludes_end? ? 0 : 1
      end

      # ranges are converted to have int endpoints and are always exclusive
      arows = arows_begin...arows_end
      acolumns = acolumns_begin...acolumns_end

      start_row = arows.begin + (arows.begin < 0 ? nrows : 0)
      start_col = acolumns.begin + (acolumns.begin < 0 ? ncolumns : 0)
      total_rows = arows.end - start_row + (arows.end < 0 ? nrows : 0)
      total_cols = acolumns.end - start_col + (acolumns.end < 0 ? ncolumns : 0)
      SubMatrix(T).new(self, {start_row, start_col}, {total_rows, total_cols})
    end

    # :ditto:
    def [](row : Int32, acolumns : Range(Int32 | Nil, Int32 | Nil))
      self[row..row, acolumns]
    end

    # :ditto:
    def [](arows : Range(Int32 | Nil, Int32 | Nil), column : Int32)
      self[arows, column..column]
    end

    # Assign value to a given submatrix
    #
    # value can be a scalar or a matrix of same size as affected submatrix
    def []=(arows : Range(Int32, Int32), acolumns : Range(Int32, Int32), value)
      submatrix = self[arows, acolumns]
      if value.is_a? Matrix
        raise IndexError.new("submatrix size must match assigned value") unless submatrix.size == value.size
        submatrix.each_index { |i, j| submatrix.unsafe_set i, j, value.unsafe_fetch(i, j) }
      else
        submatrix.each_index { |i, j| submatrix.unsafe_set i, j, value }
      end
    end

    # :ditto:
    def []=(row : Int32, acolumns : Range(Int32, Int32), value)
      self[row..row, acolumns] = value
    end

    # :ditto:
    def []=(nrows : Range(Int32, Int32), column : Int32, value)
      self[nrows, column..column] = value
    end

    # Returns estimated tolerance of equality\inequality
    #
    # This value is used by default in `#almost_eq` compare
    def tolerance
      amax = 0.0
      # TODO - better iteration?
      each do |value|
        vabs = {% if T == Complex %}value.real.abs + value.imag.abs{% else %}value.abs{% end %}
        amax = vabs if amax < vabs
      end
      amax = 1.0 if amax < 1.0
      amax * nrows * ncolumns * 10*real_type_const(EPSILON)
    end

    # Approximately compare with `other` matrix
    #
    # Returns true if all elements are within `eps` from corresponding elements of `other` matrix
    def almost_eq(other : Matrix(T), eps)
      each_with_index(all: true) do |value, row, column|
        return false if (value - other.unsafe_fetch(row, column)).abs > eps
      end
      true
    end

    # :ditto:
    # Uses `#tolerance` as an `eps` by default
    def almost_eq(other : Matrix(T))
      almost_eq other, {self.tolerance, other.tolerance}.min
    end

    # Generate matrix of given size with elements randomly distributed from range 0.0..1.0
    def self.rand(nrows, ncolumns, rng = Random::DEFAULT)
      GeneralMatrix(T).new(nrows, ncolumns) { |i, j| rng.rand }
    end

    # Generate matrix of given size with all elements equal to zero
    def self.zeros(nrows, ncolumns)
      GeneralMatrix(T).new(nrows, ncolumns, MatrixFlags.for_diag(nrows == ncolumns))
    end

    # Generate matrix of given size with all elements equal to one
    def self.ones(nrows, ncolumns)
      GeneralMatrix(T).new(nrows, ncolumns) { |i, j| 1 }.tap do |m|
        if m.square?
          m.flags = MatrixFlags::Symmetric | MatrixFlags::Hermitian
        end
      end
    end

    # Alias for `#repmat`
    def self.repmat(a : Matrix(T), nrows, ncolumns)
      a.repmat(nrows, ncolumns)
    end

    # Returns diagonal matrix of given size with all diagonal elements equal to `value`
    def self.diag(nrows, ncolumns, value : Number | Complex)
      diag(nrows, ncolumns) { value }
    end

    # Returns diagonal matrix of given size with diagonal elements taken from array `values`
    def self.diag(nrows, ncolumns, values)
      # TODO - replace to BandedMatrix
      # This code breaks `expm` function
      # BandedMatrix(T).new(nrows, ncolumns, 0, 0, MatrixFlags.for_diag(nrows == ncolumns)) do |i, j|
      #   T.new(values[i])
      # end
      #
      GeneralMatrix(T).diag(nrows, ncolumns, values)
    end

    # Returns square diagonal matrix with diagonal elements taken from array `values`
    def self.diag(values)
      diag(values.size, values.size, values)
    end

    # Returns single column matrix with elements taken from array `values`
    def self.column(values)
      GeneralMatrix(T).new(values.size, 1, values)
    end

    # Returns single row matrix with elements taken from array `values`
    def self.row(values)
      GeneralMatrix(T).new(1, values.size, values)
    end

    # Returns diagonal matrix of given size with diagonal elements equal to block value
    def self.diag(nrows, ncolumns, &block)
      BandedMatrix(T).new(nrows, ncolumns, 0, 0, MatrixFlags.for_diag(nrows == ncolumns)) do |i, j|
        yield i
      end
    end

    # Returns kroneker product of matrices
    #
    # Resulting matrix size is `{a.nrows*b.nrows, a.ncolumns*b.ncolumns}`
    def self.kron(a, b)
      GeneralMatrix(T).new(a.nrows*b.nrows, a.ncolumns*b.ncolumns) do |i, j|
        a.unsafe_fetch(i // b.nrows, j // b.ncolumns) * b.unsafe_fetch(i % b.nrows, j % b.ncolumns)
      end
    end

    # returns identity matrix of size `n`
    def self.identity(n)
      (diag(n, n) { |i| 1 }).tap do |r|
        r.assume! MatrixFlags::PositiveDefinite
      end
    end

    # :ditto:
    def self.eye(n)
      self.identity(n)
    end

    # Create row from start_val...end_val with step of delta between
    def self.arange(start_val : T, end_val : T, delta = 1.0)
      return GeneralMatrix(T).new(1, 0) unless (end_val - start_val).sign == delta.sign && delta.abs > 0
      GeneralMatrix(T).new(1, ((end_val - start_val).abs.ceil / delta.abs).ceil.to_i) do |i, j|
        start_val + j*delta
      end
    end

    # alias to `#transpose`
    def t
      transpose
    end

    # alias to `#transpose!`
    def t!
      transpose!
    end

    # alias to `#conjtranspose`
    def conjt
      conjtranspose
    end

    # alias to `#conjtranspose!`
    def conjt!
      conjtranspose!
    end

    # Returns concatenation with another matrix by `Axis::Rows` (horizontal) or `Axis::Columns` (vertical)
    def cat(other : Matrix(T), axis : Axis)
      if self.size[1 - axis.to_i] != other.size[1 - axis.to_i]
        raise ArgumentError.new("matrix size along other axis should match for concatenation")
      end
      case axis
      in Axis::Columns
        GeneralMatrix(T).new(nrows + other.nrows, ncolumns) do |row, column|
          row < nrows ? unsafe_fetch(row, column) : other.unsafe_fetch(row - nrows, column)
        end
      in Axis::Rows
        GeneralMatrix(T).new(nrows, ncolumns + other.ncolumns) do |row, column|
          column < ncolumns ? unsafe_fetch(row, column) : other.unsafe_fetch(row, column - ncolumns)
        end
      end
    end

    # Returns vertical concatenation with another matrix
    def vcat(other)
      cat other, Axis::Columns
    end

    # Returns horizontal concatenation with another matrix
    def hcat(other)
      cat other, axis: Axis::Rows
    end

    # Yields each element with index and replace it with returned value
    def map_with_index!(&block)
      each_with_index { |v, i, j| unsafe_set(i, j, yield(v, i, j)) }
      self
    end

    # Yields each element (without index) and replace it with returned value
    def map!(&block)
      map_with_index! { |v, i, j| yield(v) }
      self
    end

    # Returns result of appliyng block to every element with index
    def map_with_index(&block)
      GeneralMatrix(T).new(nrows, ncolumns) { |i, j| yield(unsafe_fetch(i, j), i, j) }
    end

    # Returns result of appliyng block to every element (without index)
    def map(&block)
      map_with_index { |v, i, j| yield(v) }
    end

    # TODO - macro magic?
    protected def map_with_index_f64(&block)
      GeneralMatrix(Float64).new(nrows, ncolumns) { |i, j| yield(unsafe_fetch(i, j), i, j) }
    end

    protected def map_with_index_complex(&block)
      GeneralMatrix(Complex).new(nrows, ncolumns) { |i, j| yield(unsafe_fetch(i, j), i, j) }
    end

    protected def map_f64(&block)
      map_with_index_f64 { |v, i, j| yield(v) }
    end

    protected def map_complex(&block)
      map_with_index_complex { |v, i, j| yield(v) }
    end

    # Works like a tril in scipy - remove all elements above k-diagonal
    def tril!(k = 0)
      oldflags = flags
      map_with_index! { |v, i, j| i < j - k ? 0 : v }
      self.flags = oldflags.tril(k <= 0, square?)
      self
    end

    # Works like a triu in scipy - remove all elements below k-diagonal
    def triu!(k = 0)
      oldflags = flags
      map_with_index! { |v, i, j| i > j - k ? 0 : v }
      self.flags = oldflags.triu(k >= 0, square?)
      self
    end

    # Returns sum of diagonal elements (trace of matrix)
    def trace
      diag.sum
    end

    # Perform inplace addition with matrix `m` multiplied to scalar `k`
    #
    # `a.add!(k, b)` is equal to `a = a + k * b`, but faster as no new matrix is allocated
    def add!(k : Number, m : Matrix)
      assert_same_size(m)
      oldflags = flags
      map_with_index! { |v, i, j| v + k*m.unsafe_fetch(i, j) }
      self.flags = oldflags.sum(m.flags.scale(k.is_a?(Complex) && k.imag != 0))
      self
    end

    # Performs inplace addition with matrix `m`
    #
    # `a.add!(b)` is equal to `a = a + b`, but faster as no new matrix is allocated
    def add!(m)
      add!(1, m)
    end

    # Perform inplace multiplication to scalar `k`
    def scale!(k : Number | Complex)
      oldflags = flags
      map_with_index! { |v, i, j| k*v }
      new_flags = oldflags.scale(k.is_a?(Complex) && k.imag != 0)
      self
    end

    # Perform `reduce` from `initial` value over given `axis`
    def reduce(axis : Axis, initial, &block)
      case axis
      in Axis::Columns
        GeneralMatrix(T).new(1, ncolumns) do |_, column|
          result = T.new(initial)
          nrows.times { |row| result = yield(result, unsafe_fetch(row, column)) }
          result
        end
      in Axis::Rows
        GeneralMatrix(T).new(nrows, 1) do |row, _|
          result = T.new(initial)
          ncolumns.times { |column| result = yield(result, unsafe_fetch(row, column)) }
          result
        end
      end
    end

    # Perform sum over given `axis`
    def sum(axis : Axis)
      reduce(axis, 0) { |memo, e| memo + e }
    end

    # Perform product over given `axis`
    def product(axis : Axis)
      reduce(axis, 1) { |memo, e| memo * e }
    end

    # Calculate maximum over given `axis`
    def max(axis : Axis)
      reduce(axis, -T::INFINITY) { |memo, e| {memo, e}.max }
    end

    # Calculate minimum over given `axis`
    def min(axis : Axis)
      reduce(axis, T::INFINITY) { |memo, e| {memo, e}.min }
    end

    # Converts complex matrix to real one if all imaginary parts are less then `eps`, returns `nil` otherwise
    #
    # Returns just matrix if it is already real
    def chop(eps = self.tolerance)
      {% if T == Complex %}
        if all? { |v| v.imag.abs < eps }
          self.to_real
        else
          nil
        end
      {% else %}
        self
      {% end %}
    end
  end

  module Aliases
    alias Mat = Matrix(Float64)
    alias Mat32 = Matrix(Float32)
    alias MatComplex = Matrix(Complex)
  end
end
