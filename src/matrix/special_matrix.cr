require "complex"
require "./matrix"
require "./general_matrix"

module LA
  abstract class Matrix(T)
    # Construct (nrows, ncolumns) matrix filled with ones at and below the kth diagonal.
    #
    # The matrix has A[i,j] == 1 for j <= i + k
    #
    # `k` - Number of subdiagonal below which matrix is filled with ones. k = 0 is the main diagonal, k < 0 subdiagonal and k > 0 superdiagonal.
    #
    # Behaviour copied from [scipy](https://docs.scipy.org/doc/scipy/reference/generated/scipy.linalg.tri.html)
    #
    # Example:
    # ```
    # Mat.tri(3, 5, 2).to_aa # => [[
    # # [1, 1, 1, 0, 0],
    # # [1, 1, 1, 1, 0],
    # # [1, 1, 1, 1, 1]]
    # Mat.tri(3, 5, -1).to_aa # => [[
    # # [0, 0, 0, 0, 0],
    # # [1, 0, 0, 0, 0],
    # # [1, 1, 0, 0, 0]]
    # ```
    def self.tri(nrows, ncolumns, k = 0)
      flags = k <= 0 ? MatrixFlags::LowerTriangular : MatrixFlags::None
      GeneralMatrix(T).new(nrows, ncolumns, flags) do |i, j|
        i >= j - k ? 1 : 0
      end
    end

    # Create a block diagonal matrix from provided matrices
    #
    # Given the inputs A, B and C, the output will have these matrices arranged on the diagonal:
    #
    # Example:
    # ```
    # m = Mat.block_diag(a, b, c)
    # ```
    # m will have following structure:
    # ```
    # [[a, 0, 0],
    #  [0, b, 0],
    #  [0, 0, c]]
    # ```
    def self.block_diag(*args)
      nrows = args.sum &.nrows
      ncolumns = args.sum &.ncolumns
      GeneralMatrix(T).new(nrows, ncolumns).tap do |result|
        row = 0
        column = 0
        args.each do |arg|
          result[row...row + arg.nrows, column...column + arg.ncolumns] = arg
          row += arg.nrows
          column += arg.ncolumns
        end
      end
    end

    # Create a [Toeplitz matrix](https://en.wikipedia.org/wiki/Toeplitz_matrix)
    #
    # `column` - first column of matrix
    #
    # `row` - first row of matrix (if nil, it is assumed that `row = column.map(&.conj)`)
    #
    # Behaviour copied from [scipy](https://docs.scipy.org/doc/scipy/reference/generated/scipy.linalg.toeplitz.html)
    #
    # Example:
    # ```
    # MatComplex.toeplitz([1, 2, 3], [1, 4, 5, 6]) # =>
    # # LA::GeneralMatrix(Float64) (3x4, None):
    # # [1.0, 4.0, 5.0, 6.0]
    # # [2.0, 1.0, 4.0, 5.0]
    # # [3.0, 2.0, 1.0, 4.0]
    # MatComplex.toeplitz([1.0, 2 + 3.i, 4 - 1.i]) # =>
    # # LA::GeneralMatrix(Complex) (3x3, Hermitian):
    # # [1.0 + 0.0i, 2.0 - 3.0i, 4.0 + 1.0i]
    # # [2.0 + 3.0i, 1.0 + 0.0i, 2.0 - 3.0i]
    # # [4.0 - 1.0i, 2.0 + 3.0i, 1.0 + 0.0i]
    # ```
    def self.toeplitz(column : Indexable | Matrix, row : Indexable | Matrix | Nil = nil)
      row = row.to_a if row.is_a? Matrix
      column = column.to_a if column.is_a? Matrix
      if row
        GeneralMatrix(T).new(column.size, row.size) do |i, j|
          k = i - j
          if k >= 0
            column[k]
          else
            row[-k]
          end
        end
      else
        GeneralMatrix(T).new(column.size, column.size, MatrixFlags::Hermitian) do |i, j|
          k = i - j
          if k >= 0
            column[k]
          else
            {% if T == Complex %}
              column[-k].conj
            {% else %}
              column[-k]
            {% end %}
          end
        end
      end
    end

    # Construct a [Circulant matrix](https://en.wikipedia.org/wiki/Circulant_matrix)
    #
    # c - first column of matrix
    #
    # Example:
    # ```
    # a = circulant([1, 2, 3])
    # a.to_aa # => [[1, 3, 2],[2, 1, 3],[3, 2, 1]])
    # ```
    # Behaviour copied from [scipy](https://docs.scipy.org/doc/scipy/reference/generated/scipy.linalg.circulant.html)
    def self.circulant(c)
      GeneralMatrix(T).new(c.size, c.size) do |i, j|
        k = i - j
        c[(k + c.size) % c.size]
      end
    end

    # Create a Leslie matrix
    #
    # Given the length n array of fecundity coefficients `f` and the length n-1 array of survival coefficients `s`, return the associated Leslie matrix.
    #
    # Behaviour copied from [scipy](https://docs.scipy.org/doc/scipy/reference/generated/scipy.linalg.leslie.html)
    #
    # Example:
    # ```
    # Mat.leslie([0.1, 2.0, 1.0, 0.1], [0.2, 0.8, 0.7]) # =>
    # # LA::GeneralMatrix(Float64) (4x4, None):
    # # [0.1, 2.0, 1.0, 0.1]
    # # [0.2, 0.0, 0.0, 0.0]
    # # [0.0, 0.8, 0.0, 0.0]
    # # [0.0, 0.0, 0.7, 0.0]
    # ```
    def self.leslie(f, s)
      GeneralMatrix(T).new(s.size + 1, f.size).tap do |matrix|
        f.each_with_index { |fi, i| matrix.unsafe_set 0, i, to_my_type(fi) }
        s.each_with_index { |si, i| matrix.unsafe_set i + 1, i, to_my_type(si) }
      end
    end

    # Create a companion matrix associated with the polynomial whose coefficients are given in `a`
    #
    # Behaviour copied from [scipy](https://docs.scipy.org/doc/scipy/reference/generated/scipy.linalg.companion.html)
    #
    # Example:
    # ```
    # Mat.companion([1, -10, 31, -30]) # =>
    # # LA::GeneralMatrix(Float64) (3x3, None):
    # # [10.0, -31.0, 30.0]
    # # [1.0, 0.0, 0.0]
    # # [0.0, 1.0, 0.0]
    # ```
    def self.companion(a)
      k = -1.0/a[0]
      GeneralMatrix(T).new(a.size - 1, a.size - 1).tap do |matrix|
        (a.size - 1).times { |i| matrix.unsafe_set 0, i, to_my_type(a[i + 1]*k) }
        (a.size - 2).times { |i| matrix.unsafe_set i + 1, i, to_my_type(1) }
      end
    end

    # Constructs an n-by-n Hadamard matrix, n must be power of 2
    #
    # Behaviour copied from [scipy](https://docs.scipy.org/doc/scipy/reference/generated/scipy.linalg.hadamard.html)
    #
    # Example:
    # ```
    # Mat.hadamard(4) # =>
    # # LA::GeneralMatrix(Float64) (4x4,  Symmetric | Hermitian):
    # # [1.0, 1.0, 1.0, 1.0]
    # # [1.0, -1.0, 1.0, -1.0]
    # # [1.0, 1.0, -1.0, -1.0]
    # # [1.0, -1.0, -1.0, 1.0]
    # ```
    def self.hadamard(n)
      # TODO - faster implementation
      raise ArgumentError.new("size must be positive") unless n > 0
      raise ArgumentError.new("size must be power of two") unless n.popcount == 1
      m = case n
          when 1
            GeneralMatrix(T).new([[1]])
          when 2
            GeneralMatrix(T).new([[1, 1], [1, -1]])
          else
            hadamard(n//2).kron(hadamard(2))
          end
      m.assume!(MatrixFlags::Symmetric | MatrixFlags::Hermitian)
      m
    end

    # Create a Hankel matrix
    # The Hankel matrix has constant anti-diagonals, with `column` as its first column and `row` as its last row.
    # If `row is nil, then row with elements equal to zero is assumed.
    # Behaviour copied from [scipy](https://docs.scipy.org/doc/scipy/reference/generated/scipy.linalg.hankel.html)
    # Examples:
    # ```
    # a = Mat.hankel([1, 17, 99])
    # a.to_aa # => [
    # # [ 1, 17, 99],
    # # [17, 99,  0],
    # # [99,  0,  0]]
    # a = Mat.hankel([1, 2, 3, 4], [4, 7, 7, 8, 9])
    # a.to_aa # => [
    # # [1, 2, 3, 4, 7],
    # # [2, 3, 4, 7, 7],
    # # [3, 4, 7, 7, 8],
    # # [4, 7, 7, 8, 9]]
    # ```
    def self.hankel(column : Indexable | Matrix, row : Indexable | Matrix | Nil = nil)
      row = row.to_a if row.is_a? Matrix
      column = column.to_a if column.is_a? Matrix
      if row
        GeneralMatrix(T).new(column.size, row.size) do |i, j|
          k = i + j
          if k < column.size
            column[k]
          else
            row[k - column.size + 1]
          end
        end
      else
        GeneralMatrix(T).new(column.size, column.size, MatrixFlags::Symmetric) do |i, j|
          k = i + j
          if k < column.size
            column[k]
          else
            0
          end
        end
      end
    end

    # Create an Helmert matrix of order n
    #
    # If `full` is true the (n x n) matrix will be returned.
    # Otherwise the submatrix that does not include the first row will be returned
    #
    # Behaviour copied from [scipy](https://docs.scipy.org/doc/scipy/reference/generated/scipy.linalg.helmert.html)
    #
    # Example:
    # ```
    # Mat.helmert(5, full: true) # =>
    # # LA::GeneralMatrix(Float64) (5x5, Orthogonal):
    # # [0.4472135954999579, 0.4472135954999579, 0.4472135954999579, 0.4472135954999579, 0.4472135954999579]
    # # [0.7071067811865476, -0.7071067811865476, 0.0, 0.0, 0.0]
    # # [0.408248290463863, 0.408248290463863, -0.816496580927726, 0.0, 0.0]
    # # [0.28867513459481287, 0.28867513459481287, 0.28867513459481287, -0.8660254037844386, 0.0]
    # # [0.22360679774997896, 0.22360679774997896, 0.22360679774997896, 0.22360679774997896, -0.8944271909999159]
    # ```
    def self.helmert(n, full = false)
      if full
        result = GeneralMatrix(T).new(n, n)
      else
        result = GeneralMatrix(T).new(n - 1, n)
      end
      # first row
      if full
        result[0, 0...n] = to_my_type(Math.sqrt(1.0/n))
        rowdelta = 1
      else
        rowdelta = 0
      end
      # rest
      (n - 1).times do |i|
        x = i + 1
        v = to_my_type(Math.sqrt(1.0/(x + x*x)))
        result.unsafe_set i + rowdelta, i + 1, -v*x
        result[i + rowdelta, 0..i] = v
      end
      result.assume!(MatrixFlags::Orthogonal) if full
      result
    end

    # Create a Hilbert matrix of order n.
    #
    # Returns the n by n matrix with entries h[i,j] = 1 / (i + j + 1).
    #
    # Behaviour copied from [scipy](https://docs.scipy.org/doc/scipy/reference/generated/scipy.linalg.hilbert.html)
    #
    # Example:
    # ```
    # Mat.hilbert(3) # =>
    # # LA::GeneralMatrix(Float64) (3x3, Symmetric | Hermitian | PositiveDefinite):
    # # [1.0, 0.5, 0.3333333333333333]
    # # [0.5, 0.3333333333333333, 0.25]
    # # [0.3333333333333333, 0.25, 0.2]
    # ```
    def self.hilbert(n)
      GeneralMatrix(T).new(n, n, MatrixFlags::Symmetric | MatrixFlags::Hermitian | MatrixFlags::PositiveDefinite) do |i, j|
        T.multiplicative_identity / (i + j + 1)
      end
    end
  end

  enum Enums::DFTScale
    None
    N
    SqrtN
  end

  abstract class Matrix(T)
    def self.dft(n, scale : DFTScale = DFTScale::None)
      {% raise "DFT matrix must be Complex" unless T == Complex %}
      j = Complex.new(0, 1)
      w = Math.exp(-2*Math::PI*j / n)
      result = Matrix(T).ones(n, n).clone
      result.each_index do |i, j|
        next if i == 0 || j == 0
        if j == 1
          result.unsafe_set(i, j, w*result.unsafe_fetch(i - 1, j))
        else
          result.unsafe_set(i, j, result.unsafe_fetch(i, 1)*result.unsafe_fetch(i, j - 1))
        end
      end
      case scale
      in .sqrt_n?
        scale = 1.0 / Math.sqrt(n)
        result.map! { |v| scale*v }
      in .n?
        scale = 1.0 / n
        result.map! { |v| scale*v }
      in .none?
        # do nothing
      end
      result
    end
  end

  enum Enums::PascalKind
    Upper
    Lower
    Symmetric
  end

  abstract class Matrix(T)
    private def self.n_choose_k(n, k)
      (1..n).product(T.multiplicative_identity) / ((1..k).product(T.multiplicative_identity) * (1..n - k).product(T.multiplicative_identity))
    end

    # Returns the n x n Pascal matrix.
    #
    # The Pascal matrix is a matrix containing the binomial coefficients as its elements.
    #
    # `kind` : see `PascalKind`
    #
    # Behaviour copied from [scipy](https://docs.scipy.org/doc/scipy/reference/generated/scipy.linalg.pascal.html)
    #
    # Example:
    # ```
    # Mat.pascal(4) # =>
    # # LA::GeneralMatrix(Float64) (4x4, Symmetric | Hermitian | PositiveDefinite):
    # # [1.0, 1.0, 1.0, 1.0]
    # # [1.0, 2.0, 3.0, 4.0]
    # # [1.0, 3.0, 6.0, 10.0]
    # # [1.0, 4.0, 10.0, 20.0]
    # Mat.pascal(4, PascalKind::Lower) # =>
    # # LA::GeneralMatrix(Float64) (4x4, LowerTriangular):
    # # [1.0, 0.0, 0.0, 0.0]
    # # [1.0, 1.0, 0.0, 0.0]
    # # [1.0, 2.0, 1.0, 0.0]
    # # [1.0, 3.0, 3.0, 1.0]
    # ```
    def self.pascal(n, kind : PascalKind = PascalKind::Symmetric)
      case kind
      in .upper?
        zeros(n, n).tap do |m|
          m.each_index do |i, j|
            next if i > j
            if i == j || i == 0
              m.unsafe_set(i, j, to_my_type(1.0))
            else
              m.unsafe_set(i, j, m.unsafe_fetch(i - 1, j - 1) + m.unsafe_fetch(i, j - 1))
            end
          end
          m.assume!(MatrixFlags::UpperTriangular)
        end
      in .lower?
        zeros(n, n).tap do |m|
          m.each_index do |i, j|
            next if i < j
            if i == j || j == 0
              m.unsafe_set(i, j, to_my_type(1.0))
            else
              m.unsafe_set(i, j, m.unsafe_fetch(i - 1, j - 1) + m.unsafe_fetch(i - 1, j))
            end
          end
          m.assume!(MatrixFlags::LowerTriangular)
        end
      in .symmetric?
        zeros(n, n).tap do |m|
          m.each_index do |i, j|
            if i == 0 || j == 0
              m.unsafe_set(i, j, to_my_type(1.0))
            else
              m.unsafe_set(i, j, m.unsafe_fetch(i - 1, j) + m.unsafe_fetch(i, j - 1))
            end
          end
          m.assume!(MatrixFlags::Symmetric | MatrixFlags::Hermitian | MatrixFlags::PositiveDefinite)
        end
      end
    end

    # Returns the inverse of the n x n Pascal matrix
    #
    # The Pascal matrix is a matrix containing the binomial coefficients as its elements.
    #
    # `kind` : see `PascalKind`
    #
    # Behaviour copied from [scipy](https://docs.scipy.org/doc/scipy/reference/generated/scipy.linalg.invpascal.html)
    #
    # Example:
    # ```
    # Mat.invpascal(4) # =>
    # # LA::GeneralMatrix(Float64) (5x5, Symmetric | Hermitian | PositiveDefinite):
    # # [5.0, -10.0, 10.0, -5.0, 1.0]
    # # [-10.0, 30.0, -35.0, 19.0, -4.0]
    # # [10.0, -35.0, 46.0, -27.0, 6.0]
    # # [-5.0, 19.0, -27.0, 17.0, -4.0]
    # # [1.0, -4.0, 6.0, -4.0, 1.0]
    # Mat.invpascal(5, PascalKind::Lower) # =>
    # # LA::GeneralMatrix(Float64) (5x5, LowerTriangular):
    # # [1.0, 0.0, 0.0, 0.0, 0.0]
    # # [-1.0, 1.0, 0.0, 0.0, 0.0]
    # # [1.0, -2.0, 1.0, 0.0, 0.0]
    # # [-1.0, 3.0, -3.0, 1.0, 0.0]
    # # [1.0, -4.0, 6.0, -4.0, 1.0]
    # ```
    def self.invpascal(n, kind : PascalKind = PascalKind::Symmetric)
      case kind
      when .symmetric?
        # TODO - better method
        result = invpascal(n, PascalKind::Upper)*invpascal(n, PascalKind::Lower)
        result.assume!(MatrixFlags::Symmetric | MatrixFlags::Hermitian | MatrixFlags::PositiveDefinite)
        result
      else
        pascal(n, kind).map_with_index! { |v, r, c| (r - c) % 2 == 0 ? v : -v }.tap do |m|
          if kind.lower?
            m.assume!(MatrixFlags::LowerTriangular)
          else
            m.assume!(MatrixFlags::UpperTriangular)
          end
        end
      end
    end

    # Returns a symmetric Fiedler matrix
    #
    # Given an sequence of numbers `values`, Fiedler matrices have the structure `f[i, j] = (values[i] - values[j]).abs`, and hence zero diagonals and nonnegative entries.
    # A Fiedler matrix has a dominant positive eigenvalue and other eigenvalues are negative.
    #
    # Behaviour copied from [scipy](https://docs.scipy.org/doc/scipy/reference/generated/scipy.linalg.fiedler.html)
    #
    # Example:
    # ```
    # Mat.fiedler([1, 4, 12, 45, 77]) # =>
    # # LA::GeneralMatrix(Float64) (5x5, Symmetric | Hermitian):
    # # [0.0, 3.0, 11.0, 44.0, 76.0]
    # # [3.0, 0.0, 8.0, 41.0, 73.0]
    # # [11.0, 8.0, 0.0, 33.0, 65.0]
    # # [44.0, 41.0, 33.0, 0.0, 32.0]
    # # [76.0, 73.0, 65.0, 32.0, 0.0]
    # ```
    def self.fiedler(values)
      GeneralMatrix(T).new(values.size, values.size, flags: MatrixFlags::Symmetric | MatrixFlags::Hermitian) do |i, j|
        (values[i] - values[j]).abs
      end
    end

    # Compute the inverse of the Hilbert matrix of order n.
    #
    # The entries in the inverse of a Hilbert matrix are integers.
    #
    # Behaviour copied from [scipy](https://docs.scipy.org/doc/scipy/reference/generated/scipy.linalg.invhilbert.html)
    #
    # Example:
    # ```
    # Mat.invhilbert(4) # => LA::GeneralMatrix(Float64) (4x4, Symmetric | Hermitian | PositiveDefinite):
    # # [16.0, -120.0, 240.0, -140.0]
    # # [-120.0, 1200.0, -2700.0, 1680.0]
    # # [240.0, -2700.0, 6480.0, -4200.0]
    # # [-140.0, 1680.0, -4200.0, 2800.0]
    #
    # Mat.invhilbert(16)[7, 7] # => 4.247509952853739e+19
    # ```
    def self.invhilbert(n)
      GeneralMatrix(T).new(n, n, flags: MatrixFlags::Symmetric | MatrixFlags::Hermitian | MatrixFlags::PositiveDefinite) do |i1, j1|
        i, j = i1 + 1, j1 + 1
        ((i + j).odd? ? -1 : 1)*(i + j - 1)*n_choose_k(n + i - 1, n - j)*n_choose_k(n + j - 1, n - i)*(n_choose_k(i + j - 2, i - 1)**2)
      end
    end
  end
end
