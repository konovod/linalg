require "../matrix/*"
require "./lapack_helper"

module LA
  module Enums
    # Enumeration of methods for solving least squares problems.
    #
    # - `Auto`: Automatically choose the best method (currently QR).
    # - `QR`: Use QR factorization (fast, good for well-conditioned problems).
    # - `Orthogonal`: Use orthogonal factorization with pivoting (good for rank-deficient problems).
    # - `SVD`: Use singular value decomposition (most robust, handles rank deficiency well).
    # - `LS`: Alias for QR.
    # - `LSY`: Alias for Orthogonal.
    # - `LSD`: Alias for SVD.
    enum LSMethod
      Auto       = 0
      QR
      Orthogonal
      SVD
      LS         = QR
      LSY        = Orthogonal
      LSD        = SVD
    end

    # Enumeration of methods for rank estimation.
    #
    # - `SVD`: Use singular value decomposition (accurate but slower).
    # - `QRP`: Use QR factorization with pivoting (faster but less accurate for ill-conditioned matrices).
    enum RankMethod
      SVD
      QRP
    end

    # Enumeration of matrix norm types.
    #
    # - `Frobenius`: Frobenius norm (square root of sum of squares of all elements).
    # - `One`: Maximum column sum (1-norm).
    # - `Inf`: Maximum row sum (infinity-norm).
    # - `MaxAbs`: Maximum absolute value of any element.
    enum MatrixNorm
      Frobenius
      One
      # Two
      Inf
      MaxAbs
    end
  end

  # Exception raised for linear algebra errors.
  class Utils::LinAlgError < Exception
  end

  # Returns the LAPACK version as a tuple {major, minor, patch}.
  #
  # LAPACK Routine:
  #   - Uses `ilaver` (LAPACK version query).
  def self.lapack_version
    LibLAPACK.ilaver(out major, out minor, out patch)
    {major, minor, patch}
  end

  # Calculate matrix inversion.
  #
  # Arguments:
  #   - matrix (Matrix(T)) : The matrix to invert.
  #   - overwrite_a (Bool) : If `true`, allows overwriting the input matrix. Default: `false`.
  #
  # Returns:
  #   - GeneralMatrix(T) : The inverted matrix.
  #
  # See `GeneralMatrix#inv!` for details of algorithm.
  def self.inv(matrix, *, overwrite_a = false)
    overwrite_a ? matrix.inv! : matrix.inv
  end

  # Solves the linear system A * X = B.
  #
  # Arguments:
  #   - a (Matrix(T)) : Coefficient matrix.
  #   - b (Matrix(T)) : Right-hand side matrix.
  #   - overwrite_a (Bool) : If `true`, allows overwriting `a`. Default: `false`.
  #   - overwrite_b (Bool) : If `true`, allows overwriting `b`. Default: `false`.
  #
  # Returns:
  #   - GeneralMatrix(T) : Solution matrix X.
  #
  # See `Matrix#solve` for details.
  def self.solve(a, b, *, overwrite_a = false, overwrite_b = false)
    a.solve(b, overwrite_a: overwrite_a, overwrite_b: overwrite_b)
  end

  # Solves the linear least squares problem min ||A*X - B||.
  #
  # Arguments:
  #   - a (Matrix(T)) : Coefficient matrix.
  #   - b (Matrix(T)) : Right-hand side matrix.
  #   - method (LSMethod) : Algorithm to use. Default: `LSMethod::Auto`.
  #   - overwrite_a (Bool) : If `true`, allows overwriting `a`. Default: `false`.
  #   - overwrite_b (Bool) : If `true`, allows overwriting `b`. Default: `false`.
  #   - cond (Number) : Cutoff for rank determination. Default: -1 (machine precision).
  #
  # Returns:
  #   - Tuple(GeneralMatrix(T), Int32, Array(T)) : Solution X, effective rank, singular values (if applicable).
  #
  # See `Matrix#lstsq` for details.
  def self.lstsq(a, b, method : LSMethod = LSMethod::Auto, *, overwrite_a = false, overwrite_b = false, cond = -1)
    a.lstsq(b, method, overwrite_a: overwrite_a, overwrite_b: overwrite_b, cond: cond)
  end

  # Solves the linear least squares problem using QR factorization.
  #
  # Arguments:
  #   - a (Matrix(T)) : Coefficient matrix.
  #   - b (Matrix(T)) : Right-hand side matrix.
  #   - overwrite_a (Bool) : If `true`, allows overwriting `a`. Default: `false`.
  #   - overwrite_b (Bool) : If `true`, allows overwriting `b`. Default: `false`.
  #   - cond (Number) : Cutoff for rank determination. Default: -1 (machine precision).
  #
  # Returns:
  #   - GeneralMatrix(T) : Solution matrix X.
  #
  # See `Matrix#solvels` for details.
  def self.solvels(a, b, *, overwrite_a = false, overwrite_b = false, cond = -1)
    a.solvels(b, overwrite_a: overwrite_a, overwrite_b: overwrite_b, cond: cond)
  end

  # Computes the singular value decomposition (SVD) of a matrix.
  #
  # Arguments:
  #   - matrix (Matrix(T)) : Input matrix.
  #   - overwrite_a (Bool) : If `true`, allows overwriting `matrix`. Default: `false`.
  #
  # Returns:
  #   - Tuple(GeneralMatrix(T), Array(T), GeneralMatrix(T)) : U, singular values, V^T.
  #
  # See `Matrix#svd` for details.
  def self.svd(matrix, *, overwrite_a = false)
    matrix.svd(overwrite_a: overwrite_a)
  end

  abstract class Matrix(T)
    # Returns the uplo character ('L' or 'U') based on matrix flags.
    private def uplo
      flags.lower_triangular? ? 'L'.ord.to_u8 : 'U'.ord.to_u8
    end

    # Ensures symmetry by copying upper triangle to lower triangle.
    private def adjust_symmetric
      f = flags
      each_upper(diagonal: false) { |v, i, j| unsafe_set(j, i, v) }
      self.flags = f
    end

    # Ensures triangular form by zeroing out the opposite triangle.
    private def adjust_triangular
      triu! if flags.upper_triangular?
      tril! if flags.lower_triangular?
    end

    # Calculate matrix inversion.
    #
    # Returns:
    #   - GeneralMatrix(T) : The inverted matrix.
    #
    # See `GeneralMatrix#inv!` for details on algorithm.
    def inv
      to_general.inv!
    end

    # Computes the singular value decomposition (SVD) of a matrix.
    #
    # Returns:
    #   - Tuple(GeneralMatrix(T), Array(T), GeneralMatrix(T)) : U, singular values, V.t.
    #
    # See `GeneralMatrix#svd` for details.
    def svd
      to_general.svd(overwrite_a: true)
    end

    # Computes only the singular values of a matrix.
    #
    # Returns:
    #   - Array(T) : Singular values.
    #
    # See `GeneralMatrix#svdvals` for details.
    def svdvals
      to_general.svdvals(overwrite_a: true)
    end

    # Solves the linear system A * X = B.
    #
    # Arguments:
    #   - b (GeneralMatrix(T)) : Right-hand side matrix.
    #   - overwrite_b (Bool) : If `true`, allows overwriting `b`. Default: `false`.
    #
    # Returns:
    #   - GeneralMatrix(T) : Solution matrix X.
    #
    # See `GeneralMatrix#solve` for details.
    def solve(b : GeneralMatrix(T), *, overwrite_b = false)
      to_general.solve(b, overwrite_b: overwrite_b)
    end

    # Solves the linear system A * X = B.
    #
    # Arguments:
    #   - b (Matrix(T)) : Right-hand side matrix.
    #
    # Returns:
    #   - GeneralMatrix(T) : Solution matrix X.
    def solve(b : self)
      to_general.solve(b.to_general, overwrite_a: true, overwrite_b: true)
    end

    # Calculates the determinant of a square matrix.
    #
    # Returns:
    #   - T : Determinant value.
    #
    # See `GeneralMatrix#det` for details.
    def det
      to_general.det(overwrite_a: true)
    end

    # Balances a square matrix to improve eigenvalue accuracy.
    #
    # Arguments:
    #   - permute (Bool) : If `true`, permutes to isolate eigenvalues. Default: `true`.
    #   - scale (Bool) : If `true`, scales to improve conditioning. Default: `true`.
    #   - separate (Bool) : If `true`, returns scaling factors separately. Default: `false`.
    #
    # Returns:
    #   - Tuple(GeneralMatrix(T), GeneralMatrix(T)) : Balanced matrix and scaling factors (or diagonal).
    #
    # LAPACK Routine:
    #   - Uses `gebal` (balance matrix).
    def balance(*, permute = true, scale = true, separate = false)
      a = to_general
      s = a.balance!(permute: permute, scale: scale, separate: separate)
      {a, s}
    end

    # Computes the Hessenberg form of a matrix.
    #
    # Arguments:
    #   - calc_q (Bool) : If `true`, also computes the orthogonal matrix Q. Default: `false`.
    #
    # Returns:
    #   - Tuple(GeneralMatrix(T), GeneralMatrix(T)) : Hessenberg matrix and optionally Q.
    #
    # See `GeneralMatrix#hessenberg!` for details.
    def hessenberg(*, calc_q = false)
      x = to_general
      x.hessenberg!(calc_q: calc_q)
    end

    # Computes the Hessenberg form in-place.
    #
    # Returns:
    #   - GeneralMatrix(T) : Hessenberg matrix.
    def hessenberg
      to_general.hessenberg!
    end

    # Calculate Moore–Penrose pseudo-inverse of a matrix.
    #
    # https://en.wikipedia.org/wiki/Moore%E2%80%93Penrose_inverse
    #
    # Implemented using SVD decomposition.
    #
    # Returns:
    #   - GeneralMatrix(T) : Pseudo-inverse matrix.
    def pinv
      # Pure Crystal implementation since LAPACK has no direct implementation of pseudo inverse
      u, s, v = self.svd
      v.conjtranspose!
      u.conjtranspose!

      s_inverse = s.map! { |e|
        if e != 0
          1.0 / e
        else
          e
        end
      }
      s_dash = self.class.diag(s_inverse)

      append_rows = 0
      append_cols = 0

      if v.ncolumns >= s_dash.nrows
        append_rows = v.ncolumns - s_dash.nrows
      else
        puts "S` #{s_dash}"
        puts "v #{v}"
        raise Exception.new("Invalid dimension, S` larger than v")
      end

      if u.nrows >= s_dash.ncolumns
        append_cols = u.nrows - s_dash.ncolumns
      else
        puts "S` #{s_dash}\nShape #{s_dash.size}"
        puts "U #{u}\nShape #{u.size}"
        raise Exception.new("Invalid dimension, S` larger than U")
      end

      s_dash.resize!(s_dash.nrows + append_rows, s_dash.ncolumns + append_cols)
      return v * s_dash * u
    end

    # Computes the norm of the matrix.
    #
    # Arguments:
    #   - kind (MatrixNorm) : Type of norm to compute. Default: `MatrixNorm::Frobenius`.
    #
    # Returns:
    #   - Float : Norm value.
    def norm(kind : MatrixNorm = MatrixNorm::Frobenius)
      result = of_real_type(0)
      case kind
      in .frobenius?
        each { |v| result += of_real_type(v*v) }
        Math.sqrt(result)
      in .one?
        sums = of_real_type(Slice, ncolumns)
        each_with_index do |v, row, col|
          sums[col] += v.abs
        end
        sums.max
      in .inf?
        sums = of_real_type(Slice, nrows)
        each_with_index do |v, row, col|
          sums[row] += v.abs
        end
        sums.max
      in .max_abs?
        each { |v| result = v.abs if result < v.abs }
        result
      end
    end
  end

  class GeneralMatrix(T) < Matrix(T)
    # Calculate matrix inversion in-place.
    #
    # Method selects optimal algorithm depending on `MatrixFlags`:
    #   - `transpose!` returned if matrix is orthogonal.
    #   - `trtri` used if matrix is triangular.
    #   - `potrf` + `potri` used if matrix is positive definite.
    #   - `hetrf` + `hetri` used if matrix is hermitian (complex only).
    #   - `sytrf` + `sytri` used if matrix is symmetric.
    #   - `getrf` + `getri` used otherwise.
    #
    # Returns:
    #   - self : On successful inversion, `self` becomes the inverse.
    #
    # Raises:
    #   - ArgumentError : If matrix is not square.
    #
    # LAPACK Routines:
    #   - `trtri` (triangular inverse)
    #   - `potrf`/`potri` (positive definite)
    #   - `hetrf`/`hetri` (hermitian)
    #   - `sytrf`/`sytri` (symmetric)
    #   - `getrf`/`getri` (general)
    def inv!
      raise ArgumentError.new("can't invert nonsquare matrix") unless square?
      return transpose! if flags.orthogonal?
      n = self.nrows
      if flags.triangular?
        lapack(trtri, uplo, 'N'.ord.to_u8, n, self, n)
        adjust_triangular
      elsif flags.positive_definite?
        lapack(potrf, uplo, n, self, n)
        lapack(potri, uplo, n, self, n)
        adjust_symmetric
      elsif {{T == Complex}} && flags.hermitian?
        {% if T == Complex %}
          ipiv = Slice(Int32).new(n)
          lapack(hetrf, uplo, n, self, n, ipiv)
          lapack(hetri, uplo, n, self, n, ipiv, worksize: [n])
          adjust_symmetric
        {% else %}
          raise "error" # to prevent type inference of nil
        {% end %}
      elsif flags.symmetric?
        ipiv = Slice(Int32).new(n)
        lapack(sytrf, uplo, n, self, n, ipiv)
        lapack(sytri, uplo, n, self, n, ipiv, worksize: [2*n])
        adjust_symmetric
      else
        ipiv = Slice(Int32).new(n)
        lapack(getrf, n, n, self, n, ipiv)
        lapack(getri, n, self, n, ipiv)
      end
      self
    end

    # Solves the linear system A * X = B.
    #
    # Method selects optimal algorithm depending on `MatrixFlags`:
    #   - `trtrs` used if matrix is triangular.
    #   - `posv` used if matrix is positive definite.
    #   - `hesv` used if matrix is hermitian (complex only).
    #   - `sysv` used if matrix is symmetric.
    #   - `gesv` used otherwise.
    #
    # Arguments:
    #   - b (GeneralMatrix(T)) : Right-hand side matrix.
    #   - overwrite_a (Bool) : If `true`, allows overwriting `self`. Default: `false`.
    #   - overwrite_b (Bool) : If `true`, allows overwriting `b`. Default: `false`.
    #
    # Returns:
    #   - GeneralMatrix(T) : Solution matrix X.
    #
    # Raises:
    #   - ArgumentError : If `nrows` mismatch or matrix not square.
    #
    # LAPACK Routines:
    #   - `trtrs` (triangular solve)
    #   - `posv` (positive definite solve)
    #   - `hesv` (hermitian solve)
    #   - `sysv` (symmetric solve)
    #   - `gesv` (general solve)
    def solve(b : self, *, overwrite_a = false, overwrite_b = false)
      raise ArgumentError.new("nrows of a and b must match") unless nrows == b.nrows
      raise ArgumentError.new("a must be square") unless square?
      a = overwrite_b ? self : self.clone
      x = overwrite_b ? b : b.clone
      n = nrows
      if flags.triangular?
        lapack(trtrs, uplo, 'N'.ord.to_u8, 'N'.ord.to_u8, n, b.nrows, a, n, x, b.nrows)
      elsif flags.positive_definite?
        lapack(posv, 'U'.ord.to_u8, n, b.ncolumns, a, n, x, b.nrows)
      elsif flags.hermitian?
        {% if T == Complex %}
          ipiv = Slice(Int32).new(n)
          lapack(hesv, uplo, n, b.ncolumns, a, n, ipiv, x, b.nrows)
        {% end %}
      elsif flags.symmetric?
        ipiv = Slice(Int32).new(n)
        lapack(sysv, uplo, n, b.ncolumns, a, n, ipiv, x, b.nrows)
      else
        ipiv = Slice(Int32).new(n)
        lapack(gesv, n, b.ncolumns, a, n, ipiv, x, b.nrows)
      end
      a.clear_flags
      x.clear_flags
      x
    end

    # Calculates the determinant of a square matrix.
    #
    # Arguments:
    #   - overwrite_a (Bool) : If `true`, allows overwriting `self`. Default: `false`.
    #
    # Returns:
    #   - T : Determinant value.
    #
    # Raises:
    #   - ArgumentError : If matrix is not square.
    #
    # LAPACK Routine:
    #   - Uses `getrf` (LU factorization).
    def det(*, overwrite_a = false)
      raise ArgumentError.new("matrix must be square") unless square?
      if flags.triangular?
        return diag.product
      end
      lru = overwrite_a ? self : self.clone
      ipiv = Slice(Int32).new(nrows)
      lapack(getrf, nrows, nrows, lru, nrows, ipiv)
      lru.clear_flags
      lru.diag.product
    end

    # Solves the linear least squares problem using QR factorization.
    #
    # Arguments:
    #   - b (GeneralMatrix(T)) : Right-hand side matrix.
    #   - overwrite_a (Bool) : If `true`, allows overwriting `self`. Default: `false`.
    #   - overwrite_b (Bool) : If `true`, allows overwriting `b`. Default: `false`.
    #   - cond (Number) : Cutoff for rank determination. Default: -1 (machine precision).
    #
    # Returns:
    #   - GeneralMatrix(T) : Solution matrix X (may include residuals).
    #
    # Raises:
    #   - ArgumentError : If `nrows` mismatch.
    #
    # LAPACK Routine:
    #   - Uses `gels` (QR least squares).
    def solvels(b : self, *, overwrite_a = false, overwrite_b = false, cond = -1)
      raise ArgumentError.new("nrows of a and b must match") unless nrows == b.nrows
      a = overwrite_a ? self : self.clone
      if ncolumns > nrows
        # make room for residuals
        x = GeneralMatrix(T).new(ncolumns, b.ncolumns) { |r, c| r < nrows ? b.unsafe_fetch(r, c) : T.new(0) }
      else
        x = overwrite_b ? b : b.clone
      end
      lapack(gels, 'N'.ord.to_u8, nrows, ncolumns, b.ncolumns, a, nrows, x, x.nrows)
      a.clear_flags
      x.clear_flags
      x
    end

    # Solves the linear least squares problem with multiple method options.
    #
    # Arguments:
    #   - b (GeneralMatrix(T)) : Right-hand side matrix.
    #   - method (LSMethod) : Algorithm to use. Default: `LSMethod::Auto`.
    #   - overwrite_a (Bool) : If `true`, allows overwriting `self`. Default: `false`.
    #   - overwrite_b (Bool) : If `true`, allows overwriting `b`. Default: `false`.
    #   - cond (Number) : Cutoff for rank determination. Default: -1 (machine precision).
    #
    # Returns:
    #   - Tuple(GeneralMatrix(T), Int32, Array(T)) : Solution X, effective rank, singular values (if applicable).
    #
    # Raises:
    #   - ArgumentError : If `nrows` mismatch.
    #
    # LAPACK Routines:
    #   - `gels` (QR least squares)
    #   - `gelsd` (SVD-based least squares)
    #   - `gelsy` (QR with pivoting)
    def lstsq(b : self, method : LSMethod = LSMethod::Auto, *, overwrite_a = false, overwrite_b = false, cond = -1)
      raise ArgumentError.new("nrows of a and b must match") unless nrows == b.nrows
      if method.auto?
        method = LSMethod::QR
      end
      a = overwrite_a ? self : self.clone
      if ncolumns > nrows
        # make room for residuals
        x = GeneralMatrix(T).new(ncolumns, b.ncolumns) { |r, c| r < nrows ? b.unsafe_fetch(r, c) : T.new(0) }
      else
        x = overwrite_b ? b : b.clone
      end
      rank = 0
      case method
      when .ls?
        lapack(gels, 'N'.ord.to_u8, nrows, ncolumns, b.ncolumns, a, nrows, x, x.nrows)
        s = of_real_type(Array, 0)
      when .lsd?
        ssize = {nrows, ncolumns}.min
        s = of_real_type(Array, ssize)
        rcond = of_real_type(cond)
        lapack(gelsd, nrows, ncolumns, b.ncolumns, a, nrows, x, x.nrows, s, rcond, rank)
      when .lsy?
        jpvt = Slice(Int32).new(ncolumns)
        rcond = of_real_type(cond)
        lapack(gelsy, nrows, ncolumns, b.ncolumns, a, nrows, x, x.nrows, jpvt, rcond, rank, worksize: [2*ncolumns])
        s = of_real_type(Array, 0)
      else
        s = of_real_type(Array, 0)
      end
      a.clear_flags
      x.clear_flags
      {x, rank, s}
    end

    # Computes the singular value decomposition (SVD) of a matrix.
    #
    # For an m×n matrix A, returns:
    #   - u : m×m orthogonal matrix
    #   - s : Array of singular values of length min(m,n)
    #   - vt : n×n orthogonal matrix
    #
    # Such that A = u * diag(s) * vt.
    #
    # Arguments:
    #   - overwrite_a (Bool) : If `true`, allows overwriting `self`. Default: `false`.
    #
    # Returns:
    #   - Tuple(GeneralMatrix(T), Array(T), GeneralMatrix(T)) : U, singular values, V^T.
    #
    # LAPACK Routine:
    #   - Uses `gesdd` (divide-and-conquer SVD).
    def svd(*, overwrite_a = false)
      a = overwrite_a ? self : self.clone
      m = nrows
      n = ncolumns
      mn = {m, n}.min
      mx = {m, n}.max
      s = of_real_type(Array, mn)
      u = GeneralMatrix(T).new(m, m)
      vt = GeneralMatrix(T).new(n, n)
      lapack(gesdd, 'A'.ord.to_u8, m, n, a, nrows, s, u, m, vt, n, worksize: [{5*mn*mn + 5*mn, 2*mx*mn + 2*mn*mn + mn}.max, 8*mn])
      a.clear_flags
      return {u, s, vt}
    end

    # Computes only the singular values of a matrix.
    #
    # Arguments:
    #   - overwrite_a (Bool) : If `true`, allows overwriting `self`. Default: `false`.
    #
    # Returns:
    #   - Array(T) : Singular values.
    #
    # LAPACK Routine:
    #   - Uses `gesdd` (SVD, singular values only).
    def svdvals(*, overwrite_a = false)
      a = overwrite_a ? self : self.clone
      m = nrows
      n = ncolumns
      mn = {m, n}.min
      mx = {m, n}.max
      s = of_real_type(Array, mn)
      lapack(gesdd, 'N'.ord.to_u8, m, n, a, nrows, s, nil, m, nil, n, worksize: [5*mn, 8*mn])
      a.clear_flags
      s
    end

    # Balances a square matrix in-place to improve eigenvalue accuracy.
    #
    # Arguments:
    #   - permute (Bool) : If `true`, performs permutations. Default: `true`.
    #   - scale (Bool) : If `true`, performs scaling. Default: `true`.
    #   - separate (Bool) : If `true`, returns scaling factors separately. Default: `false`.
    #
    # Returns:
    #   - GeneralMatrix(T) : If `separate` is true, returns scaling factors; otherwise, returns diagonal matrix.
    #
    # Raises:
    #   - ArgumentError : If matrix is not square.
    #
    # LAPACK Routine:
    #   - Uses `gebal` (balance matrix).
    def balance!(*, permute = true, scale = true, separate = false)
      raise ArgumentError.new("matrix must be square") unless square?
      n = self.nrows
      job = if permute && scale
              'B'
            elsif permute
              'P'
            elsif scale
              'S'
            else
              # don't call anything, return identity matrix
              return separate ? Matrix(T).ones(1, n) : Matrix(T).identity(n)
            end
      s = GeneralMatrix(T).new(1, n)
      ilo = 0
      ihi = 0
      lapack(gebal, job.ord.to_u8, n, self, n, ilo, ihi, s)
      separate ? s : Matrix(T).diag(s.raw)
    end

    # Computes the Hessenberg form of a matrix in-place.
    #
    # Arguments:
    #   - calc_q (Bool) : If `true`, also computes the orthogonal matrix Q. Default: `false`.
    #
    # Returns:
    #   - Tuple(GeneralMatrix(T), GeneralMatrix(T)) : Hessenberg matrix and optionally Q.
    #
    # Raises:
    #   - ArgumentError : If matrix is not square.
    #
    # LAPACK Routines:
    #   - `gebal` (balance)
    #   - `gehrd` (Hessenberg reduction)
    #   - `orghr` (generate Q)
    def hessenberg!(*, calc_q = false)
      raise ArgumentError.new("matrix must be square") unless square?
      # idea from scipy.
      # no need to calculate if size <= 2
      if nrows < 2
        q = calc_q ? Matrix(T).identity(nrows) : Matrix(T).zeros(1, 1)
        return {self, q}
      end
      {% if flag?(:darwin) %}
        raise "Hessenberg decomposition is not supported on mac"
      {% end %}
      n = nrows
      s = of_real_type(Slice, n)
      lapack(gebal, 'S'.ord.to_u8, n, self, n, ilo, ihi, s)
      clear_flags
      tau = GeneralMatrix(T).new(1, n)
      lapack(gehrd, n, ilo, ihi, self, ncolumns, tau)
      if calc_q
        q = clone
        lapack(orghr, n, ilo, ihi, q, ncolumns, tau)
        q.flags = MatrixFlags::Orthogonal
      else
        q = Matrix(T).zeros(1, 1)
      end
      triu!(-1)
      {self, q}
    end

    # Computes the Hessenberg form in-place (without Q).
    #
    # Returns:
    #   - self : Hessenberg matrix.
    def hessenberg!
      q = hessenberg!(calc_q: false)
      self
    end

    # Computes the norm of the matrix.
    #
    # Method selects optimal algorithm depending on `MatrixFlags`:
    #   - `lantr` used for triangular matrices.
    #   - `lanhe` used for hermitian matrices (complex only).
    #   - `lansy` used for symmetric matrices.
    #   - `lange` used for general matrices.
    #
    # Arguments:
    #   - kind (MatrixNorm) : Type of norm to compute. Default: `MatrixNorm::Frobenius`.
    #
    # Returns:
    #   - T : Norm value.
    #
    # LAPACK Routines:
    #   - `lantr` (triangular matrix norm)
    #   - `lanhe` (hermitian matrix norm)
    #   - `lansy` (symmetric matrix norm)
    #   - `lange` (general matrix norm)
    def norm(kind : MatrixNorm = MatrixNorm::Frobenius)
      let = case kind
            in .frobenius?
              'F'
            in .one?
              'O'
            in .inf?
              'I'
            in .max_abs?
              'M'
            end.ord.to_u8

      worksize = (kind.inf? || kind.one?) ? nrows : 0

      if flags.triangular?
        return lapack_util(lantr, worksize, let, uplo, 'N'.ord.to_u8, @nrows, @ncolumns, matrix(self), @nrows)
      end
      {% if T == Complex %}
        if flags.hermitian?
          return lapack_util(lanhe, worksize, let, uplo, @nrows, matrix(self), @nrows)
        end
      {% end %}
      if flags.symmetric?
        lapack_util(lansy, worksize, let, uplo, @nrows, matrix(self), @nrows)
      else
        lapack_util(lange, worksize, let, @nrows, @ncolumns, matrix(self), @nrows)
      end
    end

    # Alias for `#norm`.
    def abs(kind : MatrixNorm = MatrixNorm::Frobenius)
      norm(kind)
    end

    # Determines the effective rank of the matrix.
    #
    # Arguments:
    #   - eps (Number) : Tolerance for singular value comparison. Default: `self.tolerance`.
    #   - method (RankMethod) : Method to use. Default: `RankMethod::SVD`.
    #   - overwrite_a (Bool) : If `true`, allows overwriting `self`. Default: `false`.
    #
    # Returns:
    #   - Int32 : Estimated rank.
    #
    # LAPACK Routines:
    #   - SVD method uses `gesdd` (singular values)
    #   - QRP method uses QR with pivoting
    def rank(eps = self.tolerance, *, method : RankMethod = RankMethod::SVD, overwrite_a = false)
      # if matrix is triangular no check needed
      return diag.count { |v| v.abs > eps } if flags.triangular?
      case method
      in .qrp?
        a, pvt = qr_r(overwrite_a: overwrite_a, pivoting: true)
        a.diag.count { |v| v.abs > eps }
      in .svd?
        s = svdvals(overwrite_a: overwrite_a)
        s.count { |x| x.abs > eps }
      end
    end
  end
end
