require "../matrix/*"
require "./lapack_helper"

# TODO - inline docs

module LA
  module Enums
    enum LSMethod
      Auto       = 0
      QR
      Orthogonal
      SVD
      LS         = QR
      LSY        = Orthogonal
      LSD        = SVD
    end

    enum RankMethod
      SVD
      QRP
    end

    enum MatrixNorm
      Frobenius
      One
      # Two
      Inf
      MaxAbs
    end
  end

  class Utils::LinAlgError < Exception
  end

  def self.lapack_version
    LibLAPACK.ilaver(out major, out minor, out patch)
    {major, minor, patch}
  end

  # Calculate matrix inversion
  #
  # if `overwrite_a` is true, source matrix isn't needed anymore and can be overriden in process
  # See `#inv!` for details of algorithm
  def self.inv(matrix, *, overwrite_a = false)
    overwrite_a ? matrix.inv! : matrix.inv
  end

  # See `#solve`
  # `Matrix.solve(a,b)` is an alias for `a.solve(b)`
  def self.solve(a, b, *, overwrite_a = false, overwrite_b = false)
    a.solve(b, overwrite_a: overwrite_a, overwrite_b: overwrite_b)
  end

  # See `#lstsq`
  # `Matrix.lstsq(a,b)` is an alias for `a.lstsq(b)`
  def self.lstsq(a, b, method : LSMethod = LSMethod::Auto, *, overwrite_a = false, overwrite_b = false, cond = -1)
    a.lstsq(b, method, overwrite_a: overwrite_a, overwrite_b: overwrite_b, cond: cond)
  end

  # See `#solvels`
  # `Matrix.solvels(a,b)` is an alias for `a.solvels(b)`
  def self.solvels(a, b, *, overwrite_a = false, overwrite_b = false, cond = -1)
    a.solvels(b, overwrite_a: overwrite_a, overwrite_b: overwrite_b, cond: cond)
  end

  # See `#svd`
  # `Matrix.svd(a)` is an alias for `a.svd`
  def self.svd(matrix, *, overwrite_a = false)
    matrix.svd(overwrite_a: overwrite_a)
  end

  abstract class Matrix(T)
    private def uplo
      flags.lower_triangular? ? 'L'.ord.to_u8 : 'U'.ord.to_u8
    end

    private def adjust_symmetric
      f = flags
      each_upper(diagonal: false) { |v, i, j| unsafe_set(j, i, v) }
      self.flags = f
    end

    private def adjust_triangular
      triu! if flags.upper_triangular?
      tril! if flags.lower_triangular?
    end

    # Calculate matrix inversion
    # See `#inv!` for details on algorithm
    def inv
      to_general.inv!
    end

    def svd
      to_general.svd(overwrite_a: true)
    end

    def svdvals
      to_general.svdvals(overwrite_a: true)
    end

    def solve(b : GeneralMatrix(T), *, overwrite_b = false)
      to_general.solve(b, overwrite_b: overwrite_b)
    end

    def solve(b : self)
      to_general.solve(b.to_general, overwrite_a: true, overwrite_b: true)
    end

    def det
      to_general.det(overwrite_a: true)
    end

    def balance(*, permute = true, scale = true, separate = false)
      a = to_general
      s = a.balance!(permute: permute, scale: scale, separate: separate)
      {a, s}
    end

    def hessenberg(*, calc_q = false)
      x = to_general
      x.hessenberg!(calc_q: calc_q)
    end

    def hessenberg
      to_general.hessenberg!
    end

    # Calculate Moore–Penrose inverse of a matrix
    #
    # https://en.wikipedia.org/wiki/Moore%E2%80%93Penrose_inverse
    # Implemented using an svd decomposition
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
    # Calculate matrix inversion inplace
    # Method selects optimal algorithm depending on `MatrixFlags`
    # `transpose` returned if matrix is orthogonal
    # `trtri` is used if matrix is triangular
    # `potrf`, `potri` are used if matrix is positive definite
    # `hetrf`, `hetri` are used if matrix is hermitian
    # `sytrf`, `sytri` are used if matrix is symmetric
    # `getrf`, `getri` are used otherwise
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

    # Solves matrix equation `self*x = b` and returns x
    #
    # Method returns matrix of same size as `b`
    #
    # Matrix must be square, number of rows must match `b`
    #
    # if overwrite_a is true, a will be overriden in process of calculation
    #
    # if overwrite_b is true, b will be overriden in process of calculation
    #
    # Uses LAPACK routines `trtrs`, `posv`, `hesv`, `sysv`, `gesv` depending on matrix `flags`
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

    # Calculates determinant for a square matrix
    #
    # if overwrite_a is true, a will be overriden in process of calculation
    #
    # Uses `getrf` LAPACK routine
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

    # Calculates singular value decomposition for a matrix
    #
    # If you call `u,s,vt = a.svd` for a matrix m*n `a` then
    #
    #  - `u` will be m*m matrix
    #  - `s` will be Array(T) with size equal to `{m, n}.min`
    #  - `vt` will be n*n matrix
    #  - `(u*Mat.diag(a.nrows, a.ncolumns, s)*vt)` will be equal to `a` (within calculation tolerance)
    #
    # if overwrite_a is true, a will be overriden in process of calculation
    #
    # Uses `gesdd` LAPACK routine
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

    # Calculates array of singular values for a matrix
    #
    # if `overwrite_a` is true, a will be overriden in process of calculation
    #
    # Uses `gesdd` LAPACK routine
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

    def hessenberg!
      q = hessenberg!(calc_q: false)
      self
    end

    # returns matrix norm
    #
    # `kind` defines type of norm
    #
    # Uses LAPACK routines `lantr`, `lanhe`, `lansy`, `lange` depending on matrix `flags`
    #
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

    # Alias for `#norm`
    def abs(kind : MatrixNorm = MatrixNorm::Frobenius)
      norm(kind)
    end

    # determine effective rank either by SVD method or QR-factorization with pivoting
    #
    # if overwrite_a is true, a will be overriden in process of calculation
    #
    # QR method is faster, but could fail to determine rank in some cases
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
