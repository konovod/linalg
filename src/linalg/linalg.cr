require "../matrix/*"
require "./libLAPACKE"

module Linalg
  enum LSMethod
    Auto       = 0
    QR
    Orthogonal
    SVD
    LS         = QR
    LSY        = Orthogonal
    LSD        = SVD
  end

  enum MatrixNorm
    Frobenius
    One
    # Two
    Inf
    MaxAbs
  end

  class LinAlgError < Exception
  end

  def self.inv(matrix, *, overwrite_a = false)
    overwrite_a ? matrix.inv! : matrix.inv
  end

  def self.solve(a, b, *, overwrite_a = false, overwrite_b = false)
    a.solve(b, overwrite_a: overwrite_a, overwrite_b: overwrite_b)
  end

  def self.lstsq(a, b, method : LSMethod = LSMethod::Auto, *, overwrite_a = false, overwrite_b = false, cond = -1)
    a.lstsq(b, method, overwrite_a: overwrite_a, overwrite_b: overwrite_b, cond: cond)
  end

  def self.solvels(a, b, *, overwrite_a = false, overwrite_b = false, cond = -1)
    a.solvels(b, overwrite_a: overwrite_a, overwrite_b: overwrite_b, cond: cond)
  end

  def self.svd(matrix, *, overwrite_a = false)
    matrix.svd(overwrite_a: overwrite_a)
  end

  module Matrix(T)
    macro lapack_util(name, *args)
      {% if T == Float32
           typ = :s.id
         elsif T == Float64
           typ = :d.id
         elsif T == Complex
           typ = :z.id
         end %}
       LibLAPACKE.{{typ}}{{name}}(LibLAPACKE::ROW_MAJOR, {{*args}})
    end

    macro lapack(storage, name, *args)
      {% if T == Float32
           typ = :s.id
         elsif T == Float64
           typ = :d.id
         elsif T == Complex
           typ = :z.id
         end %}
       {% if T == Complex && storage.id == :or.id
            st = :un.id
          else
            st = storage
          end %}
       info = LibLAPACKE.{{typ}}{{st}}{{name}}(LibLAPACKE::ROW_MAJOR, {{*args}})
       raise LinAlgError.new("LAPACKE.{{typ}}{{storage}}{{name}} returned #{info}") if info != 0
    end

    private def uplo
      flags.lower_triangular? ? 'L'.ord : 'U'.ord
    end

    def inv!
      raise ArgumentError.new("can't invert nonsquare matrix") unless square?
      return transpose! if flags.orthogonal?
      n = @nrows
      if flags.positive_definite?
        lapack(po, trf, uplo, n, self, n)
        lapack(po, tri, uplo, n, self, n)
      elsif flags.hermitian?
        {% if T == Complex %}
        ipiv = Slice(Int32).new(n)
        lapack(he, trf, uplo, n, self, n, ipiv)
        lapack(he, tri, uplo, n, self, n, ipiv)
        {% end %}
      elsif flags.symmetric?
        ipiv = Slice(Int32).new(n)
        lapack(sy, trf, uplo, n, self, n, ipiv)
        lapack(sy, tri, uplo, n, self, n, ipiv)
      elsif flags.triangular?
        lapack(tr, tri, uplo, 'N'.ord, n, self, n)
      else
        ipiv = Slice(Int32).new(n)
        lapack(ge, trf, n, n, self, n, ipiv)
        lapack(ge, tri, n, self, n, ipiv)
      end
      self
    end

    def inv
      clone.inv!
    end

    def solve(b : self, *, overwrite_a = false, overwrite_b = false)
      raise ArgumentError.new("nrows of a and b must match") unless nrows == b.nrows
      raise ArgumentError.new("a must be square") unless square?
      a = overwrite_b ? self : self.clone
      x = overwrite_b ? b : b.clone
      n = nrows
      if flags.triangular?
        lapack(tr, trs, uplo, 'N'.ord, 'N'.ord, n, b.ncolumns, a, n, x, b.ncolumns)
      elsif flags.positive_definite?
        lapack(po, sv, uplo, n, b.ncolumns, a, n, x, b.ncolumns)
      elsif flags.hermitian?
        {% if T == Complex %}
        ipiv = Slice(Int32).new(n)
        lapack(he, sv, uplo, n, b.ncolumns, a, n, ipiv, x, b.ncolumns)
        {% end %}
      elsif flags.symmetric?
        ipiv = Slice(Int32).new(n)
        lapack(sy, sv, uplo, n, b.ncolumns, a, n, ipiv, x, b.ncolumns)
      else
        ipiv = Slice(Int32).new(n)
        lapack(ge, sv, n, b.ncolumns, a, n, ipiv, x, b.ncolumns)
      end
      a.clear_flags
      x.clear_flags
      x
    end

    def det(*, overwrite_a = false)
      raise ArgumentError.new("matrix must be square") unless square?
      lru = overwrite_a ? self : self.clone
      ipiv = Slice(Int32).new(nrows)
      lapack(ge, trf, nrows, nrows, lru, nrows, ipiv)
      lru.clear_flags
      (0...nrows).reduce(1) { |det, i| det*lru[i, i] }
    end

    def solvels(b : self, *, overwrite_a = false, overwrite_b = false, cond = -1)
      raise ArgumentError.new("nrows of a and b must match") unless nrows == b.nrows
      a = overwrite_a ? self : self.clone
      if ncolumns > nrows
        # make room for residuals
        x = GeneralMatrix(T).new(ncolumns, b.ncolumns) { |r, c| r < nrows ? b[r, c] : T.new(0) }
      else
        x = overwrite_b ? b : b.clone
      end
      lapack(ge, ls, 'N'.ord, nrows, ncolumns, b.ncolumns, a, ncolumns, x, x.ncolumns)
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
        x = GeneralMatrix(T).new(ncolumns, b.ncolumns) { |r, c| r < nrows ? b[r, c] : T.new(0) }
      else
        x = overwrite_b ? b : b.clone
      end
      rank = 0
      case method
      when .ls?
        lapack(ge, ls, 'N'.ord, nrows, ncolumns, b.ncolumns, a, ncolumns, x, x.ncolumns)
        s = {% if T == Complex %} Array(Float64) {% else %} Array(T) {% end %}.new
      when .lsd?
        ssize = {nrows, ncolumns}.min
        s = {% if T == Complex %} Array(Float64).new(ssize, 0.0) {% else %} Array(T).new(ssize, T.new(0)) {% end %}
        rcond = {% if T == Complex %} Float64 {% else %} T {% end %}.new(cond)
        lapack(ge, lsd, nrows, ncolumns, b.ncolumns, a, ncolumns, x, x.ncolumns, s, rcond, pointerof(rank))
      when .lsy?
        jpvt = Slice(Int32).new(ncolumns)
        rcond = {% if T == Complex %} Float64 {% else %} T {% end %}.new(cond)
        lapack(ge, lsy, nrows, ncolumns, b.ncolumns, a, ncolumns, x, x.ncolumns, jpvt, rcond, pointerof(rank))
        s = {% if T == Complex %} Array(Float64) {% else %} Array(T) {% end %}.new
      else
        s = {% if T == Complex %} Array(Float64) {% else %} Array(T) {% end %}.new
      end
      a.clear_flags
      x.clear_flags
      {x, rank, s}
    end

    def svd(*, overwrite_a = false)
      a = overwrite_a ? self : self.clone
      m = nrows
      n = ncolumns
      s = {% if T == Complex %} Array(Float64).new({m, n}.min, 0.0) {% else %} Array(T).new({m, n}.min, T.new(0)) {% end %}
      u = GeneralMatrix(T).new(m, m)
      vt = GeneralMatrix(T).new(n, n)
      lapack(ge, sdd, 'A'.ord, m, n, a, ncolumns, s, u, m, vt, n)
      a.clear_flags
      return {u, s, vt}
    end

    def svdvals(*, overwrite_a = false)
      a = overwrite_a ? self : self.clone
      m = nrows
      n = ncolumns
      s = {% if T == Complex %} Array(Float64).new({m, n}.min, 0.0) {% else %} Array(T).new({m, n}.min, T.new(0)) {% end %}
      lapack(ge, sdd, 'N'.ord, m, n, a, ncolumns, s, nil, m, nil, n)
      a.clear_flags
      s
    end

    def balance!(*, permute = true, scale = true, separate = false)
      raise ArgumentError.new("matrix must be square") unless square?
      n = @nrows
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
      lapack(ge, bal, job.ord, n, self, n, out ilo, out ihi, s)
      separate ? s : Matrix(T).diag(s.raw)
    end

    def balance(*, permute = true, scale = true, separate = false)
      a = clone
      s = a.balance!(permute: permute, scale: scale, separate: separate)
      {a, s}
    end

    def hessenberg!(*, calc_q = false)
      raise ArgumentError.new("matrix must be square") unless square?
      # idea from scipy.
      # no need to calculate if size <= 2
      if nrows < 2
        q = calc_q ? Matrix(T).identity(nrows) : Matrix(T).zeros(1, 1)
        return {self, q}
      end
      n = nrows
      s = {% if T == Complex %} Slice(Float64) {% else %} Slice(T) {% end %}.new(n)
      lapack(ge, bal, 'S'.ord, n, self, n, out ilo, out ihi, s)
      clear_flags
      tau = GeneralMatrix(T).new(1, n)
      lapack(ge, hrd, n, ilo, ihi, self, ncolumns, tau)
      if calc_q
        q = clone
        lapack(or, ghr, n, ilo, ihi, q, ncolumns, tau)
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

    def hessenberg(*, calc_q = false)
      x = self.clone
      x.hessenberg!(calc_q: calc_q)
    end

    def hessenberg
      clone.hessenberg!
    end

    # returns matrix norm
    def norm(kind : MatrixNorm = MatrixNorm::Frobenius)
      let = case kind
            when .frobenius?
              'F'
            when .one?
              'O'
            when .inf?
              'I'
            else
              'M'
            end.ord
      if flags.triangular?
        lapack_util(lantr, let, uplo, 'N'.ord, nrows, ncolumns, self, ncolumns)
      elsif flags.hermitian?
        {% if T == Complex %}
        lapack_util(lanhe, let, uplo, nrows,  self, ncolumns)
        {% else %}
        lapack_util(lange, let, nrows, ncolumns, self, ncolumns)
        {% end %}
      elsif flags.symmetric?
        lapack_util(lansy, let, uplo, nrows, self, ncolumns)
      else
        lapack_util(lange, let, nrows, ncolumns, self, ncolumns)
      end
    end

    def abs(kind : MatrixNorm = MatrixNorm::Frobenius)
      norm(kind)
    end
  end
end
