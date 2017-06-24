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

  class LinAlgError < Exception
  end

  def self.inv(matrix, *, overwrite_a = false)
    matrix.inv(overwrite_a: overwrite_a)
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
    macro lapack(storage, name, *args)
      {% if T == Float32
           typ = :s.id
         elsif T == Float64
           typ = :d.id
         elsif T == Complex
           typ = :z.id
         end %}
       info = LibLAPACKE.{{typ}}{{storage}}{{name}}(LibLAPACKE::ROW_MAJOR, {{*args}})
       raise LinAlgError.new("LAPACKE.{{typ}}{{storage}}{{name}} returned #{info}") if info != 0
    end

    private def uplo
      flags.lower? ? 'L'.ord : 'U'.ord
    end

    def inv!
      raise ArgumentError.new("can't invert nonsquare matrix") unless square?
      n = @rows
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
      raise ArgumentError.new("number of rows in a and b must match") unless rows == b.rows
      raise ArgumentError.new("a must be square") unless square?
      a = overwrite_b ? self : self.clone
      x = overwrite_b ? b : b.clone
      n = rows
      # TODO - gb, gt, pb, pt
      if flags.positive_definite?
        lapack(po, sv, uplo, n, b.columns, a, n, x, b.columns)
      elsif flags.hermitian?
        {% if T == Complex %}
        ipiv = Slice(Int32).new(n)
        lapack(he, sv, uplo, n, b.columns, a, n, ipiv, x, b.columns)
        {% end %}
      elsif flags.symmetric?
        ipiv = Slice(Int32).new(n)
        lapack(sy, sv, uplo, n, b.columns, a, n, ipiv, x, b.columns)
      else
        ipiv = Slice(Int32).new(n)
        lapack(ge, sv, n, b.columns, a, n, ipiv, x, b.columns)
      end
      x
    end

    def det(*, overwrite_a = false)
      raise ArgumentError.new("matrix must be square") unless square?
      lru = overwrite_a ? self : self.clone
      ipiv = Slice(Int32).new(rows)
      lapack(ge, trf, rows, rows, lru, rows, ipiv)
      (0...rows).reduce(1) { |det, i| det*lru[i, i] }
    end

    def solvels(b : self, *, overwrite_a = false, overwrite_b = false, cond = -1)
      raise ArgumentError.new("number of rows in a and b must match") unless rows == b.rows
      a = overwrite_a ? self : self.clone
      if columns > rows
        # make room for residuals
        x = GeneralMatrix(T).new(columns, b.columns) { |r, c| r < rows ? b[r, c] : T.new(0) }
      else
        x = overwrite_b ? b : b.clone
      end
      lapack(ge, ls, 'N'.ord, rows, columns, b.columns, a, columns, x, x.columns)
      x
    end

    def lstsq(b : self, method : LSMethod = LSMethod::Auto, *, overwrite_a = false, overwrite_b = false, cond = -1)
      raise ArgumentError.new("number of rows in a and b must match") unless rows == b.rows
      if method.auto?
        method = LSMethod::QR
      end
      a = overwrite_a ? self : self.clone
      if columns > rows
        # make room for residuals
        x = GeneralMatrix(T).new(columns, b.columns) { |r, c| r < rows ? b[r, c] : T.new(0) }
      else
        x = overwrite_b ? b : b.clone
      end
      rank = 0
      case method
      when .ls?
        lapack(ge, ls, 'N'.ord, rows, columns, b.columns, a, columns, x, x.columns)
        s = {% if T == Complex %} Array(Float64) {% else %} Array(T) {% end %}.new
      when .lsd?
        ssize = {rows, columns}.min
        s = {% if T == Complex %} Array(Float64).new(ssize, 0.0) {% else %} Array(T).new(ssize, T.new(0)) {% end %}
        rcond = {% if T == Complex %} Float64 {% else %} T {% end %}.new(cond)
        lapack(ge, lsd, rows, columns, b.columns, a, columns, x, x.columns, s, rcond, pointerof(rank))
      when .lsy?
        jpvt = Slice(Int32).new(columns)
        rcond = {% if T == Complex %} Float64 {% else %} T {% end %}.new(cond)
        lapack(ge, lsy, rows, columns, b.columns, a, columns, x, x.columns, jpvt, rcond, pointerof(rank))
        s = {% if T == Complex %} Array(Float64) {% else %} Array(T) {% end %}.new
      else
        s = {% if T == Complex %} Array(Float64) {% else %} Array(T) {% end %}.new
      end
      {x, rank, s}
    end

    def svd(*, overwrite_a = false)
      a = overwrite_a ? self : self.clone
      m = rows
      n = columns
      s = {% if T == Complex %} Slice(Float64) {% else %} Slice(T) {% end %}.new({m, n}.min)
      u = GeneralMatrix(T).new(m, m)
      vt = GeneralMatrix(T).new(n, n)
      lapack(ge, sdd, 'A'.ord, m, n, a, columns, s, u, m, vt, n)
      return {u, s, vt}
    end

    def svdvals(*, overwrite_a = false)
      a = overwrite_a ? self : self.clone
      m = rows
      n = columns
      s = {% if T == Complex %} Slice(Float64) {% else %} Slice(T) {% end %}.new({m, n}.min)
      lapack(ge, sdd, 'N'.ord, m, n, a, columns, s, nil, m, nil, n)
      s
    end

    def balance!(*, permute = true, scale = true, separate = false)
      raise ArgumentError.new("matrix must be square") unless square?
      n = @rows
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
      # idea from scypi.
      # no need to calculate if size <= 2
      if rows < 2
        q = calc_q ? Matrix(T).identity(rows) : Matrix(T).zeros(1, 1)
        return {self, q}
      end
      n = rows
      s = {% if T == Complex %} Slice(Float64) {% else %} Slice(T) {% end %}.new(n)
      lapack(ge, bal, 'S'.ord, n, self, n, out ilo, out ihi, s)
      tau = GeneralMatrix(T).new(1, n)
      lapack(ge, hrd, n, ilo, ihi, self, columns, tau)
      if calc_q
        q = clone
        {% if T == Complex %}
          lapack(un, ghr, n, ilo, ihi, q, columns, tau)
        {% else %}
          lapack(or, ghr, n, ilo, ihi, q, columns, tau)
        {% end %}
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
  end
end
