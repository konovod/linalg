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

  def self.solve(a, b)
    a.solve(b)
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

    def inv!
      raise ArgumentError.new("can't invert nonsquare matrix") unless square?
      raise ArgumentError.new("can't invert inplace virtual matrix") if flags.virtual?
      n = @rows
      ipiv = Slice(Int32).new(n)
      lapack(ge, trf, n, n, self, n, ipiv)
      lapack(ge, tri, n, self, n, ipiv)
      self
    end

    def inv
      clone.inv!
    end

    def solve(b : self)
      raise ArgumentError.new("number of rows in a and b must match") unless rows == b.rows
      raise ArgumentError.new("a must be square") unless square?
      x = b.clone
      n = rows
      ipiv = Slice(Int32).new(n)
      lapack(ge, sv, n, b.columns, self.clone, n, ipiv, x, b.columns)
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

    def eigvals(*, overwrite_a = false)
      raise ArgumentError.new("matrix must be square") unless square?
      a = overwrite_a ? self : self.clone
      {% if T == Complex %}
        vals = Array(T).new(rows, T.new(0,0))
        lapack(ge, ev, 'N'.ord, 'N'.ord, rows, a, rows, vals.to_unsafe.as(LibLAPACKE::DoubleComplex*), nil, rows, nil, rows)
        return vals
      {% else %}
        reals = Array(T).new(rows, T.new(0))
        imags = Array(T).new(rows, T.new(0))
        lapack(ge, ev, 'N'.ord, 'N'.ord, rows, a, rows, reals, imags, nil, rows, nil, rows)
        if imags.all? &.==(0)
          return reals
        else
          return Array(Complex).new(rows){|i| Complex.new(reals[i], imags[i])}
        end
      {% end %}
    end

    def eigs(*, left : Bool = false, overwrite_a = false)
      raise ArgumentError.new("matrix must be square") unless square?
      a = overwrite_a ? self : self.clone
      eigvectors = GeneralMatrix(T).new(rows, rows)
      {% if T == Complex %}
        vals = Array(T).new(rows, T.new(0,0))
        lapack(ge, ev, left ? 'V'.ord : 'N'.ord, left ? 'N'.ord : 'V'.ord, rows, a, rows,
                vals.to_unsafe.as(LibLAPACKE::DoubleComplex*),
                left ? eigvectors.to_unsafe.as(LibLAPACKE::DoubleComplex*) : nil, rows,
                left ? nil : eigvectors.to_unsafe.as(LibLAPACKE::DoubleComplex*), rows)
        return {vals, eigvectors}
      {% else %}
        reals = Array(T).new(rows, T.new(0))
        imags = Array(T).new(rows, T.new(0))
        lapack(ge, ev, left ? 'V'.ord : 'N'.ord, left ? 'N'.ord : 'V'.ord, rows, a, rows, reals, imags, left ? eigvectors.to_unsafe : nil, rows, left ? nil : eigvectors.to_unsafe, rows)
        if imags.all? &.==(0)
          values = reals
        else
          values = Array(Complex).new(rows){|i| Complex.new(reals[i], imags[i])}
        end
        {values, eigvectors}
      {% end %}
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
  end
end
