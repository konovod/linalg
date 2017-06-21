require "./matrix"
require "./libLAPACKE"

module Linalg
  # matrix inversion using xxxtrf / xxxtri

  enum LSMethod
    Auto       = 0
    QR
    Orthogonal
    SVD
    LS         = QR
    LSY        = Orthogonal
    LSD        = SVD
  end

  class Matrix(T)
    macro lapack(name, *args)
      {% storage = :ge.id %}
      {% if T == Float32
           typ = :s.id
         elsif T == Float64
           typ = :d.id
         elsif T == Complex
           typ = :z.id
         end %}
       info = LibLAPACKE.{{typ}}{{storage}}{{name}}(LibLAPACKE::ROW_MAJOR, {{*args}})
       raise "LAPACKE.{{typ}}{{storage}}{{name}} returned #{info}" if info != 0
    end

    def inv!
      raise ArgumentError.new("can't invert nonsquare matrix") unless square?
      n = @rows
      ipiv = Slice(Int32).new(n)
      lapack(trf, n, n, self, n, ipiv)
      lapack(tri, n, self, n, ipiv)
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
      lapack(sv, n, b.columns, self.clone, n, ipiv, x, b.columns)
      x
    end

    def det
      raise ArgumentError.new("matrix must be square") unless square?
      lru = self.clone
      ipiv = Slice(Int32).new(rows)
      lapack(trf, rows, rows, lru, rows, ipiv)
      (0...rows).reduce(1) { |det, i| det*lru[i, i] }
    end

    def solvels(b : self, method : LSMethod = LSMethod::Auto)
      raise ArgumentError.new("number of rows in a and b must match") unless rows == b.rows
      if method.auto?
        method = LSMethod::QR
      end
      a = self.clone
      if columns > rows
        # make room for residuals
        x = Matrix(T).new(columns, b.columns) { |r, c| r < rows ? b[r, c] : T.new(0) }
      else
        x = b.clone
      end
      case method
      when .ls?
        lapack(ls, 'N'.ord, rows, columns, b.columns, a, columns, x, x.columns)
      when .lsd?
        s = {% if T == Complex %} Slice(Float64) {% else %} Slice(T) {% end %}.new({rows, columns}.min)
        rcond = {% if T == Complex %} Float64 {% else %} T {% end %}.new(-1)
        lapack(lsd, rows, columns, b.columns, a, columns, x, x.columns, s, rcond, out rank1)
      when .lsy?
        jpvt = Slice(Int32).new(columns)
        rcond = {% if T == Complex %} Float64 {% else %} T {% end %}.new(-1)
        lapack(lsy, rows, columns, b.columns, a, columns, x, x.columns, jpvt, rcond, out rank2)
      end
      x
    end

    def eigvals
      raise ArgumentError.new("matrix must be square") unless square?
      a = self.clone
      {% if T == Complex %}
      {% else %}
        reals = Array(T).new(rows, T.new(0))
        imags = Array(T).new(rows, T.new(0))
        lapack(ev, 'N'.ord, 'N'.ord, rows, a, rows, reals, imags, nil, rows, nil, rows)
        if imags.all? &.==(0)
          reals
        else
          Array(Complex).new(rows){|i| Complex.new(reals[i], imags[i])}
        end
      {% end %}
    end
  end

  def self.inv(matrix)
    matrix.inv
  end

  def self.solve(a, b)
    a.solve(b)
  end

  def self.lstsq(a, b, method : LSMethod = LSMethod::Auto)
    a.solvels(b, method)
  end
end
