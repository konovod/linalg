require "./matrix"
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

  def self.inv(matrix)
    matrix.inv
  end

  def self.solve(a, b)
    a.solve(b)
  end

  def self.lstsq(a, b, method : LSMethod = LSMethod::Auto)
    a.solvels(b, method)
  end

  def self.svd(matrix)
    matrix.svd
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

    def det(*, overwrite_a = false)
      raise ArgumentError.new("matrix must be square") unless square?
      lru = overwrite_a ? self : self.clone
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

    def eigvals(*, overwrite_a = false)
      raise ArgumentError.new("matrix must be square") unless square?
      a = overwrite_a ? self : self.clone
      {% if T == Complex %}
        vals = Array(T).new(rows, T.new(0,0))
        lapack(ev, 'N'.ord, 'N'.ord, rows, a, rows, vals.to_unsafe.as(LibLAPACKE::DoubleComplex*), nil, rows, nil, rows)
        return vals
      {% else %}
        reals = Array(T).new(rows, T.new(0))
        imags = Array(T).new(rows, T.new(0))
        lapack(ev, 'N'.ord, 'N'.ord, rows, a, rows, reals, imags, nil, rows, nil, rows)
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
      eigvectors = Matrix(T).new(rows, rows)
      {% if T == Complex %}
        vals = Array(T).new(rows, T.new(0,0))
        lapack(ev, left ? 'V'.ord : 'N'.ord, left ? 'N'.ord : 'V'.ord, rows, a, rows,
                vals.to_unsafe.as(LibLAPACKE::DoubleComplex*),
                left ? eigvectors.to_unsafe.as(LibLAPACKE::DoubleComplex*) : nil, rows,
                left ? nil : eigvectors.to_unsafe.as(LibLAPACKE::DoubleComplex*), rows)
        return {vals, eigvectors}
      {% else %}
        reals = Array(T).new(rows, T.new(0))
        imags = Array(T).new(rows, T.new(0))
        lapack(ev, left ? 'V'.ord : 'N'.ord, left ? 'N'.ord : 'V'.ord, rows, a, rows, reals, imags, left ? eigvectors.to_unsafe : nil, rows, left ? nil : eigvectors.to_unsafe, rows)
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
      u = Matrix(T).new(m, m)
      vt = Matrix(T).new(n, n)
      lapack(sdd, 'A'.ord, m, n, a, columns, s, u, m, vt, n)
      return {u, s, vt}
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
      s = Matrix(T).new(1, n)
      lapack(bal, job.ord, n, self, n, out ilo, out ihi, s)
      separate ? s : Matrix(T).diag(s.raw)
    end

    def balance(*, permute = true, scale = true, separate = false)
      a = clone
      s = a.balance!(permute: permute, scale: scale, separate: separate)
      {a, s}
    end

    def lu(*, overwrite_a = false)
      a = overwrite_a ? self : self.clone
      m = rows
      n = columns
      k = {rows, columns}.min
      ipiv = Slice(Int32).new(m)
      lapack(trf, rows, columns, a, columns, ipiv)
      # TODO - better solution?
      # apply all transformation of piv to "own" piv
      piv = (1..m).to_a
      k.times do |i|
        tmp = piv[i]
        piv[i] = piv[ipiv[i] - 1]
        piv[ipiv[i] - 1] = tmp
      end
      p = Matrix(T).new(m, m)
      m.times do |i|
        p[piv[i] - 1, i] = T.new(1)
      end
      l = Matrix(T).new(m, k) do |i, j|
        case i <=> j
        when 0
          T.new(1)
        when 1
          a[i, j]
        else
          T.new(0)
        end
      end
      u = Matrix(T).new(k, n) do |i, j|
        if i <= j
          a[i, j]
        else
          T.new(0)
        end
      end
      {p, l, u}
    end

    def lu_factor!
      raise ArgumentError.new("matrix must be square") unless square?
      ipiv = Slice(Int32).new(rows)
      lapack(trf, rows, columns, self, columns, ipiv)
      LUMatrix(T).new(self, ipiv)
    end

    def lu_factor
      clone.lu_factor!
    end
  end

  enum LUTranspose
    None          = 0
    Transpose
    ConjTranspose
  end

  struct LUMatrix(T)
    @a : Matrix(T)
    @ipiv : Slice(Int32)

    # TODO - more macro magic?
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

    def initialize(@a, @ipiv)
    end

    def size
      @a.rows
    end

    # TODO - equilibration?

    def solve(b, transpose = LUTranspose::None, *, overwrite_b = false)
      raise ArgumentError.new("number of rows in a and b must match") unless @a.rows == b.rows
      trans = case transpose
              when .none?           then 'N'
              when .transpose?      then 'T'
              when .conj_transpose? then 'C'
              else                       'N'
              end.ord
      x = overwrite_b ? b : b.clone
      lapack(trs, trans, size, b.columns, @a, size, @ipiv, x, x.columns)
      x
    end
  end
end
