require "./matrix"
require "./libLAPACKE"

module LAPACK
  # matrix inversion using xxxtrf / xxxtri

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
       raise "lapack returned #{info}" if info != 0
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
  end

  def inv(matrix)
    matrix.inv
  end

  def solve(a, b)
    a.solve(b)
  end
end
