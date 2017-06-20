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
       LibLAPACKE.{{typ}}{{storage}}{{name}}(LibLAPACKE::ROW_MAJOR, {{*args}})
    end

    def inv!
      raise "can't invert nonsquare matrix" unless @rows == @columns
      n = @rows
      ipiv = Slice(Int32).new(n)
      lapack(trf, n, n, result, n, ipiv)
      lapack(tri, n, result, n, ipiv)
      self
    end

    def inv
      clone.inv!
    end
  end
end
