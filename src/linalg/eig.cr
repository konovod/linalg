require "../matrix/*"
require "./libLAPACKE"

module Linalg
  module Matrix(T)
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
  end
end
