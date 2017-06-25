require "../matrix/*"
require "./libLAPACKE"

module Linalg
  def self.qz(a, b, *, overwrite_a = false, overwrite_b = false)
    a.qz(b, overwrite_a: overwrite_a, overwrite_b: overwrite_b)
  end

  module Matrix(T)
    def schur(*, overwrite_a = false)
      raise ArgumentError.new("matrix must be square") unless square?
      a = overwrite_a ? self : clone
      {% if T == Complex %}
        w = Slice(T).new(rows, T.new(0))
        z = GeneralMatrix(T).new(*size)
        lapack(ge, es, 'V'.ord, 'N'.ord, nil, rows, a, columns, out sdim, w.to_unsafe.as(LibLAPACKE::DoubleComplex*), z, columns)
        a.clear_flags
        a.assume! MatrixFlags::Triangular
        z.assume! MatrixFlags::Orthogonal
        return {a, z}
      {% else %}
      wr = Slice(T).new(rows)
      wi = Slice(T).new(rows)
      z = GeneralMatrix(T).new(*size)
      lapack(ge, es, 'V'.ord, 'N'.ord, nil, rows, a, columns, out sdim, wr, wi, z, columns)
      a.clear_flags
      z.assume! MatrixFlags::Orthogonal
      return {a, z}
      {% end %}
    end

    def qz(b, overwrite_a = false, overwrite_b = false)
      raise ArgumentError.new("matrix must be square") unless square?
      a = overwrite_a ? self : clone
      bb = overwrite_b ? b : b.clone
      vsl = GeneralMatrix(T).new(*size)
      vsr = GeneralMatrix(T).new(*size)
      {% if T == Complex %}
        alpha = Slice(T).new(rows, T.new(0))
        beta = Slice(T).new(rows, T.new(0))
        lapack(gg, es, 'V'.ord,'V'.ord, 'N'.ord, nil,
            rows, a, columns, bb, bb.columns,
            out sdim,
            alpha.to_unsafe.as(LibLAPACKE::DoubleComplex*),
            beta.to_unsafe.as(LibLAPACKE::DoubleComplex*),
            vsl, columns, vsr, columns)
        a.clear_flags
        bb.clear_flags
        bb.assume! MatrixFlags::Triangular
      {% else %}
        alphar = Slice(T).new(rows)
        alphai = Slice(T).new(rows)
        beta = Slice(T).new(rows)
        lapack(gg, es, 'V'.ord, 'V'.ord, 'N'.ord, nil, rows, a, columns, bb, columns, out sdim,
            alphar, alphai, beta, vsl, columns,vsr, columns)
        a.clear_flags
        bb.clear_flags
      {% end %}
      vsl.assume! MatrixFlags::Orthogonal
      vsr.assume! MatrixFlags::Orthogonal
      return {a, bb, vsl, vsr}
    end
  end
end
