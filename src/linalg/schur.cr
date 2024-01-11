require "../matrix/*"
require "./lapack_helper"

# TODO - inline docs

module LA
  def self.qz(a, b, *, overwrite_a = false, overwrite_b = false)
    a.qz(b, overwrite_a: overwrite_a, overwrite_b: overwrite_b)
  end

  class Matrix(T)
    def schur
      to_general.schur(overwrite_a: true)
    end

    def qz(b : self)
      to_general.qz(b, overwrite_a: true, overwrite_b: true)
    end
  end

  class GeneralMatrix(T) < Matrix(T)
    def schur(*, overwrite_a = false)
      raise ArgumentError.new("matrix must be square") unless square?
      a = overwrite_a ? self : clone
      {% if T == Complex %}
        w = Slice(T).new(nrows, T.new(0))
        z = GeneralMatrix(T).new(*size)
        lapack(gees, 'V'.ord.to_u8, 'N'.ord.to_u8, nil, nrows, a, ncolumns, sdim, w.to_unsafe.as(LibCBLAS::ComplexDouble*), z, ncolumns, worksize: [nrows])
        a.clear_flags
        a.assume! MatrixFlags::UpperTriangular
        z.assume! MatrixFlags::Orthogonal
        return {a, z}
      {% else %}
        wr = Slice(T).new(nrows)
        wi = Slice(T).new(nrows)
        z = GeneralMatrix(T).new(*size)
        lapack(gees, 'V'.ord.to_u8, 'N'.ord.to_u8, nil, nrows, a, ncolumns, sdim, wr, wi, z, ncolumns)
        a.clear_flags
        a.detect MatrixFlags::UpperTriangular
        z.assume! MatrixFlags::Orthogonal
        return {a, z}
      {% end %}
    end

    def qz(b : self, overwrite_a = false, overwrite_b = false)
      raise ArgumentError.new("matrix must be square") unless square?
      a = overwrite_a ? self : clone
      bb = overwrite_b ? b : b.clone
      vsl = GeneralMatrix(T).new(*size)
      vsr = GeneralMatrix(T).new(*size)
      {% if T == Complex %}
        alpha = Slice(T).new(nrows, T.new(0))
        beta = Slice(T).new(nrows, T.new(0))
        lapack(gges, 'V'.ord.to_u8, 'V'.ord.to_u8, 'N'.ord.to_u8, nil,
          nrows, a, ncolumns, bb, bb.ncolumns,
          sdim,
          alpha.to_unsafe.as(LibCBLAS::ComplexDouble*),
          beta.to_unsafe.as(LibCBLAS::ComplexDouble*),
          vsl, ncolumns, vsr, ncolumns, worksize: [8*nrows])
        a.clear_flags
        bb.clear_flags
        bb.assume! MatrixFlags::UpperTriangular
      {% else %}
        alphar = Slice(T).new(nrows)
        alphai = Slice(T).new(nrows)
        beta = Slice(T).new(nrows)
        lapack(gges, 'V'.ord.to_u8, 'V'.ord.to_u8, 'N'.ord.to_u8, nil, nrows, a, ncolumns, bb, ncolumns, sdim,
          alphar, alphai, beta, vsl, ncolumns, vsr, ncolumns)
        a.clear_flags
        bb.clear_flags
        bb.detect MatrixFlags::UpperTriangular
      {% end %}
      vsl.assume! MatrixFlags::Orthogonal
      vsr.assume! MatrixFlags::Orthogonal
      return {a, bb, vsl, vsr}
    end
  end
end
