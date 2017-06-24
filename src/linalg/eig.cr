require "../matrix/*"
require "./libLAPACKE"

module Linalg
  module Matrix(T)
    def eigvals(*, overwrite_a = false)
      vals, aleft, aright = eigs(overwrite_a: overwrite_a, need_left: false, need_right: false)
      vals
    end

    def eigs(*, left = false, overwrite_a = false)
      vals, aleft, aright = eigs(overwrite_a: overwrite_a, need_left: left, need_right: !left)
      v = left ? aleft : aright
      {vals, v.not_nil!}
    end

    def eigs(*, need_left : Bool, need_right : Bool, overwrite_a = false)
      raise ArgumentError.new("matrix must be square") unless square?
      a = overwrite_a ? self : self.clone
      eigvectorsl = need_left ? GeneralMatrix(T).new(rows, rows) : nil
      eigvectorsr = need_right ? GeneralMatrix(T).new(rows, rows) : nil
      {% if T == Complex %}
        vals = Array(T).new(rows, T.new(0,0))
        lapack(ge, ev, need_left ? 'V'.ord : 'N'.ord, need_right ? 'V'.ord : 'N'.ord, rows, a, rows,
                vals.to_unsafe.as(LibLAPACKE::DoubleComplex*),
                eigvectorsl.try &.to_unsafe, rows,
                eigvectorsr.try &.to_unsafe, rows)
        return {vals, eigvectorsl, eigvectorsr}
      {% else %}
        reals = Array(T).new(rows, T.new(0))
        imags = Array(T).new(rows, T.new(0))
        lapack(ge, ev, need_left ? 'V'.ord : 'N'.ord, need_right ? 'V'.ord : 'N'.ord, rows, a, rows,
                reals, imags,
                eigvectorsl.try &.to_unsafe, rows,
                eigvectorsr.try &.to_unsafe, rows)
        if imags.all? &.==(0)
          vals = reals
        else
          vals = Array(Complex).new(rows){|i| Complex.new(reals[i], imags[i])}
        end
        return {vals, eigvectorsl, eigvectorsr}
      {% end %}
    end
  end
end
