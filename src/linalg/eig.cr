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

    private def eigsh(*, need_vectors, overwrite_a = false)
      {% if T == Complex %}
      job = need_vectors ? 'V'.ord : 'N'.ord
      a = overwrite_a ? self : clone
      vals = Array(Float64).new(rows, 0.0)
      vectors = a.clone
      support = Slice(Int32).new(2*rows)
      lapack(he, evr, job, 'A'.ord, uplo, rows, a, rows,
        'N'.ord,'N'.ord,'N'.ord,'N'.ord, -1.0,
        out nfound, vals,
        vectors, columns, support)
      a.clear_flags
      vectors.clear_flags
      {vals, vectors}
      {% end %}
    end

    private def eigs_sy(*, need_vectors, overwrite_a = false)
      {% if T != Complex %}
      job = need_vectors ? 'V'.ord : 'N'.ord
      a = overwrite_a ? self : clone
      vals = Array(T).new(rows, T.new(0))
      vectors = a.clone
      support = Slice(Int32).new(2*rows)
      lapack(sy, evr, job, 'A'.ord, uplo, rows, a, rows,
        'N'.ord,'N'.ord,'N'.ord,'N'.ord, -1.0,
        out nfound, vals,
        vectors, columns, support)
      a.clear_flags
      vectors.clear_flags
      {vals, vectors}
      {% end %}
    end

    def eigs(*, need_left : Bool, need_right : Bool, overwrite_a = false)
      raise ArgumentError.new("matrix must be square") unless square?
      # TODO -  hb, sb, st
      {% if T == Complex %}
      if flags.hermitian?
        vals, vectors = eigsh(need_vectors: need_left || need_right, overwrite_a: overwrite_a)
        if need_left && need_right
          return {vals, vectors, vectors.not_nil!.clone}
        else
          return {vals, need_left ? vectors : nil, need_right ? vectors : nil}
        end
      end
      {% else %}
      if flags.symmetric?
        vals, vectors = eigs_sy(need_vectors: need_left || need_right, overwrite_a: overwrite_a)
        if need_left && need_right
          return {vals, vectors, vectors.not_nil!.clone}
        else
          return {vals, need_left ? vectors : nil, need_right ? vectors : nil}
        end
      end
      {% end %}
      a = overwrite_a ? self : clone
      eigvectorsl = need_left ? GeneralMatrix(T).new(rows, rows) : nil
      eigvectorsr = need_right ? GeneralMatrix(T).new(rows, rows) : nil
      {% if T == Complex %}
        vals = Array(T).new(rows, T.new(0,0))
        lapack(ge, ev, need_left ? 'V'.ord : 'N'.ord, need_right ? 'V'.ord : 'N'.ord, rows, a, rows,
                vals.to_unsafe.as(LibLAPACKE::DoubleComplex*),
                eigvectorsl.try &.to_unsafe, rows,
                eigvectorsr.try &.to_unsafe, rows)
        a.clear_flags
        return {vals, eigvectorsl, eigvectorsr}
      {% else %}
        reals = Array(T).new(rows, T.new(0))
        imags = Array(T).new(rows, T.new(0))
        lapack(ge, ev, need_left ? 'V'.ord : 'N'.ord, need_right ? 'V'.ord : 'N'.ord, rows, a, rows,
                reals, imags,
                eigvectorsl.try &.to_unsafe, rows,
                eigvectorsr.try &.to_unsafe, rows)
        a.clear_flags
        if imags.all? &.==(0)
          vals = reals
        else
          vals = Array(Complex).new(rows){|i| Complex.new(reals[i], imags[i])}
        end
        return {vals, eigvectorsl, eigvectorsr}
      {% end %}
    end

    private def eigs_gen_sy(*, b : Matrix(T), need_vectors : Bool, overwrite_a = false, overwrite_b = false)
      {% if T != Complex %}
      job = need_vectors ? 'V'.ord : 'N'.ord
      a = overwrite_a ? self : clone
      bb = overwrite_b ? b : b.clone
      vals = Array(T).new(rows, T.new(0))
      lapack(sy, gvd, 1, job, uplo,
        rows, a, rows,
        bb, rows,
        vals)
      a.clear_flags
      bb.clear_flags
      {vals, a}
      {% end %}
    end

    private def eigs_gen_he(*, b : Matrix(T), need_vectors : Bool, overwrite_a = false, overwrite_b = false)
      {% if T == Complex %}
      job = need_vectors ? 'V'.ord : 'N'.ord
      a = overwrite_a ? self : clone
      bb = overwrite_b ? b : b.clone
      vals = Array(Float64).new(rows, 0.0)
      lapack(he, gvd, 1, job, uplo,
        rows, a, rows,
        bb, rows,
        vals)
      a.clear_flags
      bb.clear_flags
      {vals, a}
      {% end %}
    end

    # generalized eigenvalues problem
    def eigs(*, b : Matrix(T), need_left : Bool, need_right : Bool, overwrite_a = false, overwrite_b = false)
      raise ArgumentError.new("a matrix must be square") unless square?
      raise ArgumentError.new("b matrix must be square") unless b.square?
      {% if T == Complex %}
      if flags.hermitian? && b.flags.positive_definite?
        vals, vectors = eigs_gen_he(b: b, need_vectors: need_left || need_right, overwrite_a: overwrite_a, overwrite_b: overwrite_b)
        beta = Array(Float64).new(rows, 1.0)
        if need_left && need_right
          return {vals, beta, vectors, vectors.not_nil!.clone}
        else
          return {vals, beta, need_left ? vectors : nil, need_right ? vectors : nil}
        end
      end
      {% else %}
      if flags.symmetric? && b.flags.positive_definite?
        vals, vectors = eigs_gen_sy(b: b, need_vectors: need_left || need_right, overwrite_a: overwrite_a, overwrite_b: overwrite_b)
        beta = Array(T).new(rows, T.new(1))
        if need_left && need_right
          return {vals, beta, vectors, vectors.not_nil!.clone}
        else
          return {vals, beta, need_left ? vectors : nil, need_right ? vectors : nil}
        end
      end
      {% end %}
      a = overwrite_a ? self : clone
      bb = overwrite_b ? b : b.clone
      eigvectorsl = need_left ? GeneralMatrix(T).new(rows, rows) : nil
      eigvectorsr = need_right ? GeneralMatrix(T).new(rows, rows) : nil
      {% if T == Complex %}
        alpha = Array(T).new(rows, T.new(0,0))
        beta = Array(T).new(rows, T.new(0,0))
        lapack(gg, ev, need_left ? 'V'.ord : 'N'.ord, need_right ? 'V'.ord : 'N'.ord, rows, a, rows,
                bb, b.columns,
                alpha.to_unsafe.as(LibLAPACKE::DoubleComplex*),
                beta.to_unsafe.as(LibLAPACKE::DoubleComplex*),
                eigvectorsl.try &.to_unsafe, rows,
                eigvectorsr.try &.to_unsafe, rows)
        a.clear_flags
        bb.clear_flags
        return {alpha, beta, eigvectorsl, eigvectorsr}
      {% else %}
        alpha_reals = Array(T).new(rows, T.new(0))
        alpha_imags = Array(T).new(rows, T.new(0))
        beta = Array(T).new(rows, T.new(0))
        lapack(gg, ev, need_left ? 'V'.ord : 'N'.ord, need_right ? 'V'.ord : 'N'.ord, rows, a, rows,
                overwrite_b ? b : b.clone, b.columns,
                alpha_reals, alpha_imags, beta,
                eigvectorsl.try &.to_unsafe, rows,
                eigvectorsr.try &.to_unsafe, rows)
        a.clear_flags
        bb.clear_flags
        if alpha_imags.all? &.==(0)
          alpha = alpha_reals
        else
          alpha = Array(Complex).new(rows){|i| Complex.new(alpha_reals[i], alpha_imags[i])}
        end
        return {alpha, beta, eigvectorsl, eigvectorsr}
      {% end %}
    end
  end
end
