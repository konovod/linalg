require "../matrix/*"
require "./lapack_helper"

module LA
  def self.eigs(a, b, *, need_left : Bool, need_right : Bool, overwrite_a = false, overwrite_b = false)
    a.eigs(b: b, need_left: need_left, need_right: need_right, overwrite_a: overwrite_a, overwrite_b: overwrite_b)
  end

  abstract class Matrix(T)
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
        job = need_vectors ? 'V'.ord.to_u8 : 'N'.ord.to_u8
        a = overwrite_a ? self : clone
        vals = Array(Float64).new(nrows, 0.0)
        vectors = a.clone
        support = Slice(Int32).new(2*nrows)
        lapack(heevr, job, 'A'.ord.to_u8, uplo, nrows, a, nrows,
          0.0, 0.0, 0, 0, -1.0,
          nfound, vals,
          vectors, ncolumns, support)
        a.clear_flags
        vectors.clear_flags
        {vals, vectors}
      {% end %}
    end

    private def eigs_sy(*, need_vectors, overwrite_a = false)
      {% if T != Complex %}
        job = need_vectors ? 'V'.ord.to_u8 : 'N'.ord.to_u8
        a = overwrite_a ? self : clone
        vals = Array(T).new(nrows, T.new(0))
        vectors = a.clone
        support = Slice(Int32).new(2*nrows)
        lapack(syevr, job, 'A'.ord.to_u8, uplo, nrows, a, nrows,
          T.new(0.0), T.new(0.0), 0, 0, T.new(-1.0),
          nfound, vals,
          vectors, ncolumns, support)

        a.clear_flags
        vectors.clear_flags
        {vals, vectors}
      {% end %}
    end

    def eigs(*, need_left : Bool, need_right : Bool, overwrite_a = false)
      raise ArgumentError.new("matrix must be square") unless square?
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
      {% if flag?(:darwin) && T == Float32 %}
        raise "Eigenvectors for single precision on Darwin are not supported for now"
      {% end %}
      a = overwrite_a ? self : clone
      eigvectorsl = need_left ? GeneralMatrix(T).new(nrows, nrows) : nil
      eigvectorsr = need_right ? GeneralMatrix(T).new(nrows, nrows) : nil
      {% if T == Complex %}
        vals = Array(T).new(nrows, T.new(0, 0))
        lapack(geev, need_left ? 'V'.ord.to_u8 : 'N'.ord.to_u8, need_right ? 'V'.ord.to_u8 : 'N'.ord.to_u8, nrows, a, nrows,
          vals.to_unsafe.as(LibCBLAS::ComplexDouble*),
          eigvectorsl ? eigvectorsl.to_unsafe : Pointer(LibCBLAS::ComplexDouble).null, nrows,
          eigvectorsr ? eigvectorsr.to_unsafe : Pointer(LibCBLAS::ComplexDouble).null, nrows, worksize: [2*nrows])
        a.clear_flags
        return {vals, eigvectorsl, eigvectorsr}
      {% else %}
        reals = Array(T).new(nrows, T.new(0))
        imags = Array(T).new(nrows, T.new(0))
        lapack(geev, (need_left ? 'V' : 'N').ord.to_u8, (need_right ? 'V' : 'N').ord.to_u8, nrows, a, nrows,
          reals, imags,
          eigvectorsl ? eigvectorsl.to_unsafe : Pointer(T).null, nrows,
          eigvectorsr ? eigvectorsr.to_unsafe : Pointer(T).null, nrows)
        a.clear_flags
        if imags.all? &.==(0)
          vals = reals
        else
          vals = Array(Complex).new(nrows) { |i| Complex.new(reals[i], imags[i]) }
        end
        return {vals, eigvectorsl, eigvectorsr}
      {% end %}
    end

    private def eigs_gen_sy(*, b : Matrix(T), need_vectors : Bool, overwrite_a = false, overwrite_b = false)
      {% if T != Complex %}
        job = (need_vectors ? 'V' : 'N').ord.to_u8
        a = overwrite_a ? self : clone
        bb = overwrite_b ? b : b.clone
        vals = Array(T).new(nrows, T.new(0))
        lapack(sygvd, 1, job, uplo,
          nrows, a, nrows,
          bb, nrows,
          vals)
        a.clear_flags
        bb.clear_flags
        {vals, a}
      {% end %}
    end

    private def eigs_gen_he(*, b : Matrix(T), need_vectors : Bool, overwrite_a = false, overwrite_b = false)
      {% if T == Complex %}
        job = (need_vectors ? 'V' : 'N').ord.to_u8
        a = overwrite_a ? self : clone
        bb = overwrite_b ? b : b.clone
        vals = Array(Float64).new(nrows, 0.0)
        lapack(hegvd, 1, job, uplo,
          nrows, a, nrows,
          bb, nrows,
          vals)

        a.clear_flags
        bb.clear_flags
        {vals, a}
      {% end %}
    end

    # generalized eigenvalues problem
    def eigs(*, b : Matrix(T), need_left : Bool, need_right : Bool, overwrite_a = false, overwrite_b = false)
      raise ArgumentError.new("a matrix must be square") unless square?
      raise ArgumentError.new("b matrix must have same size as a") unless b.size == self.size
      {% if T == Complex %}
        if flags.hermitian? && b.flags.positive_definite?
          vals, vectors = eigs_gen_he(b: b, need_vectors: need_left || need_right, overwrite_a: overwrite_a, overwrite_b: overwrite_b)
          beta = Array(Float64).new(nrows, 1.0)
          if need_left && need_right
            return {vals, beta, vectors, vectors.not_nil!.clone}
          else
            return {vals, beta, need_left ? vectors : nil, need_right ? vectors : nil}
          end
        end
      {% else %}
        if flags.symmetric? && b.flags.positive_definite?
          vals, vectors = eigs_gen_sy(b: b, need_vectors: need_left || need_right, overwrite_a: overwrite_a, overwrite_b: overwrite_b)
          beta = Array(T).new(nrows, T.new(1))
          if need_left && need_right
            return {vals, beta, vectors, vectors.not_nil!.clone}
          else
            return {vals, beta, need_left ? vectors : nil, need_right ? vectors : nil}
          end
        end
      {% end %}
      a = overwrite_a ? self : clone
      bb = overwrite_b ? b : b.clone
      eigvectorsl = need_left ? GeneralMatrix(T).new(nrows, nrows) : nil
      eigvectorsr = need_right ? GeneralMatrix(T).new(nrows, nrows) : nil
      {% if T == Complex %}
        alpha = Array(T).new(nrows, T.new(0, 0))
        beta = Array(T).new(nrows, T.new(0, 0))
        lapack(ggev, (need_left ? 'V' : 'N').ord.to_u8, (need_right ? 'V' : 'N').ord.to_u8, nrows, a, nrows,
          bb, b.nrows,
          alpha.to_unsafe.as(LibCBLAS::ComplexDouble*),
          beta.to_unsafe.as(LibCBLAS::ComplexDouble*),
          eigvectorsl ? eigvectorsl.to_unsafe : Pointer(LibCBLAS::ComplexDouble).null, nrows,
          eigvectorsr ? eigvectorsr.to_unsafe : Pointer(LibCBLAS::ComplexDouble).null, nrows, worksize: [8*nrows])
        a.clear_flags
        bb.clear_flags
        return {alpha, beta, eigvectorsl, eigvectorsr}
      {% else %}
        alpha_reals = Array(T).new(nrows, T.new(0))
        alpha_imags = Array(T).new(nrows, T.new(0))
        beta = Array(T).new(nrows, T.new(0))
        lapack(ggev, (need_left ? 'V' : 'N').ord.to_u8, (need_right ? 'V' : 'N').ord.to_u8, nrows, a, nrows,
          overwrite_b ? b : b.clone, b.nrows,
          alpha_reals, alpha_imags, beta,
          eigvectorsl ? eigvectorsl.to_unsafe : Pointer(T).null, nrows,
          eigvectorsr ? eigvectorsr.to_unsafe : Pointer(T).null, nrows)
        a.clear_flags
        bb.clear_flags
        if alpha_imags.all? &.==(0)
          alpha = alpha_reals
        else
          alpha = Array(Complex).new(nrows) { |i| Complex.new(alpha_reals[i], alpha_imags[i]) }
        end
        return {alpha, beta, eigvectorsl, eigvectorsr}
      {% end %}
    end
  end
end
