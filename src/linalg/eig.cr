require "../matrix/*"
require "./lapack_helper"

module LA
  # See `GeneralMatrix#eigss(*, b : GeneralMatrix(T), need_left : Bool, need_right : Bool, overwrite_a = false, overwrite_b = false)`.
  def self.eigs(a, b, *, need_left : Bool, need_right : Bool, overwrite_a = false, overwrite_b = false)
    a.eigs(b: b, need_left: need_left, need_right: need_right, overwrite_a: overwrite_a, overwrite_b: overwrite_b)
  end

  class Matrix(T)
    # See `GeneralMatrix#eigs(*, need_left : Bool, need_right : Bool, overwrite_a = false)`.
    # This version returns only eigenvalues
    def eigvals
      to_general.eigvals(overwrite_a: true)
    end

    # See `GeneralMatrix#eigs(*, need_left : Bool, need_right : Bool, overwrite_a = false)`.
    # This version returns eigenvalues and either left or right eigenvectors (depending on `left`)
    def eigs(*, left = false)
      to_general.eigs(overwrite_a: true, left: left)
    end

    # See `GeneralMatrix#eigs(*, need_left : Bool, need_right : Bool, overwrite_a = false)`
    def eigs(*, need_left : Bool, need_right : Bool)
      to_general.eigs(overwrite_a: true, need_left: need_left, need_right: need_right)
    end

    # See `GeneralMatrix#eigss(*, b : GeneralMatrix(T), need_left : Bool, need_right : Bool, overwrite_a = false, overwrite_b = false)`.
    def eigs(*, b : self, need_left = false, need_right = false)
      to_general.eigs(b: b.to_general, overwrite_a: true, overwrite_b: true, need_left: need_left, need_right: need_right)
    end
  end

  class GeneralMatrix(T) < Matrix(T)
    # - See `GeneralMatrix#eigs(*, need_left : Bool, need_right : Bool, overwrite_a = false)`.
    #  This version returns only eigenvalues
    def eigvals(*, overwrite_a = false)
      vals, aleft, aright = eigs(overwrite_a: overwrite_a, need_left: false, need_right: false)
      vals
    end

    # - See `GeneralMatrix#eigs(*, need_left : Bool, need_right : Bool, overwrite_a = false)`.
    #  This version returns eigenvalues and either left or right eigenvectors (depending on `left`)
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
        vectors = need_vectors ? GeneralMatrix(T).new(*a.size) : GeneralMatrix(T).new(0, 0)
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
        vectors = need_vectors ? GeneralMatrix(T).new(*a.size) : GeneralMatrix(T).new(0, 0)
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

    # Computes eigenvalues and optionally left and right eigenvectors.
    #
    # Arguments:
    # - need_left: Whether to compute left eigenvectors.
    # - need_right: Whether to compute right eigenvectors.
    # - overwrite_a: If true, allows overwriting matrix `self`.
    #
    # Returns:
    #   A tuple `{values, left_vectors, right_vectors}` where:
    #   - `values` is an `Array(T)` of eigenvalues (or `Array(Complex)` if complex eigenvalues arise in real matrices).
    #   - `left_vectors` is a `GeneralMatrix(T)?` containing left eigenvectors if `need_left == true`, else `nil`.
    #   - `right_vectors` is a `GeneralMatrix(T)?` containing right eigenvectors if `need_right == true`, else `nil`.
    #
    # For Hermitian (complex) or symmetric (real) matrices, eigenvalues are real, and eigenvectors are orthonormal.
    #
    # Exceptions:
    # - `ArgumentError`: if the matrix is not square.
    #
    # Called LAPACK routines:
    # - `heevr` — for complex Hermitian matrices.
    # - `syevr` — for real symmetric matrices.
    # - `geev`  — for general complex or real non-symmetric/Hermitian matrices.
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

    # Computes generalized eigenvalues and optionally left and right eigenvectors.
    #
    # Arguments:
    # - b: Right-hand matrix in the generalized eigenvalue problem.
    # - need_left: Whether to compute left eigenvectors.
    # - need_right: Whether to compute right eigenvectors.
    # - overwrite_a: If true, allows overwriting matrix `self`.
    # - overwrite_b: If true, allows overwriting matrix `b`.
    #
    # Returns:
    #   A tuple `{alpha, beta, left_vectors, right_vectors}` where:
    #   - `alpha` is an `Array(T)` of generalized eigenvalue numerators (or `Array(Complex)` if complex values arise).
    #   - `beta` is an `Array(T)` of generalized eigenvalue denominators (`λ = alpha[i] / beta[i]`).
    #   - `left_vectors` is a `GeneralMatrix(T)?` containing left eigenvectors if `need_left == true`, else `nil`.
    #   - `right_vectors` is a `GeneralMatrix(T)?` containing right eigenvectors if `need_right == true`, else `nil`.
    #
    # For Hermitian/symmetric `A` and positive definite `B`, the problem is solved using a specialized algorithm ensuring real eigenvalues.
    #
    # Exceptions:
    # - `ArgumentError`: if `self` is not square or if `b` does not have the same size as `self`.
    #
    # Called LAPACK routines:
    # - `hegvd` — for complex Hermitian-definite generalized problems.
    # - `sygvd` — for real symmetric-definite generalized problems.
    # - `ggev`  — for general complex or real non-symmetric/Hermitian problems.
    def eigs(*, b : GeneralMatrix(T), need_left : Bool, need_right : Bool, overwrite_a = false, overwrite_b = false)
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
