require "../matrix/*"
require "./lapack_helper"

module LA
  # See `GeneralMatrix#cho_solve`
  def self.cho_solve(a, b, *, overwrite_b = false)
    a.cho_solve(b, overwrite_b: overwrite_b)
  end

  abstract class Matrix(T)
    # Computes the Cholesky decomposition of a symmetric/Hermitian positive-definite matrix.
    #
    # Converts the matrix to `GeneralMatrix` and performs in-place decomposition.
    #
    # Arguments:
    #   - lower (Bool) : If `true`, returns lower triangular factor (L). Otherwise, upper (U). Default: `false`.
    #   - dont_clean (Bool) : If `true`, does not zero out the opposite triangle. Default: `false`.
    #
    # Returns:
    #   - GeneralMatrix(T) : The lower (L) or upper (U) triangular factor such that A = L*L' or A = U'*U.
    #
    # Side Effect:
    #   - Marks the original matrix as `PositiveDefinite` if decomposition succeeds.
    def cholesky(*, lower = false, dont_clean = false)
      to_general.cholesky!(lower: lower, dont_clean: dont_clean).tap do |m|
        # if cholesky didn't raised, original matrix is positive definite
        self.assume! MatrixFlags::PositiveDefinite
      end
    end
  end

  class GeneralMatrix(T) < Matrix(T)
    #  Performs in-place Cholesky decomposition.
    #
    # Arguments:
    #   - lower (Bool) : If `true`, computes lower triangular factor (L). Otherwise, upper (U). Default: `false`.
    #   - dont_clean (Bool) : If `true`, leaves the unused triangle unmodified. Default: `false`.
    #
    # Returns:
    #   - self : On successful decomposition, `self` becomes the triangular factor.
    #
    # Raises:
    #   - ArgumentError : If matrix is not square.
    #
    # LAPACK Routine:
    #   - Uses `potrf` (real/complex positive-definite factorization).
    def cholesky!(*, lower = false, dont_clean = false)
      raise ArgumentError.new("Matrix must be square for cholesky decomposition") unless square?
      char = lower ? 'L' : 'U'
      lapack(potrf, char.ord.to_u8, nrows, self, nrows)
      if lower
        if dont_clean
          self.flags = MatrixFlags::LowerTriangular
        else
          tril!
        end
      else
        if dont_clean
          self.flags = MatrixFlags::UpperTriangular
        else
          triu!
        end
      end
      self
    end

    # Solves A * X = B given Cholesky factorization of A (A = L*L' or U'*U).
    #
    # Arguments:
    #   - b (GeneralMatrix(T)) : Right-hand side matrix.
    #   - overwrite_b (Bool) : Allows overwriting `b` with solution. Default: `false`.
    #
    # Returns:
    #   - GeneralMatrix(T) : Solution matrix X.
    #
    # Raises:
    #   - ArgumentError : If `nrows` mismatch, `self` is not square, or not triangular.
    #
    # LAPACK Routine:
    #   - Uses `potrs` (positive-definite solve).
    def cho_solve(b : self, *, overwrite_b = false)
      raise ArgumentError.new("nrows of a and b must match") unless nrows == b.nrows
      raise ArgumentError.new("a must be square") unless square?
      raise ArgumentError.new("a must be triangular") unless flags.triangular?
      x = overwrite_b ? b : b.clone
      n = nrows
      lapack(potrs, uplo, n, b.nrows, self, n, x, b.nrows)
      b.clear_flags
      x
    end
  end
end
