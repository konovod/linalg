require "../matrix/*"
require "./lapack_helper"

# TODO - inline docs

module LA
  def self.cho_solve(a, b, *, overwrite_b = false)
    a.cho_solve(b, overwrite_b: overwrite_b)
  end

  abstract class Matrix(T)
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

    def cholesky(*, lower = false, dont_clean = false)
      clone.cholesky!(lower: lower, dont_clean: dont_clean).tap do |m|
        # if cholesky didn't raised, original matrix is positive definite
        self.assume! MatrixFlags::PositiveDefinite
      end
    end

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
