require "../matrix/*"
require "./libLAPACKE"

module Linalg
  def self.cho_solve(a, b, *, overwrite_b = false)
    a.cho_solve(b, overwrite_b: overwrite_b)
  end

  module Matrix(T)
    def cholesky!(*, lower = false, dont_clean = false)
      raise ArgumentError.new("Matrix must be square for cholesky decomposition") unless square?
      raise ArgumentError.new("Matrix must be positive definite for cholesky decomposition") unless flags.positive_definite?
      char = lower ? 'L' : 'U'
      lapack(po, trf, char.ord, rows, self, rows)
      if lower
        if dont_clean
          self.flags = MatrixFlags::Triangular | MatrixFlags::Lower
        else
          tril!
        end
      else
        if dont_clean
          self.flags = MatrixFlags::Triangular
        else
          triu!
        end
      end
      self
    end

    def cholesky(*, lower = false, dont_clean = false)
      clone.cholesky!(lower: lower, dont_clean: dont_clean)
    end

    def cho_solve(b : self, *, overwrite_b = false)
      raise ArgumentError.new("number of rows in a and b must match") unless rows == b.rows
      raise ArgumentError.new("a must be square") unless square?
      raise ArgumentError.new("a must be triangular") unless flags.triangular?
      x = overwrite_b ? b : b.clone
      n = rows
      lapack(po, trs, uplo, n, b.columns, self, n, x, b.columns)
      x
    end
  end
end
