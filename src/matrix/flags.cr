module LA
  # MatrixFlags represent properties of matrix (symmetric, positive definite etc)
  # They are used to automatically select faster algorithms from LAPACK.
  # Flags are partially enforced by runtime checks, with the possibility of user override.
  # Examples:
  # ```
  # a.detect(MatrixFlags::Symmetric) # will perform a check and if matrix is symmetric, corresponding flag will be set
  # a.detect # perform check for all flags
  # a.assume!(MatrixFlags::Symmetric) # we know that matrix is symmetric, so we use `assume!` to set flag without check
  # a.transpose # will return clone of matrix, as for symmetrix matrices `a == a.transpose`
  # Simple operations correctly update flags:
  # (a + Mat.diag(*a.size)).flags.symmetric # => True
  # But direct access reset flags
  # a[1,1] = 1
  # a.flags.symmetric? # => False
  # ```
  #
  # Note that simple flag access doesn't perform check, so `a.flags.symmetric?` can be false for symmetric matrix, unless you called `detect` before.
  @[Flags]
  enum MatrixFlags
    # This flag shows that matrix is symmetric (equal to its transpose)
    Symmetric
    # This flag shows that matrix is hermitian (equal to its conjtranspose))
    Hermitian
    # This flag shows that matrix is [positive definite](https://en.wikipedia.org/wiki/Definite_matrix)
    PositiveDefinite
    # This flag shows that matrix is orthogonal (its columns and rows are orthonormal vectors)
    Orthogonal
    # This flag shows that matrix is upper triangular (all the entries below the main diagonal are zero)
    UpperTriangular
    # This flag shows that matrix is lower triangular (all the entries above the main diagonal are zero)
    LowerTriangular
    # Combination of `UpperTriangular | LowerTriangular` for internal use
    Triangular = UpperTriangular | LowerTriangular

    # TODO?
    # Hessenberg
    # Band
    # Diagonal
    # Bidiagonal
    # Tridiagonal

    # returns true if matrix is either UpperTriangular or LowerTriangular
    def triangular?
      self.upper_triangular? || self.lower_triangular?
    end

    # returns true if matrix is diagonal (both UpperTriangular and LowerTriangular)
    def diagonal?
      self.upper_triangular? && self.lower_triangular?
    end

    # :nodoc:
    def sum(f2 : MatrixFlags)
      self & f2 & (MatrixFlags::UpperTriangular |
                   MatrixFlags::LowerTriangular |
                   MatrixFlags::Symmetric |
                   MatrixFlags::Hermitian)
    end

    def add(f2 : MatrixFlags, alpha, beta)
      self.scale(alpha.is_a?(Complex) && alpha.imag != 0).sum(f2.scale(beta.is_a?(Complex) && beta.imag != 0))
    end

    # :nodoc:
    def mult(f2 : MatrixFlags)
      self & f2 & (MatrixFlags::UpperTriangular |
                   MatrixFlags::LowerTriangular |
                   MatrixFlags::Symmetric)
    end

    # :nodoc:
    def transpose
      result = self
      if triangular?
        wasup, waslo = upper_triangular?, lower_triangular?
        if wasup
          result |= MatrixFlags::LowerTriangular
        else
          result &= ~MatrixFlags::LowerTriangular
        end
        if waslo
          result |= MatrixFlags::UpperTriangular
        else
          result &= ~MatrixFlags::UpperTriangular
        end
      end
      result
    end

    # :nodoc:
    def scale(complex : Bool)
      r = self & ~MatrixFlags::Orthogonal
      complex ? r & ~MatrixFlags::Hermitian : r
    end

    # :nodoc:
    def self.for_diag(square : Bool)
      if square
        Symmetric | UpperTriangular | LowerTriangular | Hermitian
      else
        UpperTriangular | LowerTriangular
      end
    end

    # :nodoc:
    def real
      result = self & (MatrixFlags::Symmetric | MatrixFlags::Hermitian | MatrixFlags::Triangular)
      result |= MatrixFlags::Symmetric if self.hermitian?
      result
    end

    # :nodoc:
    def imag
      self & (MatrixFlags::Symmetric | MatrixFlags::Triangular)
    end

    # :nodoc:
    def tril(diagonal : Bool, square : Bool)
      if diagonal
        if self.upper_triangular?
          MatrixFlags.for_diag(square)
        else
          LowerTriangular
        end
      else
        self & (MatrixFlags::UpperTriangular | MatrixFlags::LowerTriangular)
      end
    end

    # :nodoc:
    def triu(diagonal : Bool, square : Bool)
      if diagonal
        if self.lower_triangular?
          MatrixFlags.for_diag(square)
        else
          UpperTriangular
        end
      else
        self & (MatrixFlags::UpperTriangular | MatrixFlags::LowerTriangular)
      end
    end
  end
end
