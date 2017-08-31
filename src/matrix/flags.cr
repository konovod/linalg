module LA
  @[Flags]
  enum MatrixFlags
    Symmetric
    Hermitian
    PositiveDefinite
    # Hessenberg
    # Band
    # Diagonal
    # Bidiagonal
    # Tridiagonal
    Orthogonal
    UpperTriangular
    LowerTriangular
    Triangular      = UpperTriangular | LowerTriangular

    def triangular?
      self.upper_triangular? || self.lower_triangular?
    end

    def diagonal?
      self.upper_triangular? && self.lower_triangular?
    end

    def sum(f2 : MatrixFlags)
      self & f2 & (MatrixFlags::UpperTriangular |
                   MatrixFlags::LowerTriangular |
                   MatrixFlags::Symmetric |
                   MatrixFlags::Hermitian)
    end

    def mult(f2 : MatrixFlags)
      self & f2 & (MatrixFlags::UpperTriangular |
                   MatrixFlags::LowerTriangular |
                   MatrixFlags::Symmetric)
    end

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

    def scale(complex : Bool)
      r = self & ~MatrixFlags::Orthogonal
      complex ? r & ~MatrixFlags::Hermitian : r
    end

    def self.for_diag(square : Bool)
      if square
        Symmetric | UpperTriangular | LowerTriangular | Hermitian
      else
        UpperTriangular | LowerTriangular
      end
    end

    def real
      result = self & (MatrixFlags::Symmetric | MatrixFlags::Hermitian | MatrixFlags::Triangular)
      result |= MatrixFlags::Symmetric if self.hermitian?
      result
    end

    def imag
      self & (MatrixFlags::Symmetric | MatrixFlags::Triangular)
    end

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
