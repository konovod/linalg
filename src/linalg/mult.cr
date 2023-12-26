require "../matrix/*"

module LA
  class GeneralMatrix(T) < Matrix(T)
    # :nodoc:
    macro blas(storage, name, *args)
      {% if T == Float32
           typ = :s.id
         elsif T == Float64
           typ = :d.id
         elsif T == Complex
           typ = :z.id
         end %}
       LibCBLAS.{{typ}}{{storage}}{{name}}(LibCBLAS::COL_MAJOR, {{*args}})
    end

    private macro blas_const(x)
      {% if T == Complex %}
        pointerof({{x}}).as(LibCBLAS::ComplexDouble*)
      {% else %}
        {{x}}
      {% end %}
    end

    # performs c = alpha*a*b + beta*c (BLAS routines `gemm`/`symm`/`hemm`)
    def add_mult(a, b : Matrix(T), *, alpha = 1.0, beta = 1.0)
      if a.ncolumns != b.nrows || a.nrows != nrows || b.ncolumns != ncolumns
        raise ArgumentError.new("matrix size mismatch")
      end
      aa = a.is_a?(GeneralMatrix(T)) ? a : a.to_general
      bb = b.is_a?(GeneralMatrix(T)) ? b : b.to_general
      no = LibCBLAS::CblasTranspose::CblasNoTrans
      calpha = T.new(alpha)
      cbeta = T.new(beta)
      if {{T == Complex}} && (a.flags.hermitian? || b.flags.hermitian?)
        {% if T == Complex %}
          side = a.flags.hermitian? ? LibCBLAS::CblasSide::CblasLeft : LibCBLAS::CblasSide::CblasRight
          if side == LibCBLAS::CblasSide::CblasRight
            aa, bb = bb, aa
          end
          up = LibCBLAS::CblasUplo::CblasUpper
          blas(he, mm, side, up, a.ncolumns, b.nrows,
            blas_const(calpha),
            aa, a.nrows,
            bb, b.nrows,
            blas_const(cbeta),
            self, self.nrows)
        {% else %}
          raise "error" # to prevent type inference of nil
        {% end %}
      elsif a.flags.symmetric? || b.flags.symmetric?
        side = a.flags.symmetric? ? LibCBLAS::CblasSide::CblasLeft : LibCBLAS::CblasSide::CblasRight
        if side == LibCBLAS::CblasSide::CblasRight
          aa, bb = bb, aa
        end
        up = LibCBLAS::CblasUplo::CblasUpper
        blas(sy, mm, side, up, a.ncolumns, b.nrows,
          blas_const(calpha),
          aa, a.nrows,
          bb, b.nrows,
          blas_const(cbeta),
          self, self.nrows)
      else
        blas(ge, mm, no, no, a.nrows, b.ncolumns, a.ncolumns,
          blas_const(calpha),
          aa, a.nrows,
          bb, b.nrows,
          blas_const(cbeta),
          self, self.nrows)
      end
    end

    # performs b = alpha*a*b or b = alpha*b*a (BLAS routine `trmm`)
    #
    # `a` must be a square triangular `GeneralMatrix(T)`
    #
    # if `left` is true, `alpha*a*b` is calculated, otherwise `alpha*b*a`
    def tr_mult!(a : Matrix(T), *, alpha = 1.0, left = false)
      raise ArgumentError.new("matrix size should match") if ncolumns != a.nrows
      raise "tr_mult! require square triangular matrix" unless a.square? && a.flags.triangular?
      aa = a.is_a?(GeneralMatrix(T)) ? a : a.to_general
      no = LibCBLAS::CblasTranspose::CblasNoTrans
      calpha = T.new(alpha)
      auplo = aa.flags.upper_triangular? ? LibCBLAS::CblasUplo::CblasUpper : LibCBLAS::CblasUplo::CblasLower
      side = left ? LibCBLAS::CblasSide::CblasLeft : LibCBLAS::CblasSide::CblasRight
      blas(tr, mm, side, auplo, no, LibCBLAS::CblasDiag::CblasNonUnit,
        self.nrows, self.ncolumns, blas_const(calpha),
        aa, aa.nrows,
        self, nrows)
      self
    end
  end

  abstract class Matrix(T)
    # Matrix product to given m
    #
    # Raises ArgumentError if inner dimensions do not match
    #
    # This method automatically calls optimal function depending on `MatrixFlags`.
    #
    # If one of the matrix is square and triangular - `trmm` is called
    #
    # If one of the matrix is symmetric\hermitian - `symm`/`hemm` is called
    #
    # Otherwise - `gemm` is called
    def *(m : Matrix(T))
      if ncolumns != m.nrows
        raise ArgumentError.new("matrix size should match (#{shape_str} * #{m.shape_str}")
      end
      if (square? && flags.triangular?) || (m.square? && m.flags.triangular?)
        if m.square? && m.flags.triangular?
          result = self.to_general
          result.tr_mult!(m, left: false)
        else
          result = m.to_general
          result.tr_mult!(self, left: true)
        end
      else
        result = GeneralMatrix(T).zeros(nrows, m.ncolumns)
        result.add_mult(self, m)
      end
      result.tap { |r| r.flags = self.flags.mult(m.flags) }
    end
  end
end
