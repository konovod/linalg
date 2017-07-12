require "../matrix/*"
require "./libLAPACKE"

module Linalg
  module Matrix(T)
    macro blas(storage, name, *args)
      {% if T == Float32
           typ = :s.id
         elsif T == Float64
           typ = :d.id
         elsif T == Complex
           typ = :z.id
         end %}
       LibCBLAS.{{typ}}{{storage}}{{name}}(LibCBLAS::ROW_MAJOR, {{*args}})
    end

    private macro blas_const(x)
      {% if T == Complex %}
        pointerof({{x}}).as(LibCBLAS::ComplexDouble*)
      {% else %}
        {{x}}
      {% end %}
    end

    # performs c = alpha*a*b + beta*c (BLAS routines gemm/symm/hemm/trmm)
    def inc_mult(a, b : Matrix(T), *, alpha = 1.0, beta = 1.0)
      if a.ncolumns != b.nrows || a.nrows != nrows || b.ncolumns != ncolumns
        raise ArgumentError.new("matrix size mismatch")
      end
      aa = a.is_a?(GeneralMatrix(T)) ? a : a.clone
      bb = b.is_a?(GeneralMatrix(T)) ? b : b.clone
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
        blas(he, mm, side, up, a.nrows, b.ncolumns,
          blas_const(calpha),
          aa, a.ncolumns,
          bb, b.ncolumns,
          blas_const(cbeta),
          self, self.ncolumns)
        {% else %}
        raise "" #to prevent type inference of nil
        {% end %}
      elsif a.flags.symmetric? || b.flags.symmetric?
        side = a.flags.symmetric? ? LibCBLAS::CblasSide::CblasLeft : LibCBLAS::CblasSide::CblasRight
        if side == LibCBLAS::CblasSide::CblasRight
          aa, bb = bb, aa
        end
        up = LibCBLAS::CblasUplo::CblasUpper
        blas(sy, mm, side, up, a.nrows, b.ncolumns,
          blas_const(calpha),
          aa, a.ncolumns,
          bb, b.ncolumns,
          blas_const(cbeta),
          self, self.ncolumns)
      else
        blas(ge, mm, no, no, a.nrows, b.ncolumns, a.ncolumns,
          blas_const(calpha),
          aa, a.ncolumns,
          bb, b.ncolumns,
          blas_const(cbeta),
          self, self.ncolumns)
      end
    end

    def tr_mult!(a : Matrix(T), *, alpha = 1.0, left = false)
      raise "GeneralMatrix required for tr_mult!" unless self.is_a?(GeneralMatrix(T))
      raise ArgumentError.new("matrix size should match") if ncolumns != a.nrows
      raise "tr_mult! require square triangular matrix" unless a.square? && a.flags.triangular?
      aa = a.is_a?(GeneralMatrix(T)) ? a : a.clone
      no = LibCBLAS::CblasTranspose::CblasNoTrans
      calpha = T.new(alpha)
      auplo = aa.flags.upper_triangular? ? LibCBLAS::CblasUplo::CblasUpper : LibCBLAS::CblasUplo::CblasLower
      side = left ? LibCBLAS::CblasSide::CblasLeft : LibCBLAS::CblasSide::CblasRight
      blas(tr, mm, side, auplo, no, LibCBLAS::CblasDiag::CblasNonUnit,
        self.nrows, self.ncolumns, blas_const(calpha),
        aa, aa.ncolumns,
        self, ncolumns)
      self
    end

    def *(m : Matrix(T))
      if ncolumns != m.nrows
        raise ArgumentError.new("matrix size should match ([#{nrows}x#{ncolumns}] * [#{m.nrows}x#{m.ncolumns}]")
      end
      if (square? && flags.triangular?) || (m.square? && m.flags.triangular?)
        if m.square? && m.flags.triangular?
          result = self.clone
          result.tr_mult!(m, left: false)
          result
        else
          result = m.clone
          result.tr_mult!(self, left: true)
          result
        end
      else
        result = Matrix(T).zeros(nrows, m.ncolumns)
        result.inc_mult(self, m)
        result.tap { |r| r.flags = self.flags.mult(m.flags) }
      end
    end
  end
end
