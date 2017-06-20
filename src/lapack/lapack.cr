module LAPACK
  # matrix inversion using xxxtrf / xxxtri

  #   macro call_f(name, args)
  #   LibLAPACKE.dgetrf(LibLAPACKE::ROW_MAJOR, 3, 3, matrix2, 3, ipiv)
  # end

  class Matrix(T)
    def inv
      raise "can't invert nonsquare matrix" unless @rows == @columns
      result = clone
      n = @rows
      ipiv = Slice(Int32).new(n)
      LibLAPACKE.dgetrf(LibLAPACKE::ROW_MAJOR, n, n, result, n, ipiv)
      LibLAPACKE.dgetri(LibLAPACKE::ROW_MAJOR, n, result, n, ipiv)
      result
    end
  end
end
