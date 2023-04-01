# TODO - inline docs

module LA
  abstract class Matrix(T)
    private def qr_initial(a, pivoting)
      m = a.nrows
      n = a.ncolumns
      tau = GeneralMatrix(T).new(1, n)
      if pivoting
        jpvt = Array(Int32).new(n, 0)
        lapack(geqp3, m, n, a, a.nrows, jpvt, tau, worksize: [2*n])
        a.clear_flags
      else
        jpvt = Array(Int32).new(0)
        lapack(geqrf, m, n, a, a.nrows, tau)
        a.clear_flags
      end
      {tau, jpvt}
    end

    def qr_raw(*, overwrite_a = false, pivoting = false)
      a = overwrite_a ? self : clone
      tau, pvt = qr_initial(a, pivoting)
      {a, tau, pvt}
    end

    def qr_r(*, overwrite_a = false, pivoting = false)
      a = overwrite_a ? self : clone
      tau, pvt = qr_initial(a, pivoting)
      a.triu!
      {a, pvt}
    end

    def qr(*, overwrite_a = false, pivoting = false)
      a = overwrite_a ? self : clone
      tau, pvt = qr_initial(a, pivoting)
      r = a.clone
      r.triu!
      m = a.nrows
      n = a.ncolumns
      k = {m, n}.min
      lapack(orgqr, m, n, k, a, a.nrows, tau)
      a.assume! MatrixFlags::Orthogonal
      {a, r, pvt}
    end
  end
end
