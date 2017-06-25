module Linalg
  module Matrix(T)
    private def qr_initial(a, pivoting)
      m = a.rows
      n = a.columns
      tau = GeneralMatrix(T).new(1, n)
      if pivoting
        jpvt = Array(Int32).new(n, 0)
        lapack(ge, qp3, m, n, a, a.columns, jpvt, tau)
        a.clear_flags
      else
        jpvt = Array(Int32).new(0)
        lapack(ge, qrf, m, n, a, a.columns, tau)
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
      {a, tau, pvt}
    end

    def qr(*, overwrite_a = false, pivoting = false)
      a = overwrite_a ? self : clone
      tau, pvt = qr_initial(a, pivoting)
      r = a.clone
      r.triu!
      {% if T == Complex %}
        m = a.rows
        n = a.columns
        k = {m,n}.min
        lapack(un, gqr, m, n, k, a, a.columns, tau)
      {% else %}
        m = a.rows
        n = a.columns
        k = {m,n}.min
        lapack(or, gqr, m, n, k, a, a.columns, tau)
      {% end %}
      a.assume! MatrixFlags::Orthogonal
      {a, r, pvt}
    end
  end
end
