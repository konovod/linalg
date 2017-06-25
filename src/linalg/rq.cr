module Linalg
  module Matrix(T)
    private def rq_initial(a)
      m = a.rows
      n = a.columns
      tau = GeneralMatrix(T).new(1, n)
      lapack(ge, rqf, m, n, a, a.columns, tau)
      a.clear_flags
      tau
    end

    def rq_r(*, overwrite_a = false)
      a = overwrite_a ? self : clone
      tau = rq_initial(a)
      a.triu!
      a
    end

    def rq(*, overwrite_a = false)
      a = overwrite_a ? self : clone
      tau = rq_initial(a)
      r = a.clone
      r.triu!
      {% if T == Complex %}
        m = a.rows
        n = a.columns
        k = {m,n}.min
        lapack(un, grq, m, n, k, a, a.columns, tau)
      {% else %}
        m = a.rows
        n = a.columns
        k = {m,n}.min
        lapack(or, grq, m, n, k, a, a.columns, tau)
      {% end %}
      a.assume! MatrixFlags::Orthogonal
      {r, a}
    end
  end
end
