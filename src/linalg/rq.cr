private macro decomposition(letters, trix)

private def {{letters}}_initial(a)
  m = a.rows
  n = a.columns
  tau = GeneralMatrix(T).new(1, n)
  lapack(ge, {{letters}}f, m, n, a, a.columns, tau)
  a.clear_flags
  tau
end

def {{letters}}_r(*, overwrite_a = false)
  a = overwrite_a ? self : clone
  tau = {{letters}}_initial(a)
  a.{{trix}}!
  a
end

def {{letters}}(*, overwrite_a = false)
  a = overwrite_a ? self : clone
  tau = {{letters}}_initial(a)
  r = a.clone
  r.{{trix}}!
  m = a.rows
  n = a.columns
  k = {m,n}.min
  lapack(or, g{{letters}}, m, n, k, a, a.columns, tau)
  a.assume! MatrixFlags::Orthogonal
  {r, a}
end

end

module Linalg
  module Matrix(T)
    decomposition(rq, triu)
    decomposition(lq, tril)
    decomposition(ql, tril)
  end
end
