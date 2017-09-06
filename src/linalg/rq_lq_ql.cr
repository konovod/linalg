private macro decomposition(letters, trix)

private def {{letters}}_initial(a)
  m = a.nrows
  n = a.ncolumns
  tau = GeneralMatrix(T).new(1, n)
  lapack(ge, {{letters}}f, m, n, a, a.ncolumns, tau)
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
  m = a.nrows
  n = a.ncolumns
  k = {m,n}.min
  lapack(or, g{{letters}}, m, n, k, a, a.ncolumns, tau)
  a.assume! MatrixFlags::Orthogonal
  {% if letters.id.chars[0] == 'q' %}
    {a, r}
  {% else %}
    {r, a}
  {% end %}
end

end

module LA
  abstract class Matrix(T)
    decomposition(rq, triu)
    decomposition(lq, tril)
    decomposition(ql, tril)
  end
end
